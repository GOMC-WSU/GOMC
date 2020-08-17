#ifdef GOMC_CUDA
#include "TransformParticlesCUDAKernel.cuh"

#define MIN_FORCE 1E-12
#define MAX_FORCE 30

__device__ inline double randomGPU(unsigned int counter, unsigned int step, unsigned int seed) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  RNG::ctr_type r = philox4x32(c, k);
  return (double)r[0] / UINT_MAX;
}

__device__ inline double WrapPBC(double &v, double ax) {
  if(v >= ax)
    v -= ax;
  else if(v < 0)
    v += ax;
  return v;
}

__device__ inline double UnwrapPBC(double &v, double ref, double ax, double halfax) {
  if(abs(ref - v) > halfax) {
    if(ref < halfax)
      v -= ax;
    else
      v += ax;
  }
  return v;
}

__device__ inline void ApplyRotation(double &x, double &y, double &z,
                                     double comx, double comy, double comz,
                                     double rotx, double roty, double rotz,
                                     double axx, double axy, double axz, int atomNumber)
{
  double rotLen = sqrt(rotx * rotx + roty * roty + rotz * rotz);
  double axisx = rotx * (1.0 / rotLen);
  double axisy = roty * (1.0 / rotLen);
  double axisz = rotz * (1.0 / rotLen);
  double matrix[3][3], cross[3][3], tensor[3][3];

  // build cross
  cross[0][0] = 0.0;
  cross[0][1] = -axisz;
  cross[0][2] = axisy;

  cross[1][0] = axisz;
  cross[1][1] = 0.0;
  cross[1][2] = -axisx;

  cross[2][0] = -axisy;
  cross[2][1] = axisx;
  cross[2][2] = 0.0;

  // build tensor
  for(int i=0; i<3; i++) {
    tensor[0][i] = axisx;
    tensor[1][i] = axisy;
    tensor[2][i] = axisz;
  }
  for(int i=0; i<3; i++) {
    tensor[i][0] *= axisx;
    tensor[i][1] *= axisy;
    tensor[i][2] *= axisz;
  }

  // build matrix
  double c = cos(rotLen);
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      matrix[i][j] = 0.0;
    }
    matrix[i][i] = c;
  }
  double s = sin(rotLen);
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      matrix[i][j] += s * cross[i][j] + (1-c) * tensor[i][j];
    }
  }

  // unwrap molecule
  UnwrapPBC(x, comx, axx, axx/2.0);
  UnwrapPBC(y, comy, axy, axy/2.0);
  UnwrapPBC(z, comz, axz, axz/2.0);

  // move particle to zero
  x -= comx;
  y -= comy;
  z -= comz;

  // rotate
  double newx = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z;
  double newy = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z;
  double newz = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z;

  x = newx;
  y = newy;
  z = newz;

  // move back to com
  x += comx;
  y += comy;
  z += comz;

  // wrap again
  WrapPBC(x, axx);
  WrapPBC(y, axy);
  WrapPBC(z, axz);
}

void CallTranslateParticlesGPU(VariablesCUDA *vars,
                               std::vector<uint> &moleculeIndex,
                               double t_max,
                               double *mForcex,
                               double *mForcey,
                               double *mForcez,
                               unsigned int step,
                               unsigned int seed,
                               std::vector<int> particleMol,
                               int atomCount,
                               int molCount,
                               double xAxes,
                               double yAxes,
                               double zAxes,
                               XYZArray &newMolPos,
                               XYZArray &newCOMs,
                               double lambdaBETA)
{
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)(atomCount / threadsPerBlock) + 1;
  int *gpu_particleMol;

  cudaMalloc((void**) &gpu_particleMol, particleMol.size() * sizeof(int));

  cudaMemcpy(vars->gpu_mForcex, mForcex, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcey, mForcey, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcez, mForcez, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double), cudaMemcpyHostToDevice);

  checkLastErrorCUDA(__FILE__, __LINE__);
  TranslateParticlesKernel<<<blocksPerGrid, threadsPerBlock>>>(molCount,
                                                               t_max,
                                                               vars->gpu_mForcex,
                                                               vars->gpu_mForcey,
                                                               vars->gpu_mForcez,
                                                               step,
                                                               seed,
                                                               vars->gpu_x,
                                                               vars->gpu_y,
                                                               vars->gpu_z,
                                                               gpu_particleMol,
                                                               atomCount,
                                                               xAxes,
                                                               yAxes,
                                                               zAxes,
                                                               vars->gpu_comx,
                                                               vars->gpu_comy,
                                                               vars->gpu_comz,
                                                               lambdaBETA);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
  
  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.x, vars->gpu_comx, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.y, vars->gpu_comy, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.z, vars->gpu_comz, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(gpu_particleMol);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void CallRotateParticlesGPU(VariablesCUDA *vars,
                            std::vector<uint> &moleculeIndex,
                            double r_max,
                            double *mTorquex,
                            double *mTorquey,
                            double *mTorquez,
                            unsigned int step,
                            unsigned int seed,
                            std::vector<int> particleMol,
                            int atomCount,
                            int molCount,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            XYZArray &newMolPos,
                            XYZArray &newCOMs)
{
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)(atomCount / threadsPerBlock) + 1;
  int *gpu_particleMol;

  cudaMalloc((void**) &gpu_particleMol, particleMol.size() * sizeof(int));

  cudaMemcpy(vars->gpu_mTorquex, mTorquex, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquey, mTorquey, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquez, mTorquez, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double), cudaMemcpyHostToDevice);

  RotateParticlesKernel<<<blocksPerGrid, threadsPerBlock>>>(molCount,
                                                            r_max,
                                                            vars->gpu_mTorquex,
                                                            vars->gpu_mTorquey,
                                                            vars->gpu_mTorquez,
                                                            step,
                                                            seed,
                                                            vars->gpu_x,
                                                            vars->gpu_y,
                                                            vars->gpu_z,
                                                            gpu_particleMol,
                                                            atomCount,
                                                            xAxes,
                                                            yAxes,
                                                            zAxes,
                                                            vars->gpu_comx,
                                                            vars->gpu_comy,
                                                            vars->gpu_comz);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);
  
  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(gpu_particleMol);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

__global__ void TranslateParticlesKernel(unsigned int numberOfMolecules,
                                         double t_max,
                                         double *molForcex,
                                         double *molForcey,
                                         double *molForcez,
                                         unsigned int step,
                                         unsigned int seed,
                                         double *gpu_x,
                                         double *gpu_y,
                                         double *gpu_z,
                                         int *gpu_particleMol,
                                         int atomCount,
                                         double xAxes,
                                         double yAxes,
                                         double zAxes,
                                         double *gpu_comx,
                                         double *gpu_comy,
                                         double *gpu_comz,
                                         double lambdaBETA)
{
  int atomNumber = blockIdx.x * blockDim.x + threadIdx.x;
  if(atomNumber >= atomCount) return;

  int molIndex = gpu_particleMol[atomNumber];
  bool updateCOM = atomNumber == 0 || (gpu_particleMol[atomNumber] != gpu_particleMol[atomNumber-1]);

  // This section calculates the amount of shift
  double lbfx = molForcex[molIndex] * lambdaBETA;
  double lbfy = molForcey[molIndex] * lambdaBETA;
  double lbfz = molForcez[molIndex] * lambdaBETA;
  double lbmaxx = lbfx * t_max;
  double lbmaxy = lbfy * t_max;
  double lbmaxz = lbfz * t_max;

  double shiftx, shifty, shiftz;

  if(abs(lbmaxx) > MIN_FORCE && abs(lbmaxx) < MAX_FORCE) {
    shiftx = log(exp(-1.0 * lbmaxx) + 2 * randomGPU(molIndex * 3, step, seed) * sinh(lbmaxx)) / lbfx;
  } else {
    double rr = randomGPU(molIndex * 3, step, seed) * 2.0 - 1.0;
    shiftx = t_max * rr;
  }

  if(abs(lbmaxy) > MIN_FORCE && abs(lbmaxy) < MAX_FORCE) {
    shifty = log(exp(-1.0 * lbmaxy) + 2 * randomGPU(molIndex * 3 + 1, step, seed) * sinh(lbmaxy)) / lbfy;
  } else {
    double rr = randomGPU(molIndex * 3 + 1, step, seed) * 2.0 - 1.0;
    shifty = t_max * rr;
  }

  if(abs(lbmaxz) > MIN_FORCE && abs(lbmaxz) < MAX_FORCE) {
    shiftz = log(exp(-1.0 * lbmaxz) + 2 * randomGPU(molIndex * 3 + 2, step, seed) * sinh(lbmaxz)) / lbfz;
  } else {
    double rr = randomGPU(molIndex * 3 + 2, step, seed) * 2.0 - 1.0;
    shiftz = t_max * rr;
  }

  // perform the shift on the coordinates
  gpu_x[atomNumber] += shiftx;
  gpu_y[atomNumber] += shifty;
  gpu_z[atomNumber] += shiftz;

  // rewrapping
  WrapPBC(gpu_x[atomNumber], xAxes);
  WrapPBC(gpu_y[atomNumber], yAxes);
  WrapPBC(gpu_z[atomNumber], zAxes);

  if(updateCOM) {
    gpu_comx[molIndex] += shiftx;
    gpu_comy[molIndex] += shifty;
    gpu_comz[molIndex] += shiftz;
  }
}

__global__ void RotateParticlesKernel(unsigned int numberOfMolecules,
                                      double r_max,
                                      double *molTorquex,
                                      double *molTorquey,
                                      double *molTorquez,
                                      unsigned int step,
                                      unsigned int seed,
                                      double *gpu_x,
                                      double *gpu_y,
                                      double *gpu_z,
                                      int *gpu_particleMol,
                                      int atomCount,
                                      double xAxes,
                                      double yAxes,
                                      double zAxes,
                                      double *gpu_comx,
                                      double *gpu_comy,
                                      double *gpu_comz)
{
  int atomNumber = blockIdx.x * blockDim.x + threadIdx.x;
  if(atomNumber >= atomCount) return;

  int molIndex = gpu_particleMol[atomNumber];

  // This section calculates the amount of shift
  double lbtx = molTorquex[molIndex];
  double lbty = molTorquey[molIndex];
  double lbtz = molTorquez[molIndex];
  double lbmaxx = lbtx * r_max;
  double lbmaxy = lbty * r_max;
  double lbmaxz = lbtz * r_max;

  double rotx, roty, rotz;

  if(abs(lbmaxx) > MIN_FORCE && abs(lbmaxx) < MAX_FORCE) {
    rotx = log(exp(-1.0 * lbmaxx) + 2 * randomGPU(molIndex * 3, step, seed) * sinh(lbmaxx)) / lbtx;
  } else {
    double rr = randomGPU(molIndex * 3, step, seed) * 2.0 - 1.0;
    rotx = r_max * rr;
  }

  if(abs(lbmaxy) > MIN_FORCE && abs(lbmaxy) < MAX_FORCE) {
    roty = log(exp(-1.0 * lbmaxy) + 2 * randomGPU(molIndex * 3 + 1, step, seed) * sinh(lbmaxy)) / lbty;
  } else {
    double rr = randomGPU(molIndex * 3 + 1, step, seed) * 2.0 - 1.0;
    roty = r_max * rr;
  }

  if(abs(lbmaxz) > MIN_FORCE && abs(lbmaxz) < MAX_FORCE) {
    rotz = log(exp(-1.0 * lbmaxz) + 2 * randomGPU(molIndex * 3 + 2, step, seed) * sinh(lbmaxz)) / lbtz;
  } else {
    double rr = randomGPU(molIndex * 3 + 2, step, seed) * 2.0 - 1.0;
    rotz = r_max * rr;
  }

  if(atomNumber == 0) {
    print_tuple("GPU", r_max, r_max, r_max);
  }

  // perform the rot on the coordinates
  ApplyRotation(gpu_x[atomNumber], gpu_y[atomNumber], gpu_z[atomNumber],
                gpu_comx[molIndex], gpu_comy[molIndex], gpu_comz[molIndex],
                rotx, roty, rotz, xAxes, yAxes, zAxes, atomNumber);
}

#endif
