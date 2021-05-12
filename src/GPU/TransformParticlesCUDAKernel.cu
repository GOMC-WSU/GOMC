/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifdef GOMC_CUDA
#include "TransformParticlesCUDAKernel.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "CUDAMemoryManager.cuh"
#include "Random123/boxmuller.hpp"

#define MIN_FORCE 1E-12
#define MAX_FORCE 30

__device__ inline double randomGPU(unsigned int counter, ulong step, ulong seed)
{
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  RNG::ctr_type r = philox4x64(c, k);
  return (double)r[0] / ULONG_MAX;
}

__device__ inline double randomGaussianGPU(unsigned int counter, ulong step,
                                           ulong seed, double mean, double stdDev)
{
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  RNG::ctr_type r = philox4x64(c, k);
  double2 normal2 = r123::boxmuller(r[0], r[1]);
  double shiftedVal = mean + normal2.x * stdDev;
  return  shiftedVal;
}

__device__ inline void ApplyRotation(double &x, double &y, double &z,
                                     double comx, double comy, double comz,
                                     double rotx, double roty, double rotz,
                                     double axx, double axy, double axz)
{
  double rotLen = sqrt(rotx * rotx + roty * roty + rotz * rotz);
  double axisx = rotx * (1.0 / rotLen);
  double axisy = roty * (1.0 / rotLen);
  double axisz = rotz * (1.0 / rotLen);
  double matrix[3][3], cross[3][3], tensor[3][3];
  double halfAxx = axx * 0.5;
  double halfAxy = axy * 0.5;
  double halfAxz = axz * 0.5;

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
  for(int i = 0; i < 3; i++) {
    tensor[0][i] = axisx;
    tensor[1][i] = axisy;
    tensor[2][i] = axisz;
  }
  for(int i = 0; i < 3; i++) {
    tensor[i][0] *= axisx;
    tensor[i][1] *= axisy;
    tensor[i][2] *= axisz;
  }

  // build matrix
  double s, c;
  sincos(rotLen, &s, &c);
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      matrix[i][j] = 0.0;
    }
    matrix[i][i] = c;
  }
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      matrix[i][j] += s * cross[i][j] + (1 - c) * tensor[i][j];
    }
  }

  // unwrap molecule
  UnwrapPBC(x, comx, axx, halfAxx);
  UnwrapPBC(y, comy, axy, halfAxy);
  UnwrapPBC(z, comz, axz, halfAxz);

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
                               const std::vector<int8_t> &isMoleculeInvolved,
                               double t_max,
                               double *mForcex,
                               double *mForcey,
                               double *mForcez,
                               ulong step,
                               ulong seed,
                               const std::vector<int> &particleMol,
                               int atomCount,
                               int molCount,
                               double xAxes,
                               double yAxes,
                               double zAxes,
                               XYZArray &newMolPos,
                               XYZArray &newCOMs,
                               double lambdaBETA,
                               XYZArray &t_k,
                               XYZArray &molForceRecRef)
{
  int8_t *gpu_isMoleculeInvolved;
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)(atomCount / threadsPerBlock) + 1;
  int *gpu_particleMol;

  CUMALLOC((void **) &gpu_isMoleculeInvolved,
           isMoleculeInvolved.size() * sizeof(int8_t));
  CUMALLOC((void**) &gpu_particleMol, particleMol.size() * sizeof(int));

  cudaMemcpy(vars->gpu_mForcex, mForcex, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcey, mForcey, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcez, mForcez, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecx, molForceRecRef.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecy, molForceRecRef.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecz, molForceRecRef.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_isMoleculeInvolved, &isMoleculeInvolved[0],
             isMoleculeInvolved.size() * sizeof(int8_t), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double), cudaMemcpyHostToDevice);

  checkLastErrorCUDA(__FILE__, __LINE__);
  TranslateParticlesKernel <<< blocksPerGrid, threadsPerBlock>>>(molCount,
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
      lambdaBETA,
      vars->gpu_t_k_x,
      vars->gpu_t_k_y,
      vars->gpu_t_k_z,
      gpu_isMoleculeInvolved,
      vars->gpu_mForceRecx,
      vars->gpu_mForceRecy,
      vars->gpu_mForceRecz);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.x, vars->gpu_comx, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.y, vars->gpu_comy, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.z, vars->gpu_comz, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.x, vars->gpu_t_k_x, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.y, vars->gpu_t_k_y, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.z, vars->gpu_t_k_z, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  CUFREE(gpu_isMoleculeInvolved);
  CUFREE(gpu_particleMol);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void CallRotateParticlesGPU(VariablesCUDA *vars,
                            const std::vector<int8_t> &isMoleculeInvolved,
                            double r_max,
                            double *mTorquex,
                            double *mTorquey,
                            double *mTorquez,
                            ulong step,
                            ulong seed,
                            const std::vector<int> &particleMol,
                            int atomCount,
                            int molCount,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            XYZArray &newMolPos,
                            XYZArray &newCOMs,
                            double lambdaBETA,
                            XYZArray &r_k)
{
  int8_t *gpu_isMoleculeInvolved;
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)(atomCount / threadsPerBlock) + 1;
  int *gpu_particleMol;

  CUMALLOC((void **) &gpu_isMoleculeInvolved,
           isMoleculeInvolved.size() * sizeof(int8_t));
  CUMALLOC((void**) &gpu_particleMol, particleMol.size() * sizeof(int));

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
  cudaMemcpy(gpu_isMoleculeInvolved, &isMoleculeInvolved[0],
             isMoleculeInvolved.size() * sizeof(int8_t), cudaMemcpyHostToDevice);

  RotateParticlesKernel <<< blocksPerGrid, threadsPerBlock>>>(molCount,
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
      vars->gpu_comz,
      lambdaBETA,
      vars->gpu_r_k_x,
      vars->gpu_r_k_y,
      vars->gpu_r_k_z,
      gpu_isMoleculeInvolved);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.x, vars->gpu_r_k_x, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.y, vars->gpu_r_k_y, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.z, vars->gpu_r_k_z, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  CUFREE(gpu_isMoleculeInvolved);
  CUFREE(gpu_particleMol);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

__global__ void TranslateParticlesKernel(unsigned int numberOfMolecules,
    double t_max,
    double *molForcex,
    double *molForcey,
    double *molForcez,
    ulong step,
    ulong seed,
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
    double lambdaBETA,
    double *gpu_t_k_x,
    double *gpu_t_k_y,
    double *gpu_t_k_z,
    int8_t *gpu_isMoleculeInvolved,
    double *gpu_mForceRecx,
    double *gpu_mForceRecy,
    double *gpu_mForceRecz)
{
  int atomNumber = blockIdx.x * blockDim.x + threadIdx.x;
  if(atomNumber >= atomCount) return;

  int molIndex = gpu_particleMol[atomNumber];
  if(!gpu_isMoleculeInvolved[molIndex]) return;
  bool updateMol = atomNumber == 0 || (gpu_particleMol[atomNumber] != gpu_particleMol[atomNumber - 1]);

  // This section calculates the amount of shift
  double lbfx = (molForcex[molIndex] + gpu_mForceRecx[molIndex]) * lambdaBETA;
  double lbfy = (molForcey[molIndex] + gpu_mForceRecy[molIndex]) * lambdaBETA;
  double lbfz = (molForcez[molIndex] + gpu_mForceRecz[molIndex]) * lambdaBETA;
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

  if(updateMol) {
    gpu_comx[molIndex] += shiftx;
    gpu_comy[molIndex] += shifty;
    gpu_comz[molIndex] += shiftz;

    WrapPBC(gpu_comx[molIndex], xAxes);
    WrapPBC(gpu_comy[molIndex], yAxes);
    WrapPBC(gpu_comz[molIndex], zAxes);

    gpu_t_k_x[molIndex] = shiftx;
    gpu_t_k_y[molIndex] = shifty;
    gpu_t_k_z[molIndex] = shiftz;
  }
}

__global__ void RotateParticlesKernel(unsigned int numberOfMolecules,
                                      double r_max,
                                      double *molTorquex,
                                      double *molTorquey,
                                      double *molTorquez,
                                      ulong step,
                                      ulong seed,
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
                                      double lambdaBETA,
                                      double *gpu_r_k_x,
                                      double *gpu_r_k_y,
                                      double *gpu_r_k_z,
                                      int8_t *gpu_isMoleculeInvolved)
{
  int atomNumber = blockIdx.x * blockDim.x + threadIdx.x;
  if(atomNumber >= atomCount) return;
  int molIndex = gpu_particleMol[atomNumber];
  if(!gpu_isMoleculeInvolved[molIndex]) return;
  bool updateMol = atomNumber == 0 || (gpu_particleMol[atomNumber] != gpu_particleMol[atomNumber - 1]);

  // This section calculates the amount of rotation
  double lbtx = molTorquex[molIndex] * lambdaBETA;
  double lbty = molTorquey[molIndex] * lambdaBETA;
  double lbtz = molTorquez[molIndex] * lambdaBETA;
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

  if(updateMol) {
    gpu_r_k_x[molIndex] = rotx;
    gpu_r_k_y[molIndex] = roty;
    gpu_r_k_z[molIndex] = rotz;
  }

  // perform the rotation on the coordinates
  ApplyRotation(gpu_x[atomNumber], gpu_y[atomNumber], gpu_z[atomNumber],
                gpu_comx[molIndex], gpu_comy[molIndex], gpu_comz[molIndex],
                rotx, roty, rotz, xAxes, yAxes, zAxes);
}

// CUDA implementation of MultiParticle Brownian transformation 

void BrownianMotionRotateParticlesGPU(
  VariablesCUDA *vars,
  const std::vector<unsigned int> &moleculeInvolved,
  XYZArray &mTorque,
  XYZArray &newMolPos,
  XYZArray &newCOMs,
  XYZArray &r_k,
  const XYZ &boxAxes,
  const double BETA,
  const double r_max,
  ulong step,
  ulong seed,
  const int box,
  bool isOrthogonal,
  int *kill)
{
  int atomCount = newMolPos.Count();
  int molCount = newCOMs.Count();
  int molCountInBox = moleculeInvolved.size();
  int *gpu_moleculeInvolved;
  // Each block would handle one molecule
  int threadsPerBlock = 32;
  int blocksPerGrid = molCountInBox;

  CUMALLOC((void **) &gpu_moleculeInvolved, molCountInBox * sizeof(int));

  cudaMemcpy(vars->gpu_mTorquex, mTorque.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquey, mTorque.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquez, mTorque.z, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_moleculeInvolved, &moleculeInvolved[0], molCountInBox * sizeof(int), cudaMemcpyHostToDevice);

  double3 axis = make_double3(boxAxes.x, boxAxes.y, boxAxes.z);
  double3 halfAx = make_double3(boxAxes.x / 2.0, boxAxes.y / 2.0, boxAxes.z / 2.0);

  if(isOrthogonal) {
    BrownianMotionRotateKernel<true><<< blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_startAtomIdx,
      vars->gpu_x,
      vars->gpu_y,
      vars->gpu_z,
      vars->gpu_mTorquex,
      vars->gpu_mTorquey,
      vars->gpu_mTorquez,
      vars->gpu_comx,
      vars->gpu_comy,
      vars->gpu_comz,
      vars->gpu_r_k_x,
      vars->gpu_r_k_y,
      vars->gpu_r_k_z,
      gpu_moleculeInvolved,
      vars->gpu_cell_x[box],
      vars->gpu_cell_y[box],
      vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box],
      axis,
      halfAx,
      atomCount,
      r_max,
      step,
      seed,
      BETA,
      kill);
  } else {
    BrownianMotionRotateKernel<false><<< blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_startAtomIdx,
      vars->gpu_x,
      vars->gpu_y,
      vars->gpu_z,
      vars->gpu_mTorquex,
      vars->gpu_mTorquey,
      vars->gpu_mTorquez,
      vars->gpu_comx,
      vars->gpu_comy,
      vars->gpu_comz,
      vars->gpu_r_k_x,
      vars->gpu_r_k_y,
      vars->gpu_r_k_z,
      gpu_moleculeInvolved,
      vars->gpu_cell_x[box],
      vars->gpu_cell_y[box],
      vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box],
      axis,
      halfAx,
      atomCount,
      r_max,
      step,
      seed,
      BETA,
      kill);
  }

  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.x, vars->gpu_r_k_x, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.y, vars->gpu_r_k_y, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.z, vars->gpu_r_k_z, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  CUFREE(gpu_moleculeInvolved);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

template<bool isOrthogonal>
__global__ void BrownianMotionRotateKernel(
  int *startAtomIdx,
  double *gpu_x,
  double *gpu_y,
  double *gpu_z,
  double *molTorquex,
  double *molTorquey,
  double *molTorquez,
  double *gpu_comx,
  double *gpu_comy,
  double *gpu_comz,
  double *gpu_r_k_x,
  double *gpu_r_k_y,
  double *gpu_r_k_z,
  int *moleculeInvolved,
  double *gpu_cell_x,
  double *gpu_cell_y,
  double *gpu_cell_z,
  double *gpu_Invcell_x,
  double *gpu_Invcell_y,
  double *gpu_Invcell_z,
  double3 axis,
  double3 halfAx,
  int atomCount,
  double r_max,
  ulong step,
  ulong seed,
  double BETA,
  int *kill)
{
  //Each grid take cares of one molecule
  int molIndex = moleculeInvolved[blockIdx.x];
  int startIdx = startAtomIdx[molIndex];
  int endIdx = startAtomIdx[molIndex + 1];
  int atomIdx;

  __shared__ double matrix[3][3];
  __shared__ double3 com;

  // thread 0 will setup the matrix and update the gpu_r_k
  if(threadIdx.x == 0) {
    com = make_double3(gpu_comx[molIndex], gpu_comy[molIndex], gpu_comz[molIndex]);
    // This section calculates the amount of rotation
    double stdDev = sqrt(2.0 * r_max);
    double btm_x = molTorquex[molIndex] * BETA * r_max;
    double btm_y = molTorquey[molIndex] * BETA * r_max;
    double btm_z = molTorquez[molIndex] * BETA * r_max;

    double rot_x = btm_x + randomGaussianGPU(molIndex * 3, step, seed, 0.0, stdDev);
    double rot_y = btm_y + randomGaussianGPU(molIndex * 3 + 1, step, seed, 0.0, stdDev);
    double rot_z = btm_z + randomGaussianGPU(molIndex * 3 + 2, step, seed, 0.0, stdDev);
    // update the trial torque
    gpu_r_k_x[molIndex] = rot_x;
    gpu_r_k_y[molIndex] = rot_y;
    gpu_r_k_z[molIndex] = rot_z;
    //check for bad configuration
    if(!isfinite(rot_x + rot_y + rot_z)) {
      atomicAdd(kill, 1);
    }
    // build rotation matrix
    double cross[3][3], tensor[3][3];
    double rotLen = sqrt(rot_x * rot_x + rot_y * rot_y + rot_z * rot_z);
    double axisx = rot_x * (1.0 / rotLen);
    double axisy = rot_y * (1.0 / rotLen);
    double axisz = rot_z * (1.0 / rotLen);
    // build cross
    cross[0][0] = 0.0; cross[0][1] = -axisz; cross[0][2] = axisy;
    cross[1][0] = axisz; cross[1][1] = 0.0; cross[1][2] = -axisx;
    cross[2][0] = -axisy; cross[2][1] = axisx; cross[2][2] = 0.0;
    // build tensor
    int i, j;
    for(i = 0; i < 3; ++i) {
      tensor[0][i] = axisx;
      tensor[1][i] = axisy;
      tensor[2][i] = axisz;
    }
    for(i = 0; i < 3; ++i) {
      tensor[i][0] *= axisx;
      tensor[i][1] *= axisy;
      tensor[i][2] *= axisz;
    }
    // build matrix
    double s, c;
    sincos(rotLen, &s, &c);
    for(i = 0; i < 3; ++i) {
      for(j = 0; j < 3; ++j) {
        matrix[i][j] = 0.0;
      }
      matrix[i][i] = c;
    }
    for(i = 0; i < 3; ++i) {
      for(j = 0; j < 3; ++j) {
        matrix[i][j] += s * cross[i][j] + (1 - c) * tensor[i][j];
      }
    }
  }

  __syncthreads();
  // use strid of blockDim.x, which is 32
  // each thread handles one atom rotation
  for(atomIdx = startIdx + threadIdx.x; atomIdx < endIdx; atomIdx += blockDim.x) {
    double3 coor = make_double3(gpu_x[atomIdx], gpu_y[atomIdx], gpu_z[atomIdx]);
    // unwrap molecule
    if(isOrthogonal) {
      UnwrapPBC_f3(coor, com, axis, halfAx);
    } else {
      double3 unSlant = make_double3(0.0, 0.0, 0.0);
      TransformUnSlantGPU(unSlant, coor, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
      UnwrapPBC_f3(unSlant, com, axis, halfAx);
      TransformSlantGPU(coor, unSlant, gpu_cell_x, gpu_cell_y, gpu_cell_z);
    }

    // move COM of molecule to zero
    coor.x -= com.x;
    coor.y -= com.y;
    coor.z -= com.z;
    // rotate
    double newx = matrix[0][0] * coor.x + matrix[0][1] * coor.y + matrix[0][2] * coor.z;
    double newy = matrix[1][0] * coor.x + matrix[1][1] * coor.y + matrix[1][2] * coor.z;
    double newz = matrix[2][0] * coor.x + matrix[2][1] * coor.y + matrix[2][2] * coor.z;

    // move back to com
    coor.x = newx + com.x;
    coor.y = newy + com.y;
    coor.z = newz + com.z;

    // wrap again
    if(isOrthogonal) {
      WrapPBC_f3(coor, axis);
    } else {
      double3 unSlant = make_double3(0.0, 0.0, 0.0);
      TransformUnSlantGPU(unSlant, coor, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
      WrapPBC_f3(unSlant, axis);
      TransformSlantGPU(coor, unSlant, gpu_cell_x, gpu_cell_y, gpu_cell_z);
    }
    // update the new position
    gpu_x[atomIdx] = coor.x;
    gpu_y[atomIdx] = coor.y;
    gpu_z[atomIdx] = coor.z;
  }
}

void BrownianMotionTranslateParticlesGPU(
  VariablesCUDA *vars,
  const std::vector<unsigned int> &moleculeInvolved,
  XYZArray &mForce,
  XYZArray &mForceRec,
  XYZArray &newMolPos,
  XYZArray &newCOMs,
  XYZArray &t_k,
  const XYZ &boxAxes,
  const double BETA,
  const double t_max,
  ulong step,
  ulong seed,
  const int box,
  bool isOrthogonal,
  int *kill)
{
  int atomCount = newMolPos.Count();
  int molCount = newCOMs.Count();
  int molCountInBox = moleculeInvolved.size();
  int *gpu_moleculeInvolved;
  // Each block would handle one molecule
  int threadsPerBlock = 32;
  int blocksPerGrid = molCountInBox;

  CUMALLOC((void **) &gpu_moleculeInvolved, molCountInBox * sizeof(int));

  cudaMemcpy(vars->gpu_mForcex, mForce.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcey, mForce.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcez, mForce.z, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecx, mForceRec.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecy, mForceRec.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecz, mForceRec.z, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_moleculeInvolved, &moleculeInvolved[0], molCountInBox * sizeof(int), cudaMemcpyHostToDevice);

  double3 axis = make_double3(boxAxes.x, boxAxes.y, boxAxes.z);
  double3 halfAx = make_double3(boxAxes.x / 2.0, boxAxes.y / 2.0, boxAxes.z / 2.0);

  if(isOrthogonal) {
    BrownianMotionTranslateKernel<true><<< blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_startAtomIdx,
      vars->gpu_x,
      vars->gpu_y,
      vars->gpu_z,
      vars->gpu_mForcex,
      vars->gpu_mForcey,
      vars->gpu_mForcez,
      vars->gpu_mForceRecx,
      vars->gpu_mForceRecy,
      vars->gpu_mForceRecz,
      vars->gpu_comx,
      vars->gpu_comy,
      vars->gpu_comz,
      vars->gpu_t_k_x,
      vars->gpu_t_k_y,
      vars->gpu_t_k_z,
      gpu_moleculeInvolved,
      vars->gpu_cell_x[box],
      vars->gpu_cell_y[box],
      vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box],
      axis,
      halfAx,
      atomCount,
      t_max,
      step,
      seed,
      BETA,
      kill);
  } else {
    BrownianMotionTranslateKernel<false><<< blocksPerGrid, threadsPerBlock>>>(
      vars->gpu_startAtomIdx,
      vars->gpu_x,
      vars->gpu_y,
      vars->gpu_z,
      vars->gpu_mForcex,
      vars->gpu_mForcey,
      vars->gpu_mForcez,
      vars->gpu_mForceRecx,
      vars->gpu_mForceRecy,
      vars->gpu_mForceRecz,
      vars->gpu_comx,
      vars->gpu_comy,
      vars->gpu_comz,
      vars->gpu_t_k_x,
      vars->gpu_t_k_y,
      vars->gpu_t_k_z,
      gpu_moleculeInvolved,
      vars->gpu_cell_x[box],
      vars->gpu_cell_y[box],
      vars->gpu_cell_z[box],
      vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box],
      vars->gpu_Invcell_z[box],
      axis,
      halfAx,
      atomCount,
      t_max,
      step,
      seed,
      BETA,
      kill);
  }

  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.x, vars->gpu_comx, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.y, vars->gpu_comy, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.z, vars->gpu_comz, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.x, vars->gpu_t_k_x, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.y, vars->gpu_t_k_y, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.z, vars->gpu_t_k_z, molCount * sizeof(double), cudaMemcpyDeviceToHost);
  CUFREE(gpu_moleculeInvolved);
  checkLastErrorCUDA(__FILE__, __LINE__);
}


template<bool isOrthogonal>
__global__ void BrownianMotionTranslateKernel(
  int *startAtomIdx,
  double *gpu_x,
  double *gpu_y,
  double *gpu_z,
  double *molForcex,
  double *molForcey,
  double *molForcez,
  double *molForceRecx,
  double *molForceRecy,
  double *molForceRecz,
  double *gpu_comx,
  double *gpu_comy,
  double *gpu_comz,
  double *gpu_t_k_x,
  double *gpu_t_k_y,
  double *gpu_t_k_z,
  int *moleculeInvolved,
  double *gpu_cell_x,
  double *gpu_cell_y,
  double *gpu_cell_z,
  double *gpu_Invcell_x,
  double *gpu_Invcell_y,
  double *gpu_Invcell_z,
  double3 axis,
  double3 halfAx,
  int atomCount,
  double t_max,
  ulong step,
  ulong seed,
  double BETA,
  int *kill)
{
  //Each grid take cares of one molecule
  int molIndex = moleculeInvolved[blockIdx.x];
  int startIdx = startAtomIdx[molIndex];
  int endIdx = startAtomIdx[molIndex + 1];
  int atomIdx;

  __shared__ double3 shift;

  // thread 0 will calculate the shift vector and update COM and gpu_t_k
  if(threadIdx.x == 0) {
    double3 com = make_double3(gpu_comx[molIndex], gpu_comy[molIndex], gpu_comz[molIndex]);
    // This section calculates the amount of shift
    double stdDev = sqrt(2.0 * t_max);
    double bfm_x = (molForcex[molIndex] + molForceRecx[molIndex]) * BETA * t_max;
    double bfm_y = (molForcey[molIndex] + molForceRecy[molIndex]) * BETA * t_max;
    double bfm_z = (molForcez[molIndex] + molForceRecz[molIndex]) * BETA * t_max;

    shift.x = bfm_x + randomGaussianGPU(molIndex * 3, step, seed, 0.0, stdDev);
    shift.y = bfm_y + randomGaussianGPU(molIndex * 3 + 1, step, seed, 0.0, stdDev);
    shift.z = bfm_z + randomGaussianGPU(molIndex * 3 + 2, step, seed, 0.0, stdDev);
    // update the trial translate
    gpu_t_k_x[molIndex] = shift.x;
    gpu_t_k_y[molIndex] = shift.y;
    gpu_t_k_z[molIndex] = shift.z;
    // shift COM
    com.x += shift.x;
    com.y += shift.y;
    com.z += shift.z;
    // wrap COM
    if(isOrthogonal) {
      WrapPBC_f3(com, axis);
    } else {
      double3 unSlant = make_double3(0.0, 0.0, 0.0);
      TransformUnSlantGPU(unSlant, com, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
      WrapPBC_f3(unSlant, axis);
      TransformSlantGPU(com, unSlant, gpu_cell_x, gpu_cell_y, gpu_cell_z);
    }
    //update COM
    gpu_comx[molIndex] = com.x;
    gpu_comy[molIndex] = com.y;
    gpu_comz[molIndex] = com.z;
    //check for bad configuration
    if(!isfinite(shift.x + shift.y + shift.z)) {
      atomicAdd(kill, 1);
    } else if (shift.x > halfAx.x || shift.y > halfAx.y || shift.z > halfAx.z) {
      atomicAdd(kill, 1);
    }
  }

  __syncthreads();
  // use strid of blockDim.x, which is 32
  // each thread handles one atom translation
  for(atomIdx = startIdx + threadIdx.x; atomIdx < endIdx; atomIdx += blockDim.x) {
    double3 coor = make_double3(gpu_x[atomIdx], gpu_y[atomIdx], gpu_z[atomIdx]);

    // translate the atom
    coor.x += shift.x;
    coor.y += shift.y;
    coor.z += shift.z;
    // wrap coordinate
    if(isOrthogonal) {
      WrapPBC_f3(coor, axis);
    } else {
      double3 unSlant = make_double3(0.0, 0.0, 0.0);
      TransformUnSlantGPU(unSlant, coor, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
      WrapPBC_f3(unSlant, axis);
      TransformSlantGPU(coor, unSlant, gpu_cell_x, gpu_cell_y, gpu_cell_z);
    }
    // update the new position
    gpu_x[atomIdx] = coor.x;
    gpu_y[atomIdx] = coor.y;
    gpu_z[atomIdx] = coor.z;
  }
}

#endif
