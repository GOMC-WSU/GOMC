/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifdef GOMC_CUDA
#include "CUDAMemoryManager.cuh"
#include "CalculateMinImageCUDAKernel.cuh"
#include "Random123/boxmuller.hpp"
#include "TransformParticlesCUDAKernel.cuh"

#define MIN_FORCE 1E-12
#define MAX_FORCE 30

__device__ inline double randomGPU(unsigned int counter, ulong step,
                                   ulong seed) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  RNG::ctr_type r = philox4x64(c, k);
  return r123::u01<double>(r[0]);
}

__device__ inline double3 randomCoordsGPU(unsigned int counter,
                                          unsigned int key, ulong step,
                                          ulong seed) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  c[1] = key;
  RNG::ctr_type r = philox4x64(c, k);
  double3 r01;
  r01.x = r123::u01<double>(r[0]);
  r01.y = r123::u01<double>(r[1]);
  r01.z = r123::u01<double>(r[2]);
  return r01;
}

__device__ inline double randomGaussianGPU(unsigned int counter, ulong step,
                                           ulong seed, double mean,
                                           double stdDev) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  RNG::ctr_type r = philox4x64(c, k);
  double2 normal2 = r123::boxmuller(r[0], r[1]);
  double shiftedVal = mean + normal2.x * stdDev;
  return shiftedVal;
}

__device__ inline double3 randomGaussianCoordsGPU(unsigned int counter,
                                                  unsigned int key, ulong step,
                                                  ulong seed, double mean,
                                                  double stdDev) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  c[1] = key;
  RNG::ctr_type r = philox4x64(c, k);
  double2 normal1 = r123::boxmuller(r[0], r[1]);
  double2 normal2 = r123::boxmuller(r[2], r[3]);

  double3 normals =
      make_double3(mean + normal1.x * stdDev, mean + normal1.y * stdDev,
                   mean + normal2.x * stdDev);
  return normals;
}

__device__ inline double SymRandomGPU(unsigned int counter, unsigned int key,
                                      ulong step, ulong seed) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  c[1] = key;
  RNG::ctr_type r = philox4x64(c, k);
  double r01 = r123::uneg11<double>(r[0]);
  return r01;
}

__device__ inline double3 SymRandomCoordsGPU(unsigned int counter,
                                             unsigned int key, ulong step,
                                             ulong seed) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  c[1] = key;
  RNG::ctr_type r = philox4x64(c, k);
  double3 r01;
  r01.x = r123::uneg11<double>(r[0]);
  r01.y = r123::uneg11<double>(r[1]);
  r01.z = r123::uneg11<double>(r[2]);
  return r01;
}

// Returns a uniformly random point on the unit sphere
__device__ inline double3 RandomCoordsOnSphereGPU(unsigned int counter,
                                                  unsigned int key, ulong step,
                                                  ulong seed) {
  RNG::ctr_type c = {{}};
  RNG::ukey_type uk = {{}};
  uk[0] = step;
  uk[1] = seed;
  RNG::key_type k = uk;
  c[0] = counter;
  c[1] = key;
  RNG::ctr_type r = philox4x64(c, k);
  // picking phi uniformly will cluster points at poles
  // pick u = cos(phi) uniformly instead
  // start from r[1] because I used r[0] in GetSymRandom when called in
  // multiparticle
  double u = r123::uneg11<double>(r[1]);
  // theta must be [0, 2pi) !
  double theta = 2.0 * M_PI * r123::u01<double>(r[2]);
  double sintheta, costheta;
  sincos(theta, &sintheta, &costheta);
  double rootTerm = sqrt(1.0 - u * u);
  return make_double3(rootTerm * costheta, rootTerm * sintheta, u);
}

__device__ inline void
ApplyRotation(double &x, double &y, double &z, double comx, double comy,
              double comz, double theta, double3 rotvec, double axx, double axy,
              double axz, int gpu_nonOrth, double *gpu_cell_x,
              double *gpu_cell_y, double *gpu_cell_z, double *gpu_Invcell_x,
              double *gpu_Invcell_y, double *gpu_Invcell_z) {
  double matrix[3][3], cross[3][3], tensor[3][3];

  // build cross
  cross[0][0] = 0.0;
  cross[0][1] = -rotvec.z;
  cross[0][2] = rotvec.y;

  cross[1][0] = rotvec.z;
  cross[1][1] = 0.0;
  cross[1][2] = -rotvec.x;

  cross[2][0] = -rotvec.y;
  cross[2][1] = rotvec.x;
  cross[2][2] = 0.0;

  // build tensor
  for (int i = 0; i < 3; i++) {
    tensor[0][i] = rotvec.x;
    tensor[1][i] = rotvec.y;
    tensor[2][i] = rotvec.z;
  }
  for (int i = 0; i < 3; i++) {
    tensor[i][0] *= rotvec.x;
    tensor[i][1] *= rotvec.y;
    tensor[i][2] *= rotvec.z;
  }

  // build matrix
  double s, c;
  sincos(theta, &s, &c);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix[i][j] = 0.0;
    }
    matrix[i][i] = c;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix[i][j] += s * cross[i][j] + (1 - c) * tensor[i][j];
    }
  }

  // unwrap molecule
  double3 coor = make_double3(x, y, z);
  double3 com = make_double3(comx, comy, comz);
  double3 axes = make_double3(axx, axy, axz);
  double3 halfAx = make_double3(axx * 0.5, axy * 0.5, axz * 0.5);

  if (gpu_nonOrth)
    UnwrapPBCNonOrth3(coor, com, axes, halfAx, gpu_cell_x, gpu_cell_y,
                      gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
  else
    UnwrapPBC3(coor, com, axes, halfAx);

  // move particle to zero
  coor.x -= com.x;
  coor.y -= com.y;
  coor.z -= com.z;

  // rotate
  double3 newcoor;
  newcoor.x =
      matrix[0][0] * coor.x + matrix[0][1] * coor.y + matrix[0][2] * coor.z;
  newcoor.y =
      matrix[1][0] * coor.x + matrix[1][1] * coor.y + matrix[1][2] * coor.z;
  newcoor.z =
      matrix[2][0] * coor.x + matrix[2][1] * coor.y + matrix[2][2] * coor.z;

  coor.x = newcoor.x + com.x;
  coor.y = newcoor.y + com.y;
  coor.z = newcoor.z + com.z;

  // wrap again
  if (gpu_nonOrth)
    WrapPBCNonOrth3(coor, axes, gpu_cell_x, gpu_cell_y, gpu_cell_z,
                    gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
  else
    WrapPBC3(coor, axes);

  x = coor.x;
  y = coor.y;
  z = coor.z;
}

void CallTranslateParticlesGPU(
    VariablesCUDA *vars, const std::vector<int8_t> &isMoleculeInvolved, int box,
    double t_max, double *mForcex, double *mForcey, double *mForcez,
    std::vector<int> &inForceRange, ulong step, unsigned int key, ulong seed,
    const std::vector<int> &particleMol, int atomCount, int molCount,
    double xAxes, double yAxes, double zAxes, XYZArray &newMolPos,
    XYZArray &newCOMs, double lambdaBETA, XYZArray &t_k,
    XYZArray &molForceRecRef) {
  int8_t *gpu_isMoleculeInvolved;
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)(atomCount / threadsPerBlock) + 1;
  int *gpu_particleMol;

  CUMALLOC((void **)&gpu_isMoleculeInvolved,
           isMoleculeInvolved.size() * sizeof(int8_t));
  CUMALLOC((void **)&gpu_particleMol, particleMol.size() * sizeof(int));

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
             isMoleculeInvolved.size() * sizeof(int8_t),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);

  checkLastErrorCUDA(__FILE__, __LINE__);
  TranslateParticlesKernel<<<blocksPerGrid, threadsPerBlock>>>(
      molCount, t_max, vars->gpu_mForcex, vars->gpu_mForcey, vars->gpu_mForcez,
      vars->gpu_inForceRange, step, key, seed, vars->gpu_x, vars->gpu_y,
      vars->gpu_z, gpu_particleMol, atomCount, xAxes, yAxes, zAxes,
      vars->gpu_comx, vars->gpu_comy, vars->gpu_comz, vars->gpu_cell_x[box],
      vars->gpu_cell_y[box], vars->gpu_cell_z[box], vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box], vars->gpu_Invcell_z[box], vars->gpu_nonOrth,
      lambdaBETA, vars->gpu_t_k_x, vars->gpu_t_k_y, vars->gpu_t_k_z,
      gpu_isMoleculeInvolved, vars->gpu_mForceRecx, vars->gpu_mForceRecy,
      vars->gpu_mForceRecz);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.x, vars->gpu_comx, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.y, vars->gpu_comy, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.z, vars->gpu_comz, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.x, vars->gpu_t_k_x, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.y, vars->gpu_t_k_y, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.z, vars->gpu_t_k_z, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(&inForceRange[0], vars->gpu_inForceRange, molCount * sizeof(int),
             cudaMemcpyDeviceToHost);
  CUFREE(gpu_isMoleculeInvolved);
  CUFREE(gpu_particleMol);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

void CallRotateParticlesGPU(
    VariablesCUDA *vars, const std::vector<int8_t> &isMoleculeInvolved, int box,
    double r_max, double *mTorquex, double *mTorquey, double *mTorquez,
    std::vector<int> &inForceRange, ulong step, unsigned int key, ulong seed,
    const std::vector<int> &particleMol, int atomCount, int molCount,
    double xAxes, double yAxes, double zAxes, XYZArray &newMolPos,
    XYZArray &newCOMs, double lambdaBETA, XYZArray &r_k) {
  int8_t *gpu_isMoleculeInvolved;
  int threadsPerBlock = 256;
  int blocksPerGrid = (int)(atomCount / threadsPerBlock) + 1;
  int *gpu_particleMol;

  CUMALLOC((void **)&gpu_isMoleculeInvolved,
           isMoleculeInvolved.size() * sizeof(int8_t));
  CUMALLOC((void **)&gpu_particleMol, particleMol.size() * sizeof(int));

  cudaMemcpy(vars->gpu_mTorquex, mTorquex, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquey, mTorquey, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquez, mTorquez, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_particleMol, &particleMol[0], particleMol.size() * sizeof(int),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_isMoleculeInvolved, &isMoleculeInvolved[0],
             isMoleculeInvolved.size() * sizeof(int8_t),
             cudaMemcpyHostToDevice);

  RotateParticlesKernel<<<blocksPerGrid, threadsPerBlock>>>(
      molCount, r_max, vars->gpu_mTorquex, vars->gpu_mTorquey,
      vars->gpu_mTorquez, vars->gpu_inForceRange, step, key, seed, vars->gpu_x,
      vars->gpu_y, vars->gpu_z, gpu_particleMol, atomCount, xAxes, yAxes, zAxes,
      vars->gpu_comx, vars->gpu_comy, vars->gpu_comz, vars->gpu_cell_x[box],
      vars->gpu_cell_y[box], vars->gpu_cell_z[box], vars->gpu_Invcell_x[box],
      vars->gpu_Invcell_y[box], vars->gpu_Invcell_z[box], vars->gpu_nonOrth,
      lambdaBETA, vars->gpu_r_k_x, vars->gpu_r_k_y, vars->gpu_r_k_z,
      gpu_isMoleculeInvolved);
  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.x, vars->gpu_r_k_x, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.y, vars->gpu_r_k_y, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.z, vars->gpu_r_k_z, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(&inForceRange[0], vars->gpu_inForceRange, molCount * sizeof(int),
             cudaMemcpyDeviceToHost);
  CUFREE(gpu_isMoleculeInvolved);
  CUFREE(gpu_particleMol);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

__global__ void TranslateParticlesKernel(
    unsigned int numberOfMolecules, double t_max, double *molForcex,
    double *molForcey, double *molForcez, int *gpu_inForceRange, ulong step,
    unsigned int key, ulong seed, double *gpu_x, double *gpu_y, double *gpu_z,
    int *gpu_particleMol, int atomCount, double xAxes, double yAxes,
    double zAxes, double *gpu_comx, double *gpu_comy, double *gpu_comz,
    double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
    double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
    int *gpu_nonOrth, double lambdaBETA, double *gpu_t_k_x, double *gpu_t_k_y,
    double *gpu_t_k_z, int8_t *gpu_isMoleculeInvolved, double *gpu_mForceRecx,
    double *gpu_mForceRecy, double *gpu_mForceRecz) {
  int atomNumber = blockIdx.x * blockDim.x + threadIdx.x;
  if (atomNumber >= atomCount)
    return;

  int molIndex = gpu_particleMol[atomNumber];
  if (!gpu_isMoleculeInvolved[molIndex])
    return;
  bool updateMol = atomNumber == 0 || (gpu_particleMol[atomNumber] !=
                                       gpu_particleMol[atomNumber - 1]);

  // This section calculates the amount of shift
  double lbfx = (molForcex[molIndex] + gpu_mForceRecx[molIndex]) * lambdaBETA;
  double lbfy = (molForcey[molIndex] + gpu_mForceRecy[molIndex]) * lambdaBETA;
  double lbfz = (molForcez[molIndex] + gpu_mForceRecz[molIndex]) * lambdaBETA;
  double lbmaxx = lbfx * t_max;
  double lbmaxy = lbfy * t_max;
  double lbmaxz = lbfz * t_max;

  double shiftx, shifty, shiftz;
  bool forceInRange;

  forceInRange =
      (std::abs(lbmaxx) > MIN_FORCE && std::abs(lbmaxx) < MAX_FORCE &&
       std::abs(lbmaxy) > MIN_FORCE && std::abs(lbmaxy) < MAX_FORCE &&
       std::abs(lbmaxz) > MIN_FORCE && std::abs(lbmaxz) < MAX_FORCE);

  if (forceInRange) {
    double3 randnums = randomCoordsGPU(molIndex, key, step, seed);
    shiftx = log(exp(-1.0 * lbmaxx) + 2 * randnums.x * sinh(lbmaxx)) / lbfx;
    shifty = log(exp(-1.0 * lbmaxy) + 2 * randnums.y * sinh(lbmaxy)) / lbfy;
    shiftz = log(exp(-1.0 * lbmaxz) + 2 * randnums.z * sinh(lbmaxz)) / lbfz;
  } else {
    double3 randnums = SymRandomCoordsGPU(molIndex, key, step, seed);
    shiftx = t_max * randnums.x;
    shifty = t_max * randnums.y;
    shiftz = t_max * randnums.z;
  }

  // perform the shift on the coordinates
  double3 coor =
      make_double3(gpu_x[atomNumber] + shiftx, gpu_y[atomNumber] + shifty,
                   gpu_z[atomNumber] + shiftz);

  // rewrapping
  double3 axes = make_double3(xAxes, yAxes, zAxes);
  if (gpu_nonOrth[0])
    WrapPBCNonOrth3(coor, axes, gpu_cell_x, gpu_cell_y, gpu_cell_z,
                    gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
  else
    WrapPBC3(coor, axes);

  gpu_x[atomNumber] = coor.x;
  gpu_y[atomNumber] = coor.y;
  gpu_z[atomNumber] = coor.z;

  // update the CoM just once per molecule
  if (updateMol) {
    double3 com =
        make_double3(gpu_comx[molIndex] + shiftx, gpu_comy[molIndex] + shifty,
                     gpu_comz[molIndex] + shiftz);

    if (gpu_nonOrth[0])
      WrapPBCNonOrth3(com, axes, gpu_cell_x, gpu_cell_y, gpu_cell_z,
                      gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);
    else
      WrapPBC3(com, axes);

    gpu_comx[molIndex] = com.x;
    gpu_comy[molIndex] = com.y;
    gpu_comz[molIndex] = com.z;
    gpu_t_k_x[molIndex] = shiftx;
    gpu_t_k_y[molIndex] = shifty;
    gpu_t_k_z[molIndex] = shiftz;
    gpu_inForceRange[molIndex] = forceInRange;
  }
}

__global__ void RotateParticlesKernel(
    unsigned int numberOfMolecules, double r_max, double *molTorquex,
    double *molTorquey, double *molTorquez, int *gpu_inForceRange, ulong step,
    unsigned int key, ulong seed, double *gpu_x, double *gpu_y, double *gpu_z,
    int *gpu_particleMol, int atomCount, double xAxes, double yAxes,
    double zAxes, double *gpu_comx, double *gpu_comy, double *gpu_comz,
    double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
    double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
    int *gpu_nonOrth, double lambdaBETA, double *gpu_r_k_x, double *gpu_r_k_y,
    double *gpu_r_k_z, int8_t *gpu_isMoleculeInvolved) {
  int atomNumber = blockIdx.x * blockDim.x + threadIdx.x;
  if (atomNumber >= atomCount)
    return;
  int molIndex = gpu_particleMol[atomNumber];
  if (!gpu_isMoleculeInvolved[molIndex])
    return;
  bool updateMol = atomNumber == 0 || (gpu_particleMol[atomNumber] !=
                                       gpu_particleMol[atomNumber - 1]);

  // This section calculates the amount of rotation
  double lbtx = molTorquex[molIndex] * lambdaBETA;
  double lbty = molTorquey[molIndex] * lambdaBETA;
  double lbtz = molTorquez[molIndex] * lambdaBETA;
  double lbmaxx = lbtx * r_max;
  double lbmaxy = lbty * r_max;
  double lbmaxz = lbtz * r_max;

  double rotx, roty, rotz, theta;
  double3 rotvec;
  bool forceInRange;

  forceInRange =
      (std::abs(lbmaxx) > MIN_FORCE && std::abs(lbmaxx) < MAX_FORCE &&
       std::abs(lbmaxy) > MIN_FORCE && std::abs(lbmaxy) < MAX_FORCE &&
       std::abs(lbmaxz) > MIN_FORCE && std::abs(lbmaxz) < MAX_FORCE);

  if (forceInRange) {
    double3 randnums = randomCoordsGPU(molIndex, key, step, seed);
    rotx = log(exp(-1.0 * lbmaxx) + 2 * randnums.x * sinh(lbmaxx)) / lbtx;
    roty = log(exp(-1.0 * lbmaxy) + 2 * randnums.y * sinh(lbmaxy)) / lbty;
    rotz = log(exp(-1.0 * lbmaxz) + 2 * randnums.z * sinh(lbmaxz)) / lbtz;
    theta = sqrt(rotx * rotx + roty * roty + rotz * rotz);
    rotvec = make_double3(rotx * (1.0 / theta), roty * (1.0 / theta),
                          rotz * (1.0 / theta));
  } else {
    double3 randnums = RandomCoordsOnSphereGPU(molIndex, key, step, seed);
    // These values are ignored if !forceInRange so just initialize to zero.
    rotx = 0.0;
    roty = 0.0;
    rotz = 0.0;
    theta = r_max * SymRandomGPU(molIndex, key, step, seed);
    rotvec = randnums;
  }

  if (updateMol) {
    gpu_r_k_x[molIndex] = rotx;
    gpu_r_k_y[molIndex] = roty;
    gpu_r_k_z[molIndex] = rotz;
    gpu_inForceRange[molIndex] = forceInRange;
  }

  // perform the rotation on the coordinates
  ApplyRotation(gpu_x[atomNumber], gpu_y[atomNumber], gpu_z[atomNumber],
                gpu_comx[molIndex], gpu_comy[molIndex], gpu_comz[molIndex],
                theta, rotvec, xAxes, yAxes, zAxes, *gpu_nonOrth, gpu_cell_x,
                gpu_cell_y, gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                gpu_Invcell_z);
}

// CUDA implementation of MultiParticle Brownian transformation

void BrownianMotionRotateParticlesGPU(
    VariablesCUDA *vars, const std::vector<unsigned int> &moleculeInvolved,
    XYZArray &mTorque, XYZArray &newMolPos, XYZArray &newCOMs, XYZArray &r_k,
    const XYZ &boxAxes, const double BETA, const double r_max, ulong step,
    unsigned int key, ulong seed, const int box, const bool isOrthogonal,
    int *kill) {
  int atomCount = newMolPos.Count();
  int molCount = newCOMs.Count();
  int molCountInBox = moleculeInvolved.size();
  int *gpu_moleculeInvolved;
  // Each block would handle one molecule
  int threadsPerBlock = 32;
  int blocksPerGrid = molCountInBox;

  CUMALLOC((void **)&gpu_moleculeInvolved, molCountInBox * sizeof(int));

  cudaMemcpy(vars->gpu_mTorquex, mTorque.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquey, mTorque.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mTorquez, mTorque.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_moleculeInvolved, &moleculeInvolved[0],
             molCountInBox * sizeof(int), cudaMemcpyHostToDevice);

  double3 axis = make_double3(boxAxes.x, boxAxes.y, boxAxes.z);
  double3 halfAx =
      make_double3(boxAxes.x * 0.5, boxAxes.y * 0.5, boxAxes.z * 0.5);

  if (isOrthogonal)
    BrownianMotionRotateKernel<true><<<blocksPerGrid, threadsPerBlock>>>(
        vars->gpu_startAtomIdx, vars->gpu_x, vars->gpu_y, vars->gpu_z,
        vars->gpu_mTorquex, vars->gpu_mTorquey, vars->gpu_mTorquez,
        vars->gpu_comx, vars->gpu_comy, vars->gpu_comz, vars->gpu_r_k_x,
        vars->gpu_r_k_y, vars->gpu_r_k_z, gpu_moleculeInvolved,
        vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
        vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
        vars->gpu_Invcell_z[box], axis, halfAx, atomCount, r_max, step, key,
        seed, BETA, kill);
  else
    BrownianMotionRotateKernel<true><<<blocksPerGrid, threadsPerBlock>>>(
        vars->gpu_startAtomIdx, vars->gpu_x, vars->gpu_y, vars->gpu_z,
        vars->gpu_mTorquex, vars->gpu_mTorquey, vars->gpu_mTorquez,
        vars->gpu_comx, vars->gpu_comy, vars->gpu_comz, vars->gpu_r_k_x,
        vars->gpu_r_k_y, vars->gpu_r_k_z, gpu_moleculeInvolved,
        vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
        vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
        vars->gpu_Invcell_z[box], axis, halfAx, atomCount, r_max, step, key,
        seed, BETA, kill);

  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.x, vars->gpu_r_k_x, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.y, vars->gpu_r_k_y, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(r_k.z, vars->gpu_r_k_z, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  CUFREE(gpu_moleculeInvolved);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

template <const bool isOrthogonal>
__global__ void BrownianMotionRotateKernel(
    int *startAtomIdx, double *gpu_x, double *gpu_y, double *gpu_z,
    double *molTorquex, double *molTorquey, double *molTorquez,
    double *gpu_comx, double *gpu_comy, double *gpu_comz, double *gpu_r_k_x,
    double *gpu_r_k_y, double *gpu_r_k_z, int *moleculeInvolved,
    double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
    double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
    double3 axis, double3 halfAx, int atomCount, double r_max, ulong step,
    unsigned int key, ulong seed, double BETA, int *kill) {
  // Each block takes care of one molecule
  int molIndex = moleculeInvolved[blockIdx.x];
  int startIdx = startAtomIdx[molIndex];
  int endIdx = startAtomIdx[molIndex + 1];
  int atomIdx;

  __shared__ double matrix[3][3];
  __shared__ double3 com;

  // thread 0 will set up the matrix and update the gpu_r_k
  if (threadIdx.x == 0) {
    com = make_double3(gpu_comx[molIndex], gpu_comy[molIndex],
                       gpu_comz[molIndex]);
    // This section calculates the amount of rotation
    double stdDev = sqrt(2.0 * r_max);
    double btm_x = molTorquex[molIndex] * BETA * r_max;
    double btm_y = molTorquey[molIndex] * BETA * r_max;
    double btm_z = molTorquez[molIndex] * BETA * r_max;

    double3 randnums =
        randomGaussianCoordsGPU(molIndex, key, step, seed, 0.0, stdDev);
    double rot_x = btm_x + randnums.x;
    double rot_y = btm_y + randnums.y;
    double rot_z = btm_z + randnums.z;
    // update the trial torque
    gpu_r_k_x[molIndex] = rot_x;
    gpu_r_k_y[molIndex] = rot_y;
    gpu_r_k_z[molIndex] = rot_z;
    // check for bad configuration
    if (!isfinite(rot_x + rot_y + rot_z)) {
      atomicAdd(kill, 1);
    }
    // build rotation matrix
    double cross[3][3], tensor[3][3];
    double rotLen = sqrt(rot_x * rot_x + rot_y * rot_y + rot_z * rot_z);
    double axisx = rot_x * (1.0 / rotLen);
    double axisy = rot_y * (1.0 / rotLen);
    double axisz = rot_z * (1.0 / rotLen);
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
    int i, j;
    for (i = 0; i < 3; ++i) {
      tensor[0][i] = axisx;
      tensor[1][i] = axisy;
      tensor[2][i] = axisz;
    }
    for (i = 0; i < 3; ++i) {
      tensor[i][0] *= axisx;
      tensor[i][1] *= axisy;
      tensor[i][2] *= axisz;
    }
    // build matrix
    double s, c;
    sincos(rotLen, &s, &c);
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
        matrix[i][j] = 0.0;
      }
      matrix[i][i] = c;
    }
    for (i = 0; i < 3; ++i) {
      for (j = 0; j < 3; ++j) {
        matrix[i][j] += s * cross[i][j] + (1 - c) * tensor[i][j];
      }
    }
  }

  __syncthreads();
  // use stride of blockDim.x, which is 32
  // each thread handles one atom rotation
  for (atomIdx = startIdx + threadIdx.x; atomIdx < endIdx;
       atomIdx += blockDim.x) {
    double3 coor = make_double3(gpu_x[atomIdx], gpu_y[atomIdx], gpu_z[atomIdx]);
    // unwrap molecule
    if (isOrthogonal)
      UnwrapPBC3(coor, com, axis, halfAx);
    else
      UnwrapPBCNonOrth3(coor, com, axis, halfAx, gpu_cell_x, gpu_cell_y,
                        gpu_cell_z, gpu_Invcell_x, gpu_Invcell_y,
                        gpu_Invcell_z);

    // move COM of molecule to zero
    coor.x -= com.x;
    coor.y -= com.y;
    coor.z -= com.z;
    // rotate
    double newx =
        matrix[0][0] * coor.x + matrix[0][1] * coor.y + matrix[0][2] * coor.z;
    double newy =
        matrix[1][0] * coor.x + matrix[1][1] * coor.y + matrix[1][2] * coor.z;
    double newz =
        matrix[2][0] * coor.x + matrix[2][1] * coor.y + matrix[2][2] * coor.z;

    // move back to com
    coor.x = newx + com.x;
    coor.y = newy + com.y;
    coor.z = newz + com.z;

    // wrap again
    if (isOrthogonal)
      WrapPBC3(coor, axis);
    else
      WrapPBCNonOrth3(coor, axis, gpu_cell_x, gpu_cell_y, gpu_cell_z,
                      gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);

    // update the new position
    gpu_x[atomIdx] = coor.x;
    gpu_y[atomIdx] = coor.y;
    gpu_z[atomIdx] = coor.z;
  }
}

void BrownianMotionTranslateParticlesGPU(
    VariablesCUDA *vars, const std::vector<unsigned int> &moleculeInvolved,
    XYZArray &mForce, XYZArray &mForceRec, XYZArray &newMolPos,
    XYZArray &newCOMs, XYZArray &t_k, const XYZ &boxAxes, const double BETA,
    const double t_max, ulong step, unsigned int key, ulong seed, const int box,
    const bool isOrthogonal, int *kill) {
  int atomCount = newMolPos.Count();
  int molCount = newCOMs.Count();
  int molCountInBox = moleculeInvolved.size();
  int *gpu_moleculeInvolved;
  // Each block would handle one molecule
  int threadsPerBlock = 32;
  int blocksPerGrid = molCountInBox;

  CUMALLOC((void **)&gpu_moleculeInvolved, molCountInBox * sizeof(int));

  cudaMemcpy(vars->gpu_mForcex, mForce.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcey, mForce.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForcez, mForce.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecx, mForceRec.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecy, mForceRec.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_mForceRecz, mForceRec.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_x, newMolPos.x, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_y, newMolPos.y, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_z, newMolPos.z, atomCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comx, newCOMs.x, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comy, newCOMs.y, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(vars->gpu_comz, newCOMs.z, molCount * sizeof(double),
             cudaMemcpyHostToDevice);
  cudaMemcpy(gpu_moleculeInvolved, &moleculeInvolved[0],
             molCountInBox * sizeof(int), cudaMemcpyHostToDevice);

  double3 axis = make_double3(boxAxes.x, boxAxes.y, boxAxes.z);
  double3 halfAx =
      make_double3(boxAxes.x * 0.5, boxAxes.y * 0.5, boxAxes.z * 0.5);

  if (isOrthogonal)
    BrownianMotionTranslateKernel<true><<<blocksPerGrid, threadsPerBlock>>>(
        vars->gpu_startAtomIdx, vars->gpu_x, vars->gpu_y, vars->gpu_z,
        vars->gpu_mForcex, vars->gpu_mForcey, vars->gpu_mForcez,
        vars->gpu_mForceRecx, vars->gpu_mForceRecy, vars->gpu_mForceRecz,
        vars->gpu_comx, vars->gpu_comy, vars->gpu_comz, vars->gpu_t_k_x,
        vars->gpu_t_k_y, vars->gpu_t_k_z, gpu_moleculeInvolved,
        vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
        vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
        vars->gpu_Invcell_z[box], axis, halfAx, atomCount, t_max, step, key,
        seed, BETA, kill);
  else
    BrownianMotionTranslateKernel<false><<<blocksPerGrid, threadsPerBlock>>>(
        vars->gpu_startAtomIdx, vars->gpu_x, vars->gpu_y, vars->gpu_z,
        vars->gpu_mForcex, vars->gpu_mForcey, vars->gpu_mForcez,
        vars->gpu_mForceRecx, vars->gpu_mForceRecy, vars->gpu_mForceRecz,
        vars->gpu_comx, vars->gpu_comy, vars->gpu_comz, vars->gpu_t_k_x,
        vars->gpu_t_k_y, vars->gpu_t_k_z, gpu_moleculeInvolved,
        vars->gpu_cell_x[box], vars->gpu_cell_y[box], vars->gpu_cell_z[box],
        vars->gpu_Invcell_x[box], vars->gpu_Invcell_y[box],
        vars->gpu_Invcell_z[box], axis, halfAx, atomCount, t_max, step, key,
        seed, BETA, kill);

  cudaDeviceSynchronize();
  checkLastErrorCUDA(__FILE__, __LINE__);

  cudaMemcpy(newMolPos.x, vars->gpu_x, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.y, vars->gpu_y, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newMolPos.z, vars->gpu_z, atomCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.x, vars->gpu_comx, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.y, vars->gpu_comy, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(newCOMs.z, vars->gpu_comz, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.x, vars->gpu_t_k_x, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.y, vars->gpu_t_k_y, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  cudaMemcpy(t_k.z, vars->gpu_t_k_z, molCount * sizeof(double),
             cudaMemcpyDeviceToHost);
  CUFREE(gpu_moleculeInvolved);
  checkLastErrorCUDA(__FILE__, __LINE__);
}

template <const bool isOrthogonal>
__global__ void BrownianMotionTranslateKernel(
    int *startAtomIdx, double *gpu_x, double *gpu_y, double *gpu_z,
    double *molForcex, double *molForcey, double *molForcez,
    double *molForceRecx, double *molForceRecy, double *molForceRecz,
    double *gpu_comx, double *gpu_comy, double *gpu_comz, double *gpu_t_k_x,
    double *gpu_t_k_y, double *gpu_t_k_z, int *moleculeInvolved,
    double *gpu_cell_x, double *gpu_cell_y, double *gpu_cell_z,
    double *gpu_Invcell_x, double *gpu_Invcell_y, double *gpu_Invcell_z,
    double3 axis, double3 halfAx, int atomCount, double t_max, ulong step,
    unsigned int key, ulong seed, double BETA, int *kill) {
  // Each block takes care of one molecule
  int molIndex = moleculeInvolved[blockIdx.x];
  int startIdx = startAtomIdx[molIndex];
  int endIdx = startAtomIdx[molIndex + 1];
  int atomIdx;

  __shared__ double3 shift;

  // thread 0 will calculate the shift vector and update COM and gpu_t_k
  if (threadIdx.x == 0) {
    double3 com = make_double3(gpu_comx[molIndex], gpu_comy[molIndex],
                               gpu_comz[molIndex]);
    // This section calculates the amount of shift
    double stdDev = sqrt(2.0 * t_max);
    double bfm_x =
        (molForcex[molIndex] + molForceRecx[molIndex]) * BETA * t_max;
    double bfm_y =
        (molForcey[molIndex] + molForceRecy[molIndex]) * BETA * t_max;
    double bfm_z =
        (molForcez[molIndex] + molForceRecz[molIndex]) * BETA * t_max;

    double3 randnums =
        randomGaussianCoordsGPU(molIndex, key, step, seed, 0.0, stdDev);
    shift.x = bfm_x + randnums.x;
    shift.y = bfm_y + randnums.y;
    shift.z = bfm_z + randnums.z;
    // update the trial translate
    gpu_t_k_x[molIndex] = shift.x;
    gpu_t_k_y[molIndex] = shift.y;
    gpu_t_k_z[molIndex] = shift.z;
    // shift COM
    com.x += shift.x;
    com.y += shift.y;
    com.z += shift.z;
    // wrap COM
    if (isOrthogonal)
      WrapPBC3(com, axis);
    else
      WrapPBCNonOrth3(com, axis, gpu_cell_x, gpu_cell_y, gpu_cell_z,
                      gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);

    // update COM
    gpu_comx[molIndex] = com.x;
    gpu_comy[molIndex] = com.y;
    gpu_comz[molIndex] = com.z;
    // check for bad configuration
    if (!isfinite(shift.x + shift.y + shift.z)) {
      atomicAdd(kill, 1);
    } else if (shift.x > halfAx.x || shift.y > halfAx.y || shift.z > halfAx.z) {
      atomicAdd(kill, 1);
    }
  }

  __syncthreads();
  // use stride of blockDim.x, which is 32
  // each thread handles one atom translation
  for (atomIdx = startIdx + threadIdx.x; atomIdx < endIdx;
       atomIdx += blockDim.x) {
    double3 coor = make_double3(gpu_x[atomIdx], gpu_y[atomIdx], gpu_z[atomIdx]);

    // translate the atom
    coor.x += shift.x;
    coor.y += shift.y;
    coor.z += shift.z;
    // wrap coordinate
    if (isOrthogonal)
      WrapPBC3(coor, axis);
    else
      WrapPBCNonOrth3(coor, axis, gpu_cell_x, gpu_cell_y, gpu_cell_z,
                      gpu_Invcell_x, gpu_Invcell_y, gpu_Invcell_z);

    // update the new position
    gpu_x[atomIdx] = coor.x;
    gpu_y[atomIdx] = coor.y;
    gpu_z[atomIdx] = coor.z;
  }
}

#endif
