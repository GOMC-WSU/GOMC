#pragma once
#ifdef GOMC_CUDA
#include <vector>
#include "Random123/philox.h"
typedef r123::Philox4x32 RNG;

using namespace std;

#include <cuda.h>
#include <cuda_runtime.h>
#include "VariablesCUDA.cuh"
#include "XYZArray.h"

void CallTranslateParticlesGPU(VariablesCUDA *vars,
                               vector<uint> &moleculeIndex,
                               uint moveType,
                               double t_max,
                               double *mForcex,
                               double *mForcey,
                               double *mForcez,
                               unsigned int step,
                               unsigned int seed,
                               vector<int> particleMol,
                               int atomCount,
                               int molCount,
                               double xAxes,
                               double yAxes,
                               double zAxes,
                               XYZArray &newMolPos,
                               XYZArray &newCOMs,
                               double lambdaBETA);

void CallRotateParticlesGPU(VariablesCUDA *vars,
                            vector<uint> &moleculeIndex,
                            uint moveType,
                            double r_max,
                            double *mTorquex,
                            double *mTorquey,
                            double *mTorquez,
                            unsigned int step,
                            unsigned int seed,
                            vector<int> particleMol,
                            int atomCount,
                            int molCount,
                            double xAxes,
                            double yAxes,
                            double zAxes,
                            XYZArray &newMolPos,
                            XYZArray &newCOMs);

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
                                         double lambdaBETA);

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
                                      double *gpu_comz);
#endif