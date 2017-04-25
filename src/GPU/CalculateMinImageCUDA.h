#pragma once

#ifdef GOMC_CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include "ConstantDefinitionsCUDA.h"

__device__ double MinImageSignedGpu(double raw,double ax, double halfAx) 
{
  if (raw > halfAx)
    raw -= ax;
  else if (raw < -halfAx)
    raw += ax;
  return raw;
}

__device__ bool InRcutGPU(double &distSq, double gpu_x1, double gpu_y1, 
			  double gpu_z1, double gpu_x2, double gpu_y2, 
			  double gpu_z2, double xAxes, double yAxes, 
			  double zAxes, double xHalfAxes, double yHalfAxes, 
			  double zHalfAxes)
{
  distSq = 0;
  double dx, dy, dz;
  dx = MinImageSigned(gpu_x1 - gpu_x2, xAxes, xHalfAxes);
  dy = MinImageSigned(gpu_y1 - gpu_y2, yAxes, yHalfAxes);
  dz = MinImageSigned(gpu_z1 - gpu_z2, zAxes, zHalfAxes);
  distSq = dx * dx + dy * dy + dz * dz;
  return ((gpu_rCut * gpu_rCut) > distSq);
}

__device__ bool InRcutGPU(double &distSq, double &virX, double &virY, 
			  double &virZ, double gpu_x1, double gpu_y1, 
			  double gpu_z1, double gpu_x2, double gpu_y2, 
			  double gpu_z2, double xAxes, double yAxes, 
			  double zAxes, double xHalfAxes, double yHalfAxes, 
			  double zHalfAxes)
{
  distSq = 0;
  virX = MinImageSigned(gpu_x1 - gpu_x2, xAxes, xHalfAxes);
  virY = MinImageSigned(gpu_y1 - gpu_y2, yAxes, yHalfAxes);
  virZ = MinImageSigned(gpu_z1 - gpu_z2, zAxes, zHalfAxes);
  distSq = virX * virX + virY * virY + virZ * virZ;
  return ((gpu_rCut * gpu_rCut) > distSq);
}

#endif
