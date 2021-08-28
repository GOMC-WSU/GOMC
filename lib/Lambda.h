/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef LAMBDA_H
#define LAMBDA_H

#include "BasicTypes.h" //For ulong, uint

#ifdef GOMC_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include "VariablesCUDA.cuh"
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif
//Defining lambda class to handle fractional molecule
class Lambda
{
public:
  Lambda() {}

  void Init(
  #ifdef GOMC_CUDA 
    VariablesCUDA *refVarCUDA
  #endif
  ) {
    std::fill_n(isFraction, BOX_TOTAL, false);
    std::fill_n(lambdaVDW, BOX_TOTAL, 1.0);
    std::fill_n(lambdaCoulomb, BOX_TOTAL, 1.0);
    std::fill_n(kindIndex, BOX_TOTAL, 0);
    std::fill_n(molIndex, BOX_TOTAL, 0);
#ifdef GOMC_CUDA
    varCUDA = refVarCUDA;
    // Update Lambda on GPU
    UpdateGPULambda(varCUDA, molIndex, lambdaVDW,
                    lambdaCoulomb, isFraction);
#endif
  }

  void Set(const double vdw, const double coulomb, const uint mol,
           const uint kind, const uint box);

  void UnSet(const uint sourceBox, const uint destBox);

  double GetLambdaVDW(const uint mol, const uint box) const;

  double GetLambdaCoulomb(const uint mol, const uint box) const;

  bool KindIsFractional(const uint kind, const uint box) const;

  uint GetMolIndex(const uint box) const;


protected:
  int molIndex[BOX_TOTAL]; // global molecule index
  int kindIndex[BOX_TOTAL];
  double lambdaVDW[BOX_TOTAL], lambdaCoulomb[BOX_TOTAL];
  bool isFraction[BOX_TOTAL];
#ifdef GOMC_CUDA
  VariablesCUDA *varCUDA;
#endif
};

inline void Lambda::Set(const double vdw, const double coulomb, const uint mol,
                        const uint kind, const uint box)
{
  molIndex[box] = mol;
  kindIndex[box] = kind;
  lambdaVDW[box] = vdw;
  lambdaCoulomb[box] = coulomb;
  isFraction[box] = true;
#ifdef GOMC_CUDA
  // Update Lambda on GPU
  UpdateGPULambda(varCUDA, molIndex, lambdaVDW,
                  lambdaCoulomb, isFraction);
#endif
}

inline void Lambda::UnSet(const uint sourceBox, const uint destBox)
{
  isFraction[sourceBox] = isFraction[destBox] = false;
  lambdaVDW[sourceBox] = lambdaCoulomb[sourceBox] = 1.0;
  lambdaVDW[destBox] = lambdaCoulomb[destBox] = 0.0;
  molIndex[sourceBox] = molIndex[destBox] = 0;
  kindIndex[sourceBox] = kindIndex[destBox] = 0;
#ifdef GOMC_CUDA
  // Update Lambda on GPU
  UpdateGPULambda(varCUDA, molIndex, lambdaVDW,
                  lambdaCoulomb, isFraction);
#endif
}


inline double Lambda::GetLambdaVDW(const uint mol, const uint box) const
{
  double val = 1.0;
  if(isFraction[box]) {
    if(molIndex[box] == mol) {
      val = lambdaVDW[box];
    }
  }
  return val;
}

inline double Lambda::GetLambdaCoulomb(const uint mol, const uint box) const
{
  double val = 1.0;
  if(isFraction[box]) {
    if(molIndex[box] == mol) {
      val = lambdaCoulomb[box];
    }
  }
  return val;
}

inline bool Lambda::KindIsFractional(const uint kind, const uint box) const
{
  bool result = false;
  if(isFraction[box]) {
    if(kindIndex[box] == kind) {
      result = true;
    }
  }
  return result;
}

inline uint Lambda::GetMolIndex(const uint box) const
{
  if(isFraction[box]) {
    return molIndex[box];
  } else {
    std::cout << "Error: Lambda.h, calling GetMolIndex\n";
    exit(EXIT_FAILURE);
  }
}

#endif