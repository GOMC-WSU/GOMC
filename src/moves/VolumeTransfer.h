/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef VOLUMETRANSFER_H
#define VOLUMETRANSFER_H

#include "MoveBase.h" //For uint.

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == NPT

class VolumeTransfer : public MoveBase
{
public:
  VolumeTransfer(System &sys, StaticVals const& statV);

  virtual uint Prep(const real subDraw, const real movePerc);
  virtual void CalcEn();
  virtual uint Transform();
  real GetCoeff() const;
  virtual void Accept(const uint rejectState, const uint step);
  virtual void PrintAcceptKind();
private:
  //Note: This is only used for GEMC-NVT
  uint bPick[2];
  //Note: This is only used for GEMC-NPT and NPT
  uint box;
  SystemPotential sysPotNew;
  const Forcefield& forcefield;
  BoxDimensions newDim;
  BoxDimensionsNonOrth newDimNonOrth;
  Coordinates newMolsPos;
  COM newCOMs;
  MoleculeLookup & molLookRef;
  const uint GEMC_KIND;
  const real PRESSURE;
  bool regrewGrid, isOrth;
  bool fixBox0;
};

inline VolumeTransfer::VolumeTransfer(System &sys, StaticVals const& statV)  :
  MoveBase(sys, statV), molLookRef(sys.molLookupRef),
  newMolsPos(boxDimRef, newCOMs, sys.molLookupRef,
             sys.prng, statV.mol), forcefield(statV.forcefield),
  newDim(), newDimNonOrth(),
  newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef,
          statV.mol), GEMC_KIND(statV.kindOfGEMC),
  PRESSURE(statV.pressure), regrewGrid(false)
{
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(statV.mol.count);
  isOrth = statV.isOrthogonal;
  fixBox0 = statV.fixVolBox0;
}

void VolumeTransfer::PrintAcceptKind()
{
  printf("%-37s", "% Accepted Volume-Transfer ");
  for(uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::VOL_TRANSFER));
  }
  std::cout << std::endl;
}

inline uint VolumeTransfer::Prep(const real subDraw, const real movePerc)
{
  uint state = mv::fail_state::NO_FAIL;

  if (GEMC_KIND == mv::GEMC_NVT) {
    prng.PickBoxPair(bPick[0], bPick[1], subDraw, movePerc);
  } else {
    prng.PickBox(box, subDraw, movePerc);
    if(fixBox0) {
      //For NPT-GEMC and when box0 is fixed, we cannot pick box 0
      while(box == 0) {
        //To avoid infinit loop, we dont use sunDraw
        box = prng.randIntExc(BOX_TOTAL);
      }
    }
  }

  if(isOrth)
    newDim = boxDimRef;
  else
    newDimNonOrth = *((BoxDimensionsNonOrth*)(&boxDimRef));

  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
  return state;
}

inline uint VolumeTransfer::Transform()
{
  uint state = mv::fail_state::NO_FAIL;
  //Reinit, if necessary.
  if (GEMC_KIND == mv::GEMC_NVT) {
    real max = std::min(moveSetRef.Scale(bPick[0], mv::VOL_TRANSFER),
                          moveSetRef.Scale(bPick[1], mv::VOL_TRANSFER));
    if(isOrth) {
      coordCurrRef.VolumeTransferTranslate(state, newMolsPos, newCOMs, newDim,
                                           comCurrRef, max, bPick);
    } else {
      coordCurrRef.VolumeTransferTranslate(state, newMolsPos, newCOMs,
                                           newDimNonOrth, comCurrRef, max, bPick);
    }
  } else {
    //NPT ot GEMC-NPT we change volume of one box
    XYZ scale;
    real max = moveSetRef.Scale(box, mv::VOL_TRANSFER);
    real delta = prng.Sym(max);
    if(isOrth) {
      state =  boxDimRef.ShiftVolume(newDim, scale, box, delta);
    } else {
      state = boxDimRef.ShiftVolume(newDimNonOrth, scale, box, delta);
    }

    if (state == mv::fail_state::NO_FAIL) {
      if(isOrth) {
        coordCurrRef.TranslateOneBox(newMolsPos, newCOMs, comCurrRef,
                                     newDim, box, scale);
      } else {
        coordCurrRef.TranslateOneBox(newMolsPos, newCOMs, comCurrRef,
                                     newDimNonOrth, box, scale);
      }
    }
  }
  return state;
}

inline void VolumeTransfer::CalcEn()
{
  if (GEMC_KIND == mv::GEMC_NVT) {
    if(isOrth) {
      cellList.GridAll(newDim, newMolsPos, molLookRef);
    } else {
      cellList.GridAll(newDimNonOrth, newMolsPos, molLookRef);
    }
  } else {
    if(isOrth) {
      cellList.GridBox(newDim, newMolsPos, molLookRef, box);
    } else {
      cellList.GridBox(newDimNonOrth, newMolsPos, molLookRef, box);
    }
  }

  regrewGrid = true;
  //back up cached fourier term
  calcEwald->backupMolCache();
  sysPotNew = sysPotRef;

  if (GEMC_KIND == mv::GEMC_NVT) {
    for(uint b = 0; b < 2; b++) {
      //calculate new K vectors
      if(isOrth) {
        calcEwald->RecipInit(bPick[b], newDim);
        //setup reciprocate terms
        calcEwald->BoxReciprocalSetup(bPick[b], newMolsPos);
        sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, atomForceNew, molForceNew, newDim, bPick[b]);
      } else {
        calcEwald->RecipInit(bPick[b], newDimNonOrth);
        //setup reciprocate terms
        calcEwald->BoxReciprocalSetup(bPick[b], newMolsPos);
        sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos,
                                       atomForceNew, molForceNew, newDimNonOrth,
                                        bPick[b]);
      }
      //calculate reciprocate term of electrostatic interaction
      sysPotNew.boxEnergy[bPick[b]].recip = calcEwald->BoxReciprocal(bPick[b]);
    }
  } else {
    //calculate new K vectors
    if(isOrth) {
      calcEwald->RecipInit(box, newDim);
      //setup reciprocate terms
      calcEwald->BoxReciprocalSetup(box, newMolsPos);
      sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, atomForceNew, molForceNew, newDim, box);
    } else {
      calcEwald->RecipInit(box, newDimNonOrth);
      //setup reciprocate terms
      calcEwald->BoxReciprocalSetup(box, newMolsPos);
      sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, atomForceNew, molForceNew, newDimNonOrth, box);
    }
    //calculate reciprocate term of electrostatic interaction
    sysPotNew.boxEnergy[box].recip = calcEwald->BoxReciprocal(box);
  }

  sysPotNew.Total();
}



inline real VolumeTransfer::GetCoeff() const
{
  ////Log-volume style shift -- is turned off, at present.
  //
  //return pow(newDim.volume[b_i]/boxDimRef.volume[b_i],
  //	      (real)molLookRef.NumInBox(b_i)+1) *
  //  pow(newDim.volume[b_ii]/boxDimRef.volume[b_ii],
  //	 (real)molLookRef.NumInBox(b_ii)+1);
  real coeff = 1.0;
  if (GEMC_KIND == mv::GEMC_NVT) {
    for (uint b = 0; b < 2; ++b) {
      if(isOrth)
        coeff *= pow(newDim.volume[bPick[b]] / boxDimRef.volume[bPick[b]],
                     (real)molLookRef.NumInBox(bPick[b]));
      else
        coeff *= pow(newDimNonOrth.volume[bPick[b]] / boxDimRef.volume[bPick[b]],
                     (real)molLookRef.NumInBox(bPick[b]));
    }
  } else {
    if(isOrth) {
      coeff = pow(newDim.volume[box] / boxDimRef.volume[box],
                  (real)molLookRef.NumInBox(box)) *
              exp(-BETA * PRESSURE * (newDim.volume[box] - boxDimRef.volume[box]));
    } else {
      coeff = pow(newDimNonOrth.volume[box] / boxDimRef.volume[box],
                  (real)molLookRef.NumInBox(box)) *
              exp(-BETA * PRESSURE * (newDimNonOrth.volume[box] - boxDimRef.volume[box]));
    }
  }
  return coeff;
}

inline void VolumeTransfer::Accept(const uint rejectState, const uint step)
{
  real volTransCoeff = GetCoeff();
  real uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  real accept = volTransCoeff * uBoltz;
  bool result = (rejectState == mv::fail_state::NO_FAIL) && prng() < accept;
  if (result) {
    //Set new energy.
    sysPotRef = sysPotNew;
    //Swap... next time we'll use the current members.
    //NOTE:
    //This will be less efficient for NPT, but necessary evil.
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    swap(atomForceRef, atomForceNew);
    swap(molForceRef, molForceNew);
    if(isOrth)
      boxDimRef = newDim;
    else
      *((BoxDimensionsNonOrth*)(&boxDimRef)) = newDimNonOrth;

    if (GEMC_KIND == mv::GEMC_NVT) {
      for (uint b = 0; b < 2; b++) {
        calcEwald->UpdateRecip(bPick[b]);
        calcEwald->UpdateRecipVec(bPick[b]);
      }
    } else {
      calcEwald->UpdateRecip(box);
      calcEwald->UpdateRecipVec(box);
    }

  } else if (rejectState == mv::fail_state::NO_FAIL && regrewGrid) {
    if (GEMC_KIND == mv::GEMC_NVT) {
      cellList.GridAll(boxDimRef, coordCurrRef, molLookRef);
    } else {
      cellList.GridBox(boxDimRef, coordCurrRef, molLookRef, box);
    }

    regrewGrid = false;
    calcEwald->exgMolCache();
#ifdef GOMC_CUDA
    if (GEMC_KIND == mv::GEMC_NVT) {
      //update unitcell to the original in GPU
      for (uint b = 0; b < 2; b++) {
        UpdateCellBasisCUDA(forcefield.particles->getCUDAVars(), bPick[b],
                            boxDimRef.cellBasis[bPick[b]].x,
                            boxDimRef.cellBasis[bPick[b]].y,
                            boxDimRef.cellBasis[bPick[b]].z);
        if(!isOrth) {
          BoxDimensionsNonOrth newAxes = *((BoxDimensionsNonOrth*)(&boxDimRef));
          UpdateInvCellBasisCUDA(forcefield.particles->getCUDAVars(), bPick[b],
                                 newAxes.cellBasis_Inv[bPick[b]].x,
                                 newAxes.cellBasis_Inv[bPick[b]].y,
                                 newAxes.cellBasis_Inv[bPick[b]].z);
        }
      }
    } else {
      UpdateCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                          boxDimRef.cellBasis[box].x,
                          boxDimRef.cellBasis[box].y,
                          boxDimRef.cellBasis[box].z);
      if(!isOrth) {
        BoxDimensionsNonOrth newAxes = *((BoxDimensionsNonOrth*)(&boxDimRef));
        UpdateInvCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                               newAxes.cellBasis_Inv[box].x,
                               newAxes.cellBasis_Inv[box].y,
                               newAxes.cellBasis_Inv[box].z);
      }
    }
#endif
  }

  if (GEMC_KIND == mv::GEMC_NVT) {
    moveSetRef.Update(mv::VOL_TRANSFER, result, step, bPick[0]);
    moveSetRef.Update(mv::VOL_TRANSFER, result, step, bPick[1]);
  } else {
    moveSetRef.Update(mv::VOL_TRANSFER, result, step, box);
  }

}

#endif

#endif /*VOLUMETRANSFER_H*/
