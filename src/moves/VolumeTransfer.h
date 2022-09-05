/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef VOLUMETRANSFER_H
#define VOLUMETRANSFER_H

#include "MoveBase.h" //For uint.

#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif

#if ENSEMBLE == GEMC || ENSEMBLE == NPT

class VolumeTransfer : public MoveBase {
public:
  VolumeTransfer(System &sys, StaticVals const &statV);

  virtual uint Prep(const double subDraw, const double movePerc);
  // To relax the system in NE_MTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0,
                          const uint kidx = 0) {
    return mv::fail_state::NO_FAIL;
  }
  virtual void CalcEn();
  virtual uint Transform();
  double GetCoeff() const;
  virtual void Accept(const uint rejectState, const ulong step);
  virtual void PrintAcceptKind();

private:
  // Note: This is only used for GEMC-NVT
  uint bPick[2];
  // Note: This is only used for GEMC-NPT and NPT
  uint box;
  SystemPotential sysPotNew;
  MoleculeLookup &molLookRef;
  BoxDimensions newDim;
  BoxDimensionsNonOrth newDimNonOrth;
  Coordinates newMolsPos;
  const Forcefield &forcefield;
  COM newCOMs;
  const uint GEMC_KIND;
  const double PRESSURE;
  bool regrewGrid, isOrth;
  bool fixBox0;
};

inline VolumeTransfer::VolumeTransfer(System &sys, StaticVals const &statV)
    : MoveBase(sys, statV), molLookRef(sys.molLookupRef), newDim(),
      newDimNonOrth(),
      newMolsPos(boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
      forcefield(statV.forcefield),
      newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef, statV.mol),
      GEMC_KIND(statV.kindOfGEMC), PRESSURE(statV.pressure), regrewGrid(false) {
  newMolsPos.Init(sys.coordinates.Count());
  newCOMs.Init(statV.mol.count);
  isOrth = statV.isOrthogonal;
  fixBox0 = statV.fixVolBox0;
}

void VolumeTransfer::PrintAcceptKind() {
  printf("%-37s", "% Accepted Volume-Transfer ");
  for (uint b = 0; b < BOX_TOTAL; b++) {
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(b, mv::VOL_TRANSFER));
  }
  std::cout << std::endl;
}

inline uint VolumeTransfer::Prep(const double subDraw, const double movePerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_VOL_TRANSFER);
  uint state = mv::fail_state::NO_FAIL;

  if (GEMC_KIND == mv::GEMC_NVT) {
    prng.PickBoxPair(bPick[0], bPick[1], subDraw, movePerc);
  } else {
    prng.PickBox(box, subDraw, movePerc);
    if (fixBox0) {
      // For NPT-GEMC and when box0 is fixed, we cannot pick box 0
      while (box == 0) {
        // To avoid infinite loop, we don't use sunDraw
        box = prng.randIntExc(BOX_TOTAL);
      }
    }
  }

  if (isOrth)
    newDim = boxDimRef;
  else
    newDimNonOrth = *(static_cast<BoxDimensionsNonOrth *>(&boxDimRef));

  coordCurrRef.CopyRange(newMolsPos, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(newCOMs, 0, 0, comCurrRef.Count());
  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_VOL_TRANSFER);
  return state;
}

inline uint VolumeTransfer::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_VOL_TRANSFER);
  uint state = mv::fail_state::NO_FAIL;
  // Reinit, if necessary.
  if (GEMC_KIND == mv::GEMC_NVT) {
    double max = std::min(moveSetRef.Scale(bPick[0], mv::VOL_TRANSFER),
                          moveSetRef.Scale(bPick[1], mv::VOL_TRANSFER));
    if (isOrth) {
      coordCurrRef.VolumeTransferTranslate(state, newMolsPos, newCOMs, newDim,
                                           comCurrRef, max, bPick);
    } else {
      coordCurrRef.VolumeTransferTranslate(
          state, newMolsPos, newCOMs, newDimNonOrth, comCurrRef, max, bPick);
    }
  } else {
    // NPT or GEMC-NPT we change volume of one box
    XYZ scale;
    double max = moveSetRef.Scale(box, mv::VOL_TRANSFER);
    double delta = prng.Sym(max);
    if (isOrth) {
      state = boxDimRef.ShiftVolume(newDim, scale, box, delta);
    } else {
      state = boxDimRef.ShiftVolume(newDimNonOrth, scale, box, delta);
    }

    if (state == mv::fail_state::NO_FAIL) {
      if (isOrth) {
        coordCurrRef.TranslateOneBox(newMolsPos, newCOMs, comCurrRef, newDim,
                                     box, scale);
      } else {
        coordCurrRef.TranslateOneBox(newMolsPos, newCOMs, comCurrRef,
                                     newDimNonOrth, box, scale);
      }
    }
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_VOL_TRANSFER);
  return state;
}

inline void VolumeTransfer::CalcEn() {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_VOL_TRANSFER);
  if (GEMC_KIND == mv::GEMC_NVT) {
    if (isOrth) {
      cellList.GridAll(newDim, newMolsPos, molLookRef);
    } else {
      cellList.GridAll(newDimNonOrth, newMolsPos, molLookRef);
    }
  } else {
    if (isOrth) {
      cellList.GridBox(newDim, newMolsPos, molLookRef, box);
    } else {
      cellList.GridBox(newDimNonOrth, newMolsPos, molLookRef, box);
    }
  }

  regrewGrid = true;
  // back up cached Fourier term
  calcEwald->backupMolCache();
  sysPotNew = sysPotRef;

  if (GEMC_KIND == mv::GEMC_NVT) {
    for (uint b = 0; b < 2; b++) {
      // calculate new K vectors
      if (isOrth) {
        calcEwald->RecipInit(bPick[b], newDim);
        // setup reciprocal terms
        calcEwald->BoxReciprocalSetup(bPick[b], newMolsPos);
        sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newDim, bPick[b]);
      } else {
        calcEwald->RecipInit(bPick[b], newDimNonOrth);
        // setup reciprocal terms
        calcEwald->BoxReciprocalSetup(bPick[b], newMolsPos);
        sysPotNew =
            calcEnRef.BoxInter(sysPotNew, newMolsPos, newDimNonOrth, bPick[b]);
      }
      // calculate reciprocal term of electrostatic interaction
      sysPotNew.boxEnergy[bPick[b]].recip =
          calcEwald->BoxReciprocal(bPick[b], true);
    }
  } else {
    // calculate new K vectors
    if (isOrth) {
      calcEwald->RecipInit(box, newDim);
      // setup reciprocal terms
      calcEwald->BoxReciprocalSetup(box, newMolsPos);
      sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newDim, box);
    } else {
      calcEwald->RecipInit(box, newDimNonOrth);
      // setup reciprocal terms
      calcEwald->BoxReciprocalSetup(box, newMolsPos);
      sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newDimNonOrth, box);
    }
    // calculate reciprocal term of electrostatic interaction
    sysPotNew.boxEnergy[box].recip = calcEwald->BoxReciprocal(box, true);
  }

  sysPotNew.Total();
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_VOL_TRANSFER);
}

inline double VolumeTransfer::GetCoeff() const {
  ////Log-volume style shift -- is turned off, at present.
  //
  // return pow(newDim.volume[b_i]/boxDimRef.volume[b_i],
  //          (double)molLookRef.NumInBox(b_i)+1) *
  //  pow(newDim.volume[b_ii]/boxDimRef.volume[b_ii],
  //     (double)molLookRef.NumInBox(b_ii)+1);
  double coeff = 1.0;
  if (GEMC_KIND == mv::GEMC_NVT) {
    for (uint b = 0; b < 2; ++b) {
      if (isOrth)
        coeff *= pow(newDim.volume[bPick[b]] / boxDimRef.volume[bPick[b]],
                     (double)molLookRef.NumInBox(bPick[b]));
      else
        coeff *=
            pow(newDimNonOrth.volume[bPick[b]] / boxDimRef.volume[bPick[b]],
                (double)molLookRef.NumInBox(bPick[b]));
    }
  } else {
    if (isOrth) {
      coeff =
          pow(newDim.volume[box] / boxDimRef.volume[box],
              (double)molLookRef.NumInBox(box)) *
          exp(-BETA * PRESSURE * (newDim.volume[box] - boxDimRef.volume[box]));
    } else {
      coeff = pow(newDimNonOrth.volume[box] / boxDimRef.volume[box],
                  (double)molLookRef.NumInBox(box)) *
              exp(-BETA * PRESSURE *
                  (newDimNonOrth.volume[box] - boxDimRef.volume[box]));
    }
  }
  return coeff;
}

inline void VolumeTransfer::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_VOL_TRANSFER);
  double volTransCoeff = GetCoeff();
  double uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  double accept = volTransCoeff * uBoltz;
  bool result = (rejectState == mv::fail_state::NO_FAIL) && prng() < accept;
  if (result) {
    // Set new energy.
    sysPotRef = sysPotNew;
    // Swap... next time we'll use the current members.
    // NOTE:
    // This will be less efficient for NPT, but necessary evil.
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    if (isOrth)
      boxDimRef = newDim;
    else
      *(static_cast<BoxDimensionsNonOrth *>(&boxDimRef)) = newDimNonOrth;

    if (GEMC_KIND == mv::GEMC_NVT) {
      for (uint b = 0; b < 2; b++) {
        calcEwald->UpdateRecip(bPick[b]);
        calcEwald->UpdateRecipVec(bPick[b]);
      }
    } else {
      calcEwald->UpdateRecip(box);
      calcEwald->UpdateRecipVec(box);
    }
    // No need to update the velocity

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
      // update unitcell to the original in GPU
      for (uint b = 0; b < 2; b++) {
        UpdateCellBasisCUDA(forcefield.particles->getCUDAVars(), bPick[b],
                            boxDimRef.cellBasis[bPick[b]].x,
                            boxDimRef.cellBasis[bPick[b]].y,
                            boxDimRef.cellBasis[bPick[b]].z);
        if (!isOrth) {
          // In this case, boxDimRef is really an object of type
          // BoxDimensionsNonOrth,
          // so cast and copy the additional data to the GPU
          const BoxDimensionsNonOrth *NonOrthAxes =
              static_cast<const BoxDimensionsNonOrth *>(&boxDimRef);
          UpdateInvCellBasisCUDA(forcefield.particles->getCUDAVars(), bPick[b],
                                 NonOrthAxes->cellBasis_Inv[bPick[b]].x,
                                 NonOrthAxes->cellBasis_Inv[bPick[b]].y,
                                 NonOrthAxes->cellBasis_Inv[bPick[b]].z);
        }
      }
    } else {
      UpdateCellBasisCUDA(
          forcefield.particles->getCUDAVars(), box, boxDimRef.cellBasis[box].x,
          boxDimRef.cellBasis[box].y, boxDimRef.cellBasis[box].z);
      if (!isOrth) {
        // In this case, boxDimRef is really an object of type
        // BoxDimensionsNonOrth,
        // so cast and copy the additional data to the GPU
        const BoxDimensionsNonOrth *NonOrthAxes =
            static_cast<const BoxDimensionsNonOrth *>(&boxDimRef);
        UpdateInvCellBasisCUDA(forcefield.particles->getCUDAVars(), box,
                               NonOrthAxes->cellBasis_Inv[box].x,
                               NonOrthAxes->cellBasis_Inv[box].y,
                               NonOrthAxes->cellBasis_Inv[box].z);
      }
    }
#endif
  }

  if (GEMC_KIND == mv::GEMC_NVT) {
    moveSetRef.Update(mv::VOL_TRANSFER, result, bPick[0]);
    moveSetRef.Update(mv::VOL_TRANSFER, result, bPick[1]);
  } else {
    moveSetRef.Update(mv::VOL_TRANSFER, result, box);
  }
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_VOL_TRANSFER);
}

#endif

#endif /*VOLUMETRANSFER_H*/
