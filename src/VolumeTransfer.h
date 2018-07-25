/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
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

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  double GetCoeff() const;
  virtual void Accept(const uint rejectState, const uint step);
  virtual void PrintAcceptKind();
private:
  //Note: This is only used for GEMC-NPT
  uint bPick[BOX_TOTAL], subPick, subPickT[BOX_TOTAL];
  SystemPotential sysPotNew;
  const Forcefield& forcefield;
  BoxDimensions newDim;
  BoxDimensionsNonOrth newDimNonOrth;
  Coordinates newMolsPos;
  COM newCOMs;
  MoleculeLookup & molLookRef;
  const uint GEMC_KIND;
  const double PRESSURE;
  bool regrewGrid, isOrth;
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
}

void VolumeTransfer::PrintAcceptKind() {
  printf("%-37s", "% Accepted Volume-Transfer ");
  for(uint b = 0; b < BOX_TOTAL; b++) {
    uint index = mv::GetMoveSubIndex(mv::VOL_TRANSFER, b);
    printf("%10.5f ", 100.0 * moveSetRef.GetAccept(index));
  }
  std::cout << std::endl;
}

inline uint VolumeTransfer::Prep(const double subDraw, const double movPerc)
{
  uint state = mv::fail_state::NO_FAIL;
  if (GEMC_KIND == mv::GEMC_NVT) {
    subPick = mv::GetMoveSubIndex(mv::VOL_TRANSFER);
    for (uint b = 0; b < BOX_TOTAL; b++) {
      bPick[b] = b;
    }
  }
  if (GEMC_KIND == mv::GEMC_NPT) {
#if ENSEMBLE == NPT
    prng.PickBox(bPick[0], subDraw, movPerc);
#else
    prng.PickBoxPair(bPick[0], bPick[1], subDraw, movPerc);
#endif
    for (uint b = 0; b < BOX_TOTAL; b++) {
      subPickT[bPick[b]] = mv::GetMoveSubIndex(mv::VOL_TRANSFER, bPick[b]);
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
    double max = moveSetRef.Scale(subPick);
    if(isOrth)
      coordCurrRef.VolumeTransferTranslate(state, newMolsPos, newCOMs, newDim,
                                           comCurrRef, max);
    else
      coordCurrRef.VolumeTransferTranslate(state, newMolsPos, newCOMs,
                                           newDimNonOrth, comCurrRef, max);
  } else {
    XYZ scale[BOX_TOTAL];
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (state == mv::fail_state::NO_FAIL) {
        if ((bPick[b] == 0) && fixBox0)
          continue;

        double max = moveSetRef.Scale(subPickT[bPick[b]]);
        double delta = prng.Sym(max);
        if(isOrth)
          state =  boxDimRef.ShiftVolume(newDim, scale[bPick[b]],
                                         bPick[b], delta);
        else
          state =  boxDimRef.ShiftVolume(newDimNonOrth, scale[bPick[b]],
                                         bPick[b], delta);
      }
    }

    if (state == mv::fail_state::NO_FAIL) {
      for (uint b = 0; b < BOX_TOTAL; b++) {
        if ((bPick[b] == 0) && fixBox0)
          continue;

        if(isOrth)
          coordCurrRef.TranslateOneBox(newMolsPos, newCOMs, comCurrRef,
                                       newDim, bPick[b], scale[bPick[b]]);
        else
          coordCurrRef.TranslateOneBox(newMolsPos, newCOMs, comCurrRef,
                                       newDimNonOrth, bPick[b], scale[bPick[b]]);
      }
    }
  }
  return state;
}

inline void VolumeTransfer::CalcEn()
{
  if(isOrth)
    cellList.GridAll(newDim, newMolsPos, molLookRef);
  else
    cellList.GridAll(newDimNonOrth, newMolsPos, molLookRef);

  regrewGrid = true;

  //back up cached fourier term
  calcEwald->exgMolCache();
  sysPotNew = sysPotRef;
  for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
    if ((bPick[b] == 0) && fixBox0)
      continue;

    //calculate new K vectors
    if(isOrth)
      calcEwald->RecipInit(bPick[b], newDim);
    else
      calcEwald->RecipInit(bPick[b], newDimNonOrth);
    //setup reciprocate terms
    calcEwald->BoxReciprocalSetup(bPick[b], newMolsPos);
    //calculate LJ interaction and real term of electrostatic interaction
    if(isOrth)
      sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newCOMs, newDim,
                                     bPick[b]);
    else
      sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolsPos, newCOMs,
                                     newDimNonOrth, bPick[b]);
    //calculate reciprocate term of electrostatic interaction
    sysPotNew.boxEnergy[bPick[b]].recip = calcEwald->BoxReciprocal(bPick[b]);
  }
  sysPotNew.Total();
}



inline double VolumeTransfer::GetCoeff() const
{
  ////Log-volume style shift -- is turned off, at present.
  //
  //return pow(newDim.volume[b_i]/boxDimRef.volume[b_i],
  //	      (double)molLookRef.NumInBox(b_i)+1) *
  //  pow(newDim.volume[b_ii]/boxDimRef.volume[b_ii],
  //	 (double)molLookRef.NumInBox(b_ii)+1);
  double coeff = 1.0;
  if (GEMC_KIND == mv::GEMC_NVT) {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      if(isOrth)
        coeff *= pow(newDim.volume[b] / boxDimRef.volume[b],
                     (double)molLookRef.NumInBox(b));
      else
        coeff *= pow(newDimNonOrth.volume[b] / boxDimRef.volume[b],
                     (double)molLookRef.NumInBox(b));
    }
  } else {
    for (uint b = 0; b < BOX_TOTAL; ++b) {
      if ((bPick[b] == 0) && fixBox0)
        continue;

      if(isOrth)
        coeff *= pow(newDim.volume[b] / boxDimRef.volume[b],
                     (double)molLookRef.NumInBox(b)) *
                 exp(-BETA * PRESSURE * (newDim.volume[b] - boxDimRef.volume[b]));
      else
        coeff *= pow(newDimNonOrth.volume[b] / boxDimRef.volume[b],
                     (double)molLookRef.NumInBox(b)) *
                 exp(-BETA * PRESSURE * (newDimNonOrth.volume[b] - boxDimRef.volume[b]));
    }

  }
  return coeff;
}

inline void VolumeTransfer::Accept(const uint rejectState, const uint step)
{
  double volTransCoeff = GetCoeff();
  double uBoltz = exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));
  double accept = volTransCoeff * uBoltz;
  bool result = (rejectState == mv::fail_state::NO_FAIL) && prng() < accept;
  if (result) {
    //Set new energy.
    sysPotRef = sysPotNew;
    //Swap... next time we'll use the current members.
    //NOTE:
    //This will be less efficient for NPT, but necessary evil.
    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);
    if(isOrth)
      boxDimRef = newDim;
    else
      *((BoxDimensionsNonOrth*)(&boxDimRef)) = newDimNonOrth;

    for (uint b = 0; b < BOX_TOTAL; b++) {
      if ((bPick[b] == 0) && fixBox0)
        continue;

      calcEwald->UpdateRecip(b);
      calcEwald->UpdateRecipVec(b);
    }
  } else if (rejectState == mv::fail_state::NO_FAIL && regrewGrid) {
    cellList.GridAll(boxDimRef, coordCurrRef, molLookRef);
    regrewGrid = false;
    calcEwald->exgMolCache();
#ifdef GOMC_CUDA
    //update unitcell to the original in GPU
    for (uint box = 0; box < BOX_TOTAL; box++) {

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
    subPick = mv::GetMoveSubIndex(mv::VOL_TRANSFER, 0);
    moveSetRef.Update(result, subPick, step);
    subPick = mv::GetMoveSubIndex(mv::VOL_TRANSFER, 1);
    moveSetRef.Update(result, subPick, step);
  }
  if (GEMC_KIND == mv::GEMC_NPT) {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if ((bPick[b] == 0) && fixBox0)
        continue;

      subPickT[bPick[b]] = mv::GetMoveSubIndex(mv::VOL_TRANSFER, bPick[b]);
      moveSetRef.Update(result, subPickT[bPick[b]], step);
    }
  }

}

#endif

#endif /*VOLUMETRANSFER_H*/
