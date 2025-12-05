/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef NEMTMC_H
#define NEMTMC_H

#if ENSEMBLE == GCMC || ENSEMBLE == GEMC

#include <numeric>

#include "MoveBase.h"
#include "TrialMol.h"

using std::vector;

class NEMTMC : public MoveBase {
public:
  NEMTMC(System &sys, StaticVals const &statV)
      : MoveBase(sys, statV), lambdaRef(sys.lambdaRef), systemRef(sys),
        ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
        backUpCoordinate(boxDimRef, backUpCOM, molLookRef, sys.prng, statV.mol),
        backUpCOM(boxDimRef, backUpCoordinate, molLookRef, statV.mol),
        backUpMoveSetting(sys.boxDimRef), propagationMove(NULL),
        r123wrapper(sys.r123wrapper) {
    if (statV.neMTMCVal.enable) {
      conformationProb = statV.neMTMCVal.conformationProb;
      lambdaLimit = statV.neMTMCVal.lambdaLimit;
      backUpCoordinate.Init(sys.coordinates.Count());
      backUpCOM.Init(statV.mol.count);

      MPEnable = statV.neMTMCVal.MPEnable;
      BrownianDynamicEnable = statV.neMTMCVal.MPBEnable;
      stepCounter = 0;
      relaxSteps = statV.neMTMCVal.relaxSteps;
      lambdaWindow = statV.neMTMCVal.lambdaVDW.size() - 1;
      lambdaCoulomb = statV.neMTMCVal.lambdaCoulomb;
      lambdaVDW = statV.neMTMCVal.lambdaVDW;
      uint totKind = molRef.GetKindsCount();
      kCount.resize(BOX_TOTAL);
      for (uint b = 0; b < BOX_TOTAL; b++) {
        kCount[b].resize(totKind);
      }
    }
  }

  virtual uint Prep(const double subDraw, const double movPerc);
  // To relax the system in NE_MTMC move
  virtual uint PrepNEMTMC(const uint box, const uint midx = 0,
                          const uint kidx = 0) {
    return mv::fail_state::NO_FAIL;
  }
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const ulong step);
  virtual void PrintAcceptKind();

private:
  double GetCoeff() const;
  uint GetBoxPairAndMol(const double subDraw, const double movPerc);
  void ShiftMolToSourceBox();
  void ShiftMolToDestBox();
  void UpdateEnergy();
  void CalcEnNEMT(uint lambdaIdxOldS, uint lambdaIdxNewS);

  void Propagation();
  void RelaxingTransform(uint box);
  void AddWork();

  Lambda &lambdaRef;
  System &systemRef;
  Random123Wrapper &r123wrapper;
  uint sourceBox, destBox;
  uint pStartNEMT, pLenNEMT;
  uint molIndex, kindIndex;
  int lambdaIdxOld, lambdaIdxNew;
  uint relaxSteps, lambdaWindow;
  vector<vector<uint>> kCount;

  double work;
  double conformationProb, lambdaLimit;
  double W_tc, W_recip;
  double correctDiffSource, correctDiffDest, selfDiffSource, selfDiffDest;
  double recipDiffSource, recipDiffDest, tcDiffSource, tcDiffDest;
  double molInSourceBox, molInDestBox;
  vector<double> lambdaCoulomb, lambdaVDW;

  // variable needs for relaxing
  bool MPEnable, BrownianDynamicEnable;
  ulong stepCounter;
  MoveBase *propagationMove;

  cbmc::TrialMol oldMolNEMT, newMolNEMT;
  Energy oldEnergy[BOX_TOTAL], newEnergy[BOX_TOTAL];
  Forcefield const &ffRef;
  MoleculeLookup &molLookRef;
  SystemPotential sysPotNew;
  // Backup data
  SystemPotential backUpPotential;
  Coordinates backUpCoordinate;
  COM backUpCOM;
  MoveSettings backUpMoveSetting;
};

void NEMTMC::PrintAcceptKind() {
  for (uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted nonEq Mol-Transfer ",
           molRef.kinds[k].name.c_str());
    for (uint b = 0; b < BOX_TOTAL; b++) {
      if (moveSetRef.GetTrial(b, mv::NE_MTMC, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::NE_MTMC, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint NEMTMC::GetBoxPairAndMol(const double subDraw,
                                     const double movPerc) {
  // Need to call a function to pick a molecule that is not fixed but cannot be
  // swapped between boxes. (beta != 1, beta !=2)
  uint state = prng.PickMolAndBoxPair2(molIndex, kindIndex, sourceBox, destBox,
                                       subDraw, movPerc);
#if ENSEMBLE == GCMC
  if (state == mv::fail_state::NO_MOL_OF_KIND_IN_BOX && sourceBox == mv::BOX1) {
    std::cout << "Error: There are no molecules of kind "
              << molRef.kinds[kindIndex].name << " left in reservoir.\n";
    exit(EXIT_FAILURE);
  }
#endif

  if (state == mv::fail_state::NO_FAIL) {
    pStartNEMT = pLenNEMT = 0;
    molRef.GetRangeStartLength(pStartNEMT, pLenNEMT, molIndex);
  }
  return state;
}

inline uint NEMTMC::Prep(const double subDraw, const double movPerc) {
  GOMC_EVENT_START(1, GomcProfileEvent::PREP_NEMTMC);
  // back up the stats and data that needs to be reset if move gets rejected
  coordCurrRef.CopyRange(backUpCoordinate, 0, 0, coordCurrRef.Count());
  comCurrRef.CopyRange(backUpCOM, 0, 0, comCurrRef.Count());
  backUpPotential = sysPotRef;
  backUpMoveSetting.SetValues(moveSetRef);
  stepCounter = 0;
  work = 0.0;

  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    newMolNEMT = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMolNEMT = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMolNEMT.SetCoords(coordCurrRef, pStartNEMT);
    // Unwrap the old coordinates for use with new coordinates after wrapping
    XYZArray mol(pLenNEMT);
    coordCurrRef.CopyRange(mol, pStartNEMT, 0, pLenNEMT);
    boxDimRef.UnwrapPBC(mol, sourceBox, comCurrRef.Get(molIndex));
    boxDimRef.WrapPBC(mol, destBox);
    // Later it will shift to random COM
    newMolNEMT.SetCoords(mol, 0);
    // store number of molecules in box before shifting molecule
    molInSourceBox =
        (double)molLookRef.NumKindInBoxSwappable(kindIndex, sourceBox);
    molInDestBox = (double)molLookRef.NumKindInBoxSwappable(kindIndex, destBox);
    for (uint b = 0; b < BOX_TOTAL; b++) {
      for (uint k = 0; k < molRef.GetKindsCount(); ++k) {
        kCount[b][k] = molLookRef.NumKindInBox(k, b);
        if (b == sourceBox && k == kindIndex) {
          // consider the fraction molecule as a different molecule kind
          --kCount[b][k];
        }
      }
    }
  }

  GOMC_EVENT_STOP(1, GomcProfileEvent::PREP_NEMTMC);
  return state;
}

inline uint NEMTMC::Transform() {
  GOMC_EVENT_START(1, GomcProfileEvent::TRANS_NEMTMC);
  // Start with full interaction in sourceBox, zero interaction in destBox
  // Since we have the lambda for growing molecule, in sourceBox, lambdaWindow
  // index corresponds to full interaction (lambda = 1) and (lambdaWindow - X)
  // is for destBox
  lambdaIdxOld = lambdaWindow;
  lambdaIdxNew = lambdaIdxOld - 1;
  // Like insert in vacuum, we use this just to create a random configuration
  lambdaRef.Set(0.0, 0.0, molIndex, kindIndex, sourceBox);
  lambdaRef.Set(0.0, 0.0, molIndex, kindIndex, destBox);
  // Start growing the fractional molecule in destBox
  molRef.kinds[kindIndex].BuildIDNew(newMolNEMT, molIndex);
  ShiftMolToDestBox();
  molRef.kinds[kindIndex].BuildIDOld(oldMolNEMT, molIndex);

  Energy originalBond = calcEnRef.MoleculeIntra(newMolNEMT);

  while (lambdaIdxNew >= 0) {
    // Set the interaction in sourceBox and destBox
    lambdaRef.Set(lambdaVDW[lambdaIdxNew], lambdaCoulomb[lambdaIdxNew],
                  molIndex, kindIndex, sourceBox);
    lambdaRef.Set(lambdaVDW[lambdaWindow - lambdaIdxNew],
                  lambdaCoulomb[lambdaWindow - lambdaIdxNew], molIndex,
                  kindIndex, destBox);
    // Calculate the old and new energy in sourceBox and destBox
    CalcEnNEMT(lambdaIdxOld, lambdaIdxNew);
    // Calculate the nonEq work
    AddWork();
    // Update the system Energy
    UpdateEnergy();
    // Relax the system (propagation)
    Propagation();
    // Go to the next lambda
    --lambdaIdxNew;
    --lambdaIdxOld;
  };
  // We have not transferred the bonded energy. We apply it at the end. Note: We
  // might call Regrowth, so we are considering the newest conformation.
  sysPotRef.boxEnergy[sourceBox] -= calcEnRef.MoleculeIntra(oldMolNEMT);
  sysPotRef.boxEnergy[destBox] += originalBond;
  sysPotRef.Total();

  GOMC_EVENT_STOP(1, GomcProfileEvent::TRANS_NEMTMC);
  return mv::fail_state::NO_FAIL;
}

inline void NEMTMC::CalcEn() { return; }

inline void NEMTMC::CalcEnNEMT(uint lambdaIdxOldS, uint lambdaIdxNewS) {
  GOMC_EVENT_START(1, GomcProfileEvent::CALC_EN_NEMTMC);
  W_tc = 1.0;
  W_recip = 1.0;
  correctDiffDest = correctDiffSource = 0.0;
  selfDiffDest = selfDiffSource = 0.0;
  recipDiffDest = recipDiffSource = 0.0;
  tcDiffDest = tcDiffSource = 0.0;
  oldEnergy[sourceBox].Zero();
  newEnergy[sourceBox].Zero();
  oldEnergy[destBox].Zero();
  newEnergy[destBox].Zero();

  double lambdaOld_VDW_S = lambdaVDW[lambdaIdxOldS];
  double lambdaNew_VDW_S = lambdaVDW[lambdaIdxNewS];
  double lambdaOld_VDW_D = lambdaVDW[lambdaWindow - lambdaIdxOldS];
  double lambdaNew_VDW_D = lambdaVDW[lambdaWindow - lambdaIdxNewS];
  double lambdaOld_Coulomb_S = lambdaCoulomb[lambdaIdxOldS];
  double lambdaNew_Coulomb_S = lambdaCoulomb[lambdaIdxNewS];
  double lambdaOld_Coulomb_D = lambdaCoulomb[lambdaWindow - lambdaIdxOldS];
  double lambdaNew_Coulomb_D = lambdaCoulomb[lambdaWindow - lambdaIdxNewS];

  // Calculating long range correction
  if (ffRef.useLRC) {
    // Calculate LRC difference for lambdaNew and lambdaOld
    tcDiffSource =
        calcEnRef.MoleculeTailChange(sourceBox, kindIndex, kCount[sourceBox],
                                     lambdaOld_VDW_S, lambdaNew_VDW_S);
    tcDiffDest = calcEnRef.MoleculeTailChange(
        destBox, kindIndex, kCount[destBox], lambdaOld_VDW_D, lambdaNew_VDW_D);
    W_tc = ((tcDiffSource + tcDiffDest));
  }

  ShiftMolToSourceBox();
  // calculate inter energy for lambda new and old in source Box
  calcEnRef.SingleMoleculeInter(oldEnergy[sourceBox], newEnergy[sourceBox],
                                lambdaOld_VDW_S, lambdaNew_VDW_S,
                                lambdaOld_Coulomb_S, lambdaNew_Coulomb_S,
                                molIndex, sourceBox);

  ShiftMolToDestBox();
  // calculate inter energy for lambda new and old in dest Box
  calcEnRef.SingleMoleculeInter(
      oldEnergy[destBox], newEnergy[destBox], lambdaOld_VDW_D, lambdaNew_VDW_D,
      lambdaOld_Coulomb_D, lambdaNew_Coulomb_D, molIndex, destBox);

  // Calculate self and correction difference for lambdaNew and lambdaOld
  // For electrostatic we use linear scaling
  double coefDiffS = lambdaNew_Coulomb_S - lambdaOld_Coulomb_S;
  double coefDiffD = lambdaNew_Coulomb_D - lambdaOld_Coulomb_D;
  correctDiffSource = coefDiffS * calcEwald->SwapCorrection(oldMolNEMT);
  correctDiffDest = coefDiffD * calcEwald->SwapCorrection(newMolNEMT);
  selfDiffSource = coefDiffS * calcEwald->SwapSelf(oldMolNEMT);
  selfDiffDest = coefDiffD * calcEwald->SwapSelf(newMolNEMT);
  // calculate Reciprocal Difference in source and dest box
  recipDiffSource =
      calcEwald->ChangeLambdaRecip(oldMolNEMT.GetCoords(), lambdaOld_Coulomb_S,
                                   lambdaNew_Coulomb_S, molIndex, sourceBox);
  recipDiffDest =
      calcEwald->ChangeLambdaRecip(newMolNEMT.GetCoords(), lambdaOld_Coulomb_D,
                                   lambdaNew_Coulomb_D, molIndex, destBox);

  // need to contribute the self and correction energy
  W_recip = ((recipDiffSource + recipDiffDest + correctDiffSource +
              correctDiffDest + selfDiffSource + selfDiffDest));
  GOMC_EVENT_STOP(1, GomcProfileEvent::CALC_EN_NEMTMC);
}

inline double NEMTMC::GetCoeff() const {
#if ENSEMBLE == GEMC
  return (molInSourceBox) / (molInDestBox + 1.0) * boxDimRef.volume[destBox] *
         boxDimRef.volInv[sourceBox];
#elif ENSEMBLE == GCMC
  if (sourceBox == mv::BOX0) { // Deletion case
    if (ffRef.isFugacity) {
      return (molInSourceBox)*boxDimRef.volInv[sourceBox] /
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return (molInSourceBox)*boxDimRef.volInv[sourceBox] *
             exp(-BETA * molRef.kinds[kindIndex].chemPot);
    }
  } else { // Insertion case
    if (ffRef.isFugacity) {
      return boxDimRef.volume[destBox] / (molInDestBox + 1.0) *
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return boxDimRef.volume[destBox] / (molInDestBox + 1.0) *
             exp(BETA * molRef.kinds[kindIndex].chemPot);
    }
  }
#endif
}

inline void NEMTMC::Accept(const uint rejectState, const ulong step) {
  GOMC_EVENT_START(1, GomcProfileEvent::ACC_NEMTMC);
  double molTransCoeff = GetCoeff();
  bool result = prng() < (molTransCoeff * exp(-BETA * work));
  // If we didn't skip the move calculation
  if (rejectState == mv::fail_state::NO_FAIL) {
    if (result) {
      // Set full interaction in destBox, zero interaction in sourceBox
      lambdaRef.UnSet(destBox, sourceBox);
    } else {
      // Set full interaction in sourceBox, zero interaction in destBox
      lambdaRef.UnSet(sourceBox, destBox);
      // Shift the molecule back
      ShiftMolToSourceBox();
      backUpCoordinate.CopyRange(coordCurrRef, 0, 0, coordCurrRef.Count());
      backUpCOM.CopyRange(comCurrRef, 0, 0, comCurrRef.Count());
      cellList.GridAll(boxDimRef, coordCurrRef, molLookRef);
      for (uint b = 0; b < BOXES_WITH_U_NB; ++b) {
        calcEwald->BoxReciprocalSums(b, coordCurrRef);
        calcEwald->UpdateRecip(b);
      }
      sysPotRef = backUpPotential;
      // We already messed up the forces, need to recalculate it
      backUpMoveSetting.SetSingleMoveAccepted(sourceBox);
      backUpMoveSetting.SetSingleMoveAccepted(destBox);
      // reset the trial/acceptance to original values if the move is rejected
      // moveSetRef.SetValues(backUpMoveSetting);
    }
  }
  // reset the key, because next move has different step number uk[0] = step!
  r123wrapper.SetKey(0);
  // reset the trial/acceptance to original values if the move is rejected
  moveSetRef.SetValues(backUpMoveSetting);
  moveSetRef.Update(mv::NE_MTMC, result, destBox, kindIndex);
  GOMC_EVENT_STOP(1, GomcProfileEvent::ACC_NEMTMC);
}

inline void NEMTMC::ShiftMolToSourceBox() {
  cellList.RemoveMol(molIndex, destBox, coordCurrRef);
  // Set coordinates, new COM; shift index to new box's list
  oldMolNEMT.GetCoords().CopyRange(coordCurrRef, 0, pStartNEMT, pLenNEMT);
  comCurrRef.SetNew(molIndex, sourceBox);
  molLookRef.ShiftMolBox(molIndex, destBox, sourceBox, kindIndex);
  cellList.AddMol(molIndex, sourceBox, coordCurrRef);
}

inline void NEMTMC::ShiftMolToDestBox() {
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  // Set coordinates, new COM; shift index to new box's list
  newMolNEMT.GetCoords().CopyRange(coordCurrRef, 0, pStartNEMT, pLenNEMT);
  comCurrRef.SetNew(molIndex, destBox);
  molLookRef.ShiftMolBox(molIndex, sourceBox, destBox, kindIndex);
  cellList.AddMol(molIndex, destBox, coordCurrRef);
}

void NEMTMC::AddWork() {
  double W1 = newEnergy[sourceBox].Total() - oldEnergy[sourceBox].Total();
  double W2 = newEnergy[destBox].Total() - oldEnergy[destBox].Total();
  work += (W1 + W2 + W_tc + W_recip);
}

inline void NEMTMC::UpdateEnergy() {
  // Add tail corrections
  sysPotRef.boxEnergy[sourceBox].tailCorrection += tcDiffSource;
  sysPotRef.boxEnergy[destBox].tailCorrection += tcDiffDest;
  // Add rest of energy.
  sysPotRef.boxEnergy[sourceBox] -= oldEnergy[sourceBox];
  sysPotRef.boxEnergy[sourceBox] += newEnergy[sourceBox];
  sysPotRef.boxEnergy[destBox] -= oldEnergy[destBox];
  sysPotRef.boxEnergy[destBox] += newEnergy[destBox];
  // Add correction energy
  sysPotRef.boxEnergy[sourceBox].correction += correctDiffSource;
  sysPotRef.boxEnergy[destBox].correction += correctDiffDest;
  // Add self energy
  sysPotRef.boxEnergy[sourceBox].self += selfDiffSource;
  sysPotRef.boxEnergy[destBox].self += selfDiffDest;
  // Add Reciprocal energy
  sysPotRef.boxEnergy[sourceBox].recip += recipDiffSource;
  sysPotRef.boxEnergy[destBox].recip += recipDiffDest;

  calcEwald->UpdateRecip(sourceBox);
  calcEwald->UpdateRecip(destBox);
  moveSetRef.SetSingleMoveAccepted(sourceBox);
  moveSetRef.SetSingleMoveAccepted(destBox);

  // Retotal
  sysPotRef.Total();
}

inline void NEMTMC::Propagation() {
  // Relax the source Box
  ShiftMolToSourceBox();
  if (sourceBox < BOXES_WITH_U_NB) {
    RelaxingTransform(sourceBox);
    oldMolNEMT.SetCoords(coordCurrRef, pStartNEMT);
  }

  // Relax the destination Box
  ShiftMolToDestBox();
  if (destBox < BOXES_WITH_U_NB) {
    RelaxingTransform(destBox);
    newMolNEMT.SetCoords(coordCurrRef, pStartNEMT);
  }
}

inline void NEMTMC::RelaxingTransform(uint box) {
  for (uint s = 0; s < relaxSteps; s++) {
    uint rejectState = mv::fail_state::NO_FAIL;
    // With conformationProb probability, we perform sample conformation
    // using CD-CBMC
    bool sampleConf = (prng() <= conformationProb);

    if (sampleConf) {
      // sample conformation
      // Check lambda to decide to perform Regrowth or IntraSwap
      if (lambdaRef.GetLambdaVDW(molIndex, box) <= lambdaLimit) {
        // Perform IntraSwap or Regrowth move
        propagationMove = systemRef.GetMoveObject(
            ((prng() < 0.5) ? mv::INTRA_SWAP : mv::REGROWTH));
      } else {
        // Perform Regrowth move
        propagationMove = systemRef.GetMoveObject(mv::REGROWTH);
      }
      rejectState = propagationMove->PrepNEMTMC(box, molIndex, kindIndex);

    } else if (MPEnable || BrownianDynamicEnable) {
      // Relax the system using Brownian Motion or Force-Based MP move
      // Change the key number, otherwise we will perform the same move!
      // stepCounter resets every time in Prep()
      r123wrapper.SetKey(s + stepCounter);
      // Use multiparticle/brownian dynamic to propagate
      if (BrownianDynamicEnable) {
        propagationMove = systemRef.GetMoveObject(mv::MULTIPARTICLE_BM);
      } else {
        propagationMove = systemRef.GetMoveObject(mv::MULTIPARTICLE);
      }
      // prepare the move with box
      rejectState = propagationMove->PrepNEMTMC(box);

    } else {
      // Use random translation and rotation to propagate
      // Randomly pick a molecule in Box
      uint pStart = 0;
      uint pLen = 0;
      uint m = 0;
      uint mk = 0;
      rejectState = prng.PickMol(m, mk, box);
      if (rejectState == mv::fail_state::NO_FAIL) {
        molRef.GetRangeStartLength(pStart, pLen, m);
        if (pLen == 1) {
          // We do displacement if we have single site atom
          propagationMove = systemRef.GetMoveObject(mv::DISPLACE);
        } else {
          // get the displace/rotate move to propagate with 50% probability
          propagationMove = systemRef.GetMoveObject(
              ((prng() < 0.5) ? mv::ROTATE : mv::DISPLACE));
        }
        // prepare the move with box, picked molkind and index
        rejectState = propagationMove->PrepNEMTMC(box, m, mk);
      }
    }

    // Transform, CalcEn, and Accept/Reject
    if (rejectState == mv::fail_state::NO_FAIL) {
      rejectState = propagationMove->Transform();
    }
    if (rejectState == mv::fail_state::NO_FAIL) {
      propagationMove->CalcEn();
    }
    propagationMove->Accept(rejectState, s + stepCounter);
    // moveSetRef.AdjustMoves(s + stepCounter);
  }
  // To advance the number of steps during intermediate states
  stepCounter += relaxSteps;
}

#endif
#endif
