#ifndef INTRAMOLECULEEXCHANGE1_H
#define INTRAMOLECULEEXCHANGE1_H

#include "MoveBase.h"
#include "TrialMol.h"
#include "GeomLib.h"
#include <cmath>

using std::vector;
using namespace geom;

// Intra Molecule Exchange Move:
// KindA is small kind. KindB is large kind
// Orientation center of cavA is random.
// Orientation center of cavB is on COM of kindL, aligned with kindL backbon.
// Delete the exchangeRatio kindS from cavA, and 1 kindL from cavB.
// Insert the exchangeRatio kindS to cavB and 1 kindL inside the cavA.

class IntraMoleculeExchange1 : public MoveBase
{
public:

  IntraMoleculeExchange1(System &sys, StaticVals const& statV) :
    ffRef(statV.forcefield), molLookRef(sys.molLookupRef), MoveBase(sys, statV),
    cavity(statV.intraMemcVal.subVol), cavA(3), invCavA(3),
    perAdjust(statV.GetPerAdjust()), cavB(3), invCavB(3)
  {
    enableID = statV.intraMemcVal.enable;
    trial.resize(BOX_TOTAL);
    accepted.resize(BOX_TOTAL);

    if(enableID) {
      if(molLookRef.GetNumCanSwapKind() < 2) {
        std::cout << "Error: MEMC move cannot be applied to pure systems or" <<
                  " systems, where only one molecule type is allowed to be swapped.\n";
        exit(EXIT_FAILURE);
      }

      if(cavity.x >= cavity.y)
        cavity.y = cavity.x;
      else
        cavity.x = cavity.y;

      volCav = cavity.x * cavity.y * cavity.z;
      exchangeRatioVec = statV.intraMemcVal.exchangeRatio;

      SetMEMC(statV);

      //checking the acceptance statistic for each kind
      counter = 0;
      molInCavCount = 0;
      lastAccept = 0.0;
      exDiff = 1;
      for(uint b = 0; b < BOX_TOTAL; b++) {
        trial[b].resize(molRef.GetKindsCount() * molRef.GetKindsCount(), 0.0);
        accepted[b].resize(molRef.GetKindsCount() * molRef.GetKindsCount(), 0.0);
      }
    }
  }

  virtual uint Prep(const real subDraw, const real movPerc);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const uint step);
  virtual void PrintAcceptKind();
  //This function carries out actions based on the internal acceptance state and
  //molecule kind
  void AcceptKind(const uint rejectState, const uint kind, const uint box);

protected:

  virtual void AdjustExRatio();
  virtual void SetMEMC(StaticVals const& statV);
  virtual void SetExchangeData();
  void ShiftMol(const uint n, const bool kindA);
  void RecoverMol(const uint n, const bool kindA);
  virtual uint PickMolInCav();
  virtual real GetCoeff() const;
  uint GetBoxPairAndMol(const real subDraw, const real movPerc);

  bool enableID;
  uint sourceBox;
  uint largeBB[2];
  uint perAdjust, molInCavCount, counter;
  uint numInCavA, numInCavB, numSCavA, numSCavB, kindS, kindL;
  vector<uint> pStartA, pLenA, pStartB, pLenB;
  vector<uint> molIndexA, kindIndexA, molIndexB, kindIndexB;
  vector<cbmc::TrialMol> oldMolA, newMolA, oldMolB, newMolB;
  vector< vector<uint> > molInCav;
  //To store total sets of exchange pairs
  vector<uint> exchangeRatioVec, kindSVec, kindLVec;
  vector< vector<uint> > largeBBVec;
  //For move acceptance of each molecule kind
  std::vector< std::vector<uint> > trial, accepted;

  int exDiff, exchangeRatio;
  real volCav, lastAccept;
  real W_recip, recipDiffA, recipDiffB, correctDiff;

  XYZ centerA, centerB, cavity;
  XYZArray cavA, invCavA, cavB, invCavB;

  MoleculeLookup & molLookRef;
  Forcefield const& ffRef;
};

inline void IntraMoleculeExchange1::AcceptKind(const uint rejectState, const uint kind,
    const uint box)
{
  trial[box][kind]++;
  if(rejectState)
    accepted[box][kind]++;
}

void IntraMoleculeExchange1::PrintAcceptKind()
{
  for(uint k = 0; k < kindLVec.size(); k++) {
    uint ks = kindSVec[k];
    uint kl = kindLVec[k];
    uint index = ks * molRef.GetKindsCount() + kl;
    printf("%-22s %5s - %-5s ", "% Accepted Intra-MEMC ", molRef.kinds[kl].name.c_str(),
           molRef.kinds[ks].name.c_str());
    for(uint b = 0; b < BOX_TOTAL; b++) {
      if(trial[b][index] > 0)
        printf("%10.5f ", (real)(100.0 * accepted[b][index] / trial[b][index]));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline void IntraMoleculeExchange1::SetMEMC(StaticVals const& statV)
{
  for(uint t = 0; t < exchangeRatioVec.size(); t++) {
    kindS = kindL = largeBB[0] = largeBB[1] = -1;
    for(uint k = 0; k < molLookRef.GetNumCanMoveKind(); k++) {
      uint kind = molLookRef.GetCanMoveKind(k);
      if(molRef.kinds[kind].name == statV.intraMemcVal.largeKind[t]) {
        kindL = kind;
      } else if(molRef.kinds[kind].name == statV.intraMemcVal.smallKind[t]) {
        kindS = kind;
      }
    }

    if(kindS == -1) {
      printf("Error: Residue name %s was not found in PDB file as small molecule kind to be exchanged or it is fixed in its position.\n",
             statV.intraMemcVal.smallKind[t].c_str());
      exit(EXIT_FAILURE);
    }

    if(kindL == -1) {
      printf("Error: Residue name %s was not found in PDB file as large molecule kind to be exchanged or it is fixed in its position.\n",
             statV.intraMemcVal.largeKind[t].c_str());
      exit(EXIT_FAILURE);
    }

    for(uint i = 0; i < molRef.kinds[kindL].NumAtoms(); i++) {
      if(molRef.kinds[kindL].atomNames[i] == statV.intraMemcVal.largeBBAtom1[t]) {
        largeBB[0] = i;
      }
      if(molRef.kinds[kindL].atomNames[i] == statV.intraMemcVal.largeBBAtom2[t]) {
        largeBB[1] = i;
      }
    }

    for(uint i = 0; i < 2; i++) {
      if(largeBB[i] == -1) {
        printf("Error: Atom name %s or %s was not found in %s residue.\n",
               statV.intraMemcVal.largeBBAtom1[t].c_str(),
               statV.intraMemcVal.largeBBAtom2[t].c_str(),
               statV.intraMemcVal.largeKind[t].c_str());
        exit(EXIT_FAILURE);
      }
    }

    if(statV.intraMemcVal.MEMC1 || statV.intraMemcVal.MEMC2) {
      if(molRef.kinds[kindL].NumAtoms() > 1) {
        if(largeBB[0] == largeBB[1]) {
          printf("Error: Atom names in large molecule backbone cannot be same!\n");
          exit(EXIT_FAILURE);
        }
      }
    }
    kindSVec.push_back(kindS);
    kindLVec.push_back(kindL);
    vector<uint> temp(largeBB, largeBB + 2);
    largeBBVec.push_back(temp);
  }
}

inline void IntraMoleculeExchange1::AdjustExRatio()
{
  if(((counter + 1) % perAdjust) == 0) {
    uint exMax = ceil((float)molInCavCount / (float)perAdjust);
    uint exMin = ceil((float)exMax / 2.0);
    if(exMin == 0)
      exMin = 1;

    uint index = kindS + kindL * molRef.GetKindsCount();
    real currAccept = (real)(accepted[sourceBox][index]) / (real)(trial[sourceBox][index]);
    if(abs(currAccept - lastAccept) >= 0.05 * currAccept) {
      if(currAccept > lastAccept) {
        exchangeRatio += exDiff;
      } else {
        exDiff *= -1;
        exchangeRatio += exDiff;
      }
      lastAccept = currAccept;
      if(exchangeRatio < exMin)
        exchangeRatio  = exMin;
      if(exchangeRatio > exMax)
        exchangeRatio = exMax;
    }
    molInCavCount = 0;
    counter = 0;
    printf("Average Mol In Cavity: %d. Exchange Ratio: %d \n", exMax,
           exchangeRatio);
  }
}

inline void IntraMoleculeExchange1::SetExchangeData()
{
  uint exType = prng.randIntExc(exchangeRatioVec.size());
  kindS = kindSVec[exType];
  kindL = kindLVec[exType];
  exchangeRatio = exchangeRatioVec[exType];
  largeBB[0] = largeBBVec[exType][0];
  largeBB[1] = largeBBVec[exType][1];
}

inline uint IntraMoleculeExchange1::PickMolInCav()
{
  uint state = mv::fail_state::NO_FAIL;
  //pick a random location in source box
  XYZ axis = boxDimRef.GetAxis(sourceBox);
  XYZ temp(prng.randExc(axis.x), prng.randExc(axis.y), prng.randExc(axis.z));
  //Use to find small molecule
  centerA = temp;
  //Pick random vector and find two vectors that are perpendicular to V1
  SetBasis(cavA, prng.RandomUnitVect());
  //Calculate inverse matrix for cav, here Inv = transpose
  TransposeMatrix(invCavA, cavA);
  //Find the small molecule kind in the cavityA
  if(calcEnRef.FindMolInCavity(molInCav, centerA, cavity, invCavA,
                               sourceBox, kindS, exchangeRatio)) {
    //printf("MolS in cav: %d.\n", molInCav[kindS].size());
    //Find the exchangeRatio number of molecules kindS in cavity
    numInCavA = exchangeRatio;
    numSCavA = molInCav[kindS].size();
    for(uint n = 0; n < numInCavA; n++) {
      //pick random exchangeRatio number of kindS in cavity
      uint picked = prng.randIntExc(molInCav[kindS].size());
      molIndexA.push_back(molInCav[kindS][picked]);
      kindIndexA.push_back(molRef.GetMolKind(molIndexA[n]));
      molInCav[kindS].erase(molInCav[kindS].begin() + picked);
    }

    //pick a molecule from kindL in the same box
    numInCavB = 1;
    state = prng.PickMol(kindL, kindIndexB, molIndexB, numInCavB, sourceBox);
  } else {
    //reject the move
    state = mv::fail_state::NO_MOL_OF_KIND_IN_BOX;
  }

  //keep track of number of small molecule in cavity
  molInCavCount += numSCavA;
  counter++;

  //After picking a large molecule, set the cavityB, and count kindS in cavityB
  if(state == mv::fail_state::NO_FAIL) {
    //Use to insert kindS in this centerB
    centerB = comCurrRef.Get(molIndexB[0]);
    //Set the V1 to the the backbone of the large molecule
    if(molRef.NumAtoms(kindL) == 1) {
      SetBasis(cavB, prng.RandomUnitVect());
    } else {
      uint start = molRef.MolStart(molIndexB[0]) + largeBB[0];
      uint end = molRef.MolStart(molIndexB[0]) + largeBB[1];
      SetBasis(cavB, boxDimRef.MinImage(coordCurrRef.Difference(start, end), sourceBox));
    }
    //Calculate inverse matrix for cav. Here Inv = Transpose
    TransposeMatrix(invCavB, cavB);
    //find how many of KindS exist in this centerB (COM of kindL)
    calcEnRef.FindMolInCavity(molInCav, centerB, cavity, invCavB,
                              sourceBox, kindS, exchangeRatio);
    numSCavB = molInCav[kindS].size();
  }
  return state;
}


inline uint IntraMoleculeExchange1::GetBoxPairAndMol(const real subDraw,
    const real movPerc)
{
  uint state = mv::fail_state::NO_FAIL;
  overlap = false;

#if ENSEMBLE == GCMC
  sourceBox = mv::BOX0;
#else
  prng.PickBox(sourceBox, subDraw, movPerc);
#endif

  SetExchangeData();

  molIndexA.clear();
  kindIndexA.clear();
  molIndexB.clear();
  kindIndexB.clear();

  newMolA.clear();
  oldMolA.clear();
  newMolB.clear();
  oldMolB.clear();

  //Pick the exchange number of kindS in cavity and a molecule of kindL
  //kindA = kindS, kindB = kindL
  state = PickMolInCav();

  if(state == mv::fail_state::NO_FAIL) {
    pStartA.clear();
    pStartB.clear();
    pStartA.resize(numInCavA);
    pStartB.resize(numInCavB);
    pLenA.clear();
    pLenB.clear();
    pLenA.resize(numInCavA);
    pLenB.resize(numInCavB);

    for(uint n = 0; n < numInCavA; n++) {
      pStartA[n] = pLenA[n] = 0;
      molRef.GetRangeStartLength(pStartA[n], pLenA[n], molIndexA[n]);
    }

    for(uint n = 0; n < numInCavB; n++) {
      pStartB[n] = pLenB[n] = 0;
      molRef.GetRangeStartLength(pStartB[n], pLenB[n], molIndexB[n]);
    }
  }

  return state;
}


inline uint IntraMoleculeExchange1::Prep(const real subDraw,
    const real movPerc)
{
  //AdjustExRatio();
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if(state == mv::fail_state::NO_FAIL) {
    //transfering type A from source
    for(uint n = 0; n < numInCavA; n++) {
      newMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
                                       sourceBox));
      oldMolA.push_back(cbmc::TrialMol(molRef.kinds[kindIndexA[n]], boxDimRef,
                                       sourceBox));
    }

    for(uint n = 0; n < numInCavB; n++) {
      //transfering type B from source
      newMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
                                       sourceBox));
      oldMolB.push_back(cbmc::TrialMol(molRef.kinds[kindIndexB[n]], boxDimRef,
                                       sourceBox));
    }

    //set the old coordinate after unwrap them
    for(uint n = 0; n < numInCavA; n++) {
      XYZArray molA(pLenA[n]);
      coordCurrRef.CopyRange(molA, pStartA[n], 0, pLenA[n]);
      boxDimRef.UnwrapPBC(molA, sourceBox, comCurrRef.Get(molIndexA[n]));
      oldMolA[n].SetCoords(molA, 0);
      //copy cavA matrix to slant the old trial of molA
      oldMolA[n].SetCavMatrix(cavA);
      //set coordinate of moleA to newMolA, later it will shift to centerB
      newMolA[n].SetCoords(molA, 0);
      //copy cavB matrix to slant the new trial of molA
      newMolA[n].SetCavMatrix(cavB);
    }

    for(uint n = 0; n < numInCavB; n++) {
      XYZArray molB(pLenB[n]);
      coordCurrRef.CopyRange(molB, pStartB[n], 0, pLenB[n]);
      boxDimRef.UnwrapPBC(molB, sourceBox, comCurrRef.Get(molIndexB[n]));
      oldMolB[n].SetCoords(molB, 0);
      //copy cavB matrix to slant the old trial of molB
      oldMolB[n].SetCavMatrix(cavB);
      //set coordinate of moleB to newMolB, later it will shift to centerA
      newMolB[n].SetCoords(molB, 0);
      //copy cavA matrix to slant the new trial of molB
      newMolB[n].SetCavMatrix(cavA);
    }

    //SetSeed(has cavity, COM is fixed, rotate around Backbone)
    for(uint n = 0; n < numInCavB; n++) {
      //Inserting molB from centerB to the centerA
      newMolB[n].SetSeed(centerA, cavity, true, true, true);
      // Set the Backbone of large molecule to be inserted
      newMolB[n].SetBackBone(largeBB);
      //perform rotational trial move for oldMolB
      oldMolB[n].SetSeed(centerB, cavity, true, true, true);
      // Set the Backbone of large molecule to be deleted
      oldMolB[n].SetBackBone(largeBB);
    }

    for(uint n = 0; n < numInCavA; n++) {
      //Inserting molA from cavity(centerA) to the cavityB(centerB)
      newMolA[n].SetSeed(centerB, cavity, true, false, false);
      //perform trial move in cavity in sourceBox for oldMolA
      oldMolA[n].SetSeed(centerA, cavity, true, false, false);
    }
  }

  return state;
}


inline uint IntraMoleculeExchange1::Transform()
{
  //Calc old energy before deleting
  for(uint n = 0; n < numInCavA; n++) {
    cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexA[n]].BuildIDOld(oldMolA[n], molIndexA[n]);
    //Add bonded energy because we dont considered in DCRotate.cpp
    oldMolA[n].AddEnergy(calcEnRef.MoleculeIntra(oldMolA[n], molIndexA[n]));
  }

  //Calc old energy before deleting
  for(uint n = 0; n < numInCavB; n++) {
    cellList.RemoveMol(molIndexB[n], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexB[n]].BuildIDOld(oldMolB[n], molIndexB[n]);
    //Add bonded energy because we dont considered in DCRotate.cpp
    oldMolB[n].AddEnergy(calcEnRef.MoleculeIntra(oldMolB[n], molIndexB[n]));
  }

  //Insert kindL to cavity of  center A
  for(uint n = 0; n < numInCavB; n++) {
    molRef.kinds[kindIndexB[n]].BuildIDNew(newMolB[n], molIndexB[n]);
    ShiftMol(n, false);
    cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
    //Add bonded energy because we dont considered in DCRotate.cpp
    newMolB[n].AddEnergy(calcEnRef.MoleculeIntra(newMolB[n], molIndexB[n]));
    overlap |= newMolB[n].HasOverlap();
  }

  //Insert kindS to cavity of center B
  for(uint n = 0; n < numInCavA; n++) {
    molRef.kinds[kindIndexA[n]].BuildIDNew(newMolA[n], molIndexA[n]);
    ShiftMol(n, true);
    cellList.AddMol(molIndexA[n], sourceBox, coordCurrRef);
    //Add bonded energy because we dont considered in DCRotate.cpp
    newMolA[n].AddEnergy(calcEnRef.MoleculeIntra(newMolA[n], molIndexA[n]));
    overlap |= newMolA[n].HasOverlap();
  }

  return mv::fail_state::NO_FAIL;
}


inline void IntraMoleculeExchange1::CalcEn()
{
  W_recip = 1.0;
  recipDiffA = 0.0, recipDiffB = 0.0;
  correctDiff = 0.0;

  if(!overlap) {
    recipDiffA = calcEwald->SwapRecip(newMolA, oldMolA);
    recipDiffB = calcEwald->SwapRecip(newMolB, oldMolB);

    //No need to contribute the self and correction energy since insertion
    //and deletion are rigid body
    W_recip = exp(-1.0 * ffRef.beta * (recipDiffA + recipDiffB));
  }
}

inline real IntraMoleculeExchange1::GetCoeff() const
{
  real ratioF =  num::Factorial(numSCavA) * num::Factorial(numSCavB) /
                   (num::Factorial(numSCavA - exchangeRatio) *
                    num::Factorial(numSCavB + exchangeRatio));

  return ratioF;
}

inline void IntraMoleculeExchange1::ShiftMol(const uint n, const bool typeA)
{
  if(typeA) {
    //update coordinate of molecule typeA
    newMolA[n].GetCoords().CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], sourceBox);
  } else {
    //update coordinate of molecule typeA
    newMolB[n].GetCoords().CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], sourceBox);
  }
}

inline void IntraMoleculeExchange1::RecoverMol(const uint n, const bool typeA)
{
  if(typeA) {
    XYZArray molA(pLenA[n]);
    oldMolA[n].GetCoords().CopyRange(molA, 0, 0, pLenA[n]);
    boxDimRef.WrapPBC(molA, sourceBox);

    molA.CopyRange(coordCurrRef, 0, pStartA[n], pLenA[n]);
    comCurrRef.SetNew(molIndexA[n], sourceBox);
  } else {
    XYZArray molB(pLenB[n]);
    oldMolB[n].GetCoords().CopyRange(molB, 0, 0, pLenB[n]);
    boxDimRef.WrapPBC(molB, sourceBox);

    molB.CopyRange(coordCurrRef, 0, pStartB[n], pLenB[n]);
    comCurrRef.SetNew(molIndexB[n], sourceBox);
  }
}

inline void IntraMoleculeExchange1::Accept(const uint rejectState,
    const uint step)
{
  bool result;
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    real molTransCoeff = GetCoeff();
    real Wrat = W_recip;

    for(uint n = 0; n < numInCavA; n++) {
      Wrat *= newMolA[n].GetWeight() / oldMolA[n].GetWeight();
    }

    for(uint n = 0; n < numInCavB; n++) {
      Wrat *= newMolB[n].GetWeight() / oldMolB[n].GetWeight();
    }

    if(!overlap) {
      result = prng() < molTransCoeff * Wrat;
    } else {
      result = false;
    }

    if(result) {
      //Add rest of energy.
      for(uint n = 0; n < numInCavB; n++) {
        sysPotRef.boxEnergy[sourceBox] += newMolB[n].GetEnergy();
        sysPotRef.boxEnergy[sourceBox] -= oldMolB[n].GetEnergy();
      }

      for(uint n = 0; n < numInCavA; n++) {
        sysPotRef.boxEnergy[sourceBox] -= oldMolA[n].GetEnergy();
        sysPotRef.boxEnergy[sourceBox] += newMolA[n].GetEnergy();
      }

      //Add Reciprocal energy
      sysPotRef.boxEnergy[sourceBox].recip += recipDiffA + recipDiffB;

      //Add correction energy
      sysPotRef.boxEnergy[sourceBox].correction += correctDiff;

      //Update reciprocal
      calcEwald->UpdateRecip(sourceBox);

      //molA and molB already added to cellist

      //Retotal
      sysPotRef.Total();
    } else {
      //transfer molA
      for(uint n = 0; n < numInCavA; n++) {
        cellList.RemoveMol(molIndexA[n], sourceBox, coordCurrRef);
        RecoverMol(n, true);
        cellList.AddMol(molIndexA[n], sourceBox, coordCurrRef);
      }
      //transfer molB
      for(uint n = 0; n < numInCavB; n++) {
        cellList.RemoveMol(molIndexB[n], sourceBox, coordCurrRef);
        RecoverMol(n, false);
        cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
      }
    }
  } else {
    //else we didn't even try because we knew it would fail
    result = false;
  }

  moveSetRef.Update(mv::INTRA_MEMC, result, step, sourceBox);
  //If we consider total aceeptance of S->L and L->S
  AcceptKind(result, kindS * molRef.GetKindsCount() + kindL, sourceBox);
  AcceptKind(result, kindL * molRef.GetKindsCount() + kindS, sourceBox);

}

#endif
