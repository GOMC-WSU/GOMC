#ifndef INTRAMOLECULEEXCHANGE2_H
#define INTRAMOLECULEEXCHANGE2_H

#include "TrialMol.h"
#include "GeomLib.h"
#include "IntraMoleculeExchange1.h"
#include <cmath>

using std::vector;
using namespace geom;

// Intra Molecule Exchange Move:
// KindA is small kind. KindB is large kind
// Orientation center of cavA is on COM of kindS, aligned with kindS backbone.
// Orientation center of cavB is on COM of kindL, aligned with kindL backbon.
// Delete the exchangeRatio kindS from cavA, and 1 kindL from cavB.
// Insert the exchangeRatio kindS to cavB and 1 kindL inside the cavA.

class IntraMoleculeExchange2 : public IntraMoleculeExchange1
{
public:

  IntraMoleculeExchange2(System &sys, StaticVals const& statV) :
    IntraMoleculeExchange1(sys, statV)
  {
    if(enableID) {
      SetMEMC(statV);
    }
  }

  virtual uint Prep(const real subDraw, const real movPerc);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const uint step);

protected:

  virtual void AdjustExRatio();
  virtual void SetMEMC(StaticVals const& statV);
  virtual void SetExchangeData();
  virtual uint PickMolInCav();
  virtual real GetCoeff() const;

  uint smallBB[2];
  //To store total sets of exchange pairs
  vector< vector<uint> > smallBBVec;
};

inline void IntraMoleculeExchange2::SetMEMC(StaticVals const& statV)
{
  for(uint t = 0; t < exchangeRatioVec.size(); t++) {
    smallBB[0] = smallBB[1] = -1;
    for(uint i = 0; i < molRef.kinds[kindSVec[t]].NumAtoms(); i++) {
      if(molRef.kinds[kindSVec[t]].atomNames[i] == statV.intraMemcVal.smallBBAtom1[t]) {
        smallBB[0] = i;
      }
      if(molRef.kinds[kindSVec[t]].atomNames[i] == statV.intraMemcVal.smallBBAtom2[t]) {
        smallBB[1] = i;
      }
    }

    for(uint i = 0; i < 2; i++) {
      if(smallBB[i] == -1) {
        printf("Error: Atom name %s or %s was not found in %s residue.\n",
               statV.intraMemcVal.smallBBAtom1[t].c_str(),
               statV.intraMemcVal.smallBBAtom2[t].c_str(),
               statV.intraMemcVal.smallKind[t].c_str());
        exit(EXIT_FAILURE);
      }
    }

    if(molRef.kinds[kindSVec[t]].NumAtoms() > 1) {
      if(smallBB[0] == smallBB[1]) {
        printf("Error: Atom names in small molecule backbone cannot be same!\n");
        exit(EXIT_FAILURE);
      }
    }
    vector<uint> temp(smallBB, smallBB + 2);
    smallBBVec.push_back(temp);
  }
}

inline void IntraMoleculeExchange2::AdjustExRatio()
{
  if(((counter + 1) % perAdjust) == 0) {
    uint exMax = ceil((float)molInCavCount / (float)perAdjust);
    uint exMin = 1;

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
    printf("ExchangeRatio: %d, Average kindS in cavity: %d \n", exchangeRatio,
           exMax);
  }
}

inline void IntraMoleculeExchange2::SetExchangeData()
{
  uint exType = prng.randIntExc(exchangeRatioVec.size());
  kindS = kindSVec[exType];
  kindL = kindLVec[exType];
  exchangeRatio = exchangeRatioVec[exType];
  largeBB[0] = largeBBVec[exType][0];
  largeBB[1] = largeBBVec[exType][1];
  smallBB[0] = smallBBVec[exType][0];
  smallBB[1] = smallBBVec[exType][1];
}

inline uint IntraMoleculeExchange2::PickMolInCav()
{
  uint state = mv::fail_state::NO_FAIL;
  //pick a random small kind in dense phase and use the COM as cavity center
  uint pickedS, pickedKS;
  state = prng.PickMol(kindS, pickedKS, pickedS, sourceBox);
  if(state == mv::fail_state::NO_FAIL) {
    centerA = comCurrRef.Get(pickedS);
    //If we want to orient the cavity with backbone of picked small mol
    if(molRef.NumAtoms(kindS) == 1) {
      SetBasis(cavA, prng.RandomUnitVect());
    } else {
      uint start = molRef.MolStart(pickedS) + smallBB[0];
      uint end = molRef.MolStart(pickedS) + smallBB[1];
      SetBasis(cavA, boxDimRef.MinImage(coordCurrRef.Difference(start, end), sourceBox));
    }
    //Calculate inverse matrix for cav here Inv = transpose
    TransposeMatrix(invCavA, cavA);

    //Find the small molecule kind in the cavityA
    if((exchangeRatio == 1) || calcEnRef.FindMolInCavity(molInCav, centerA,
        cavity, invCavA, sourceBox, kindS, exchangeRatio)) {
      //Find the exchangeRatio number of molecules kindS in cavity
      numInCavA = exchangeRatio;
      //add the random picked small molecule to the list.
      molIndexA.push_back(pickedS);
      kindIndexA.push_back(pickedKS);
      if(exchangeRatio == 1) {
        numSCavA = 1;
      } else {
        numSCavA = molInCav[kindS].size();
        //delete the picked small molecule from list
        for(uint s = 0; s < numSCavA; s++) {
          if(pickedS == molInCav[kindS][s])
            molInCav[kindS].erase(molInCav[kindS].begin() + s);
        }
      }

      for(uint n = 1; n < numInCavA; n++) {
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
  }

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
    if(exchangeRatio == 1) {
      numSCavB = 0;
    } else {
      //find how many of KindS exist in this centerB (COM of kindL)
      calcEnRef.FindMolInCavity(molInCav, centerB, cavity, invCavB,
                                sourceBox, kindS, exchangeRatio);
      numSCavB = molInCav[kindS].size();
    }
  }
  return state;
}


inline uint IntraMoleculeExchange2::Prep(const real subDraw,
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
      if(n == 0) {
        //Inserting molA from cavity(centerA) to the cavityB(centerB)
        //COM is fixed, rotation around backboe
        newMolA[n].SetSeed(centerB, cavity, true, true, true);
        // Set the Backbone of small molecule to be inserted
        newMolA[n].SetBackBone(smallBB);
        //perform trial move in cavity in sourceBox for oldMolA
        //COM is fixed, rotation around backboe
        oldMolA[n].SetSeed(centerA, cavity, true, true, true);
        // Set the Backbone of small molecule to be deleted
        oldMolA[n].SetBackBone(smallBB);
      } else {
        //Inserting molA from cavity(centerA) to the cavityB(centerB)
        newMolA[n].SetSeed(centerB, cavity, true, false, false);
        //perform trial move in cavity in sourceBox for oldMolA
        oldMolA[n].SetSeed(centerA, cavity, true, false, false);
      }
    }
  }

  return state;
}


inline uint IntraMoleculeExchange2::Transform()
{
  ///Remove the fixed COM kindS at the end because we insert it at first
  for(uint n = numInCavA; n > 0; n--) {
    cellList.RemoveMol(molIndexA[n - 1], sourceBox, coordCurrRef);
    molRef.kinds[kindIndexA[n - 1]].BuildIDOld(oldMolA[n - 1], molIndexA[n - 1]);
    //Add bonded energy because we dont considered in DCRotate.cpp
    oldMolA[n - 1].AddEnergy(calcEnRef.MoleculeIntra(oldMolA[n - 1], molIndexA[n - 1]));
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

inline void IntraMoleculeExchange2::CalcEn()
{
  IntraMoleculeExchange1::CalcEn();
}

inline real IntraMoleculeExchange2::GetCoeff() const
{
  real ratioF =  num::Factorial(numSCavA - 1) * num::Factorial(numSCavB) /
                   (num::Factorial(numSCavA - exchangeRatio) *
                    num::Factorial(numSCavB + exchangeRatio - 1));

  return ratioF;
}

inline void IntraMoleculeExchange2::Accept(const uint rejectState,
    const uint step)
{
  IntraMoleculeExchange1::Accept(rejectState, step);
}

#endif
