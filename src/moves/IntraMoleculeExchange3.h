#ifndef INTRAMOLECULEEXCHANGE3_H
#define INTRAMOLECULEEXCHANGE3_H

#include "TrialMol.h"
#include "GeomLib.h"
#include "IntraMoleculeExchange1.h"
#include <cmath>

using std::vector;
using namespace geom;

// Intra Molecule Exchange Move:
// KindA is small kind. KindB is large kind
// Orientation center of cavA is on COM of kindS, random orientation.
// Orientation center of cavB is on COM of kindL, aligned with kindL backbon.
// Delete the exchangeRatio kindS from cavA, and 1 kindL from cavB.
// Insert the exchangeRatio kindS to cavB and 1 kindL inside the cavA.
//Use CD-CBMC to build kindL

class IntraMoleculeExchange3 : public IntraMoleculeExchange1
{
public:

  IntraMoleculeExchange3(System &sys, StaticVals const& statV) :
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
};

inline void IntraMoleculeExchange3::SetMEMC(StaticVals const& statV)
{
  for(uint t = 0; t < exchangeRatioVec.size(); t++) {
    if(largeBBVec[t][0] != largeBBVec[t][1]) {
      printf("Error: In ME-3 move, atom name of backbone should be same.\n");
      printf("Atom names in backbone was set to %s or %s in %s residue.\n",
             statV.intraMemcVal.largeBBAtom1[t].c_str(),
             statV.intraMemcVal.largeBBAtom2[t].c_str(),
             statV.intraMemcVal.largeKind[t].c_str());
      exit(EXIT_FAILURE);
    }
  }
}

inline void IntraMoleculeExchange3::AdjustExRatio()
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

inline void IntraMoleculeExchange3::SetExchangeData()
{
  IntraMoleculeExchange1::SetExchangeData();
}

inline uint IntraMoleculeExchange3::PickMolInCav()
{
  uint state = mv::fail_state::NO_FAIL;
  //pick a random small kind in dense phase and use the COM as cavity center
  uint pickedS, pickedKS;
  state = prng.PickMol(kindS, pickedKS, pickedS, sourceBox);
  if(state == mv::fail_state::NO_FAIL) {
    centerA = comCurrRef.Get(pickedS);
    //Pick random vector and find two vectors that are perpendicular to V1
    SetBasis(cavA, prng.RandomUnitVect());
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
    //use the predefine atom in kindL as the centerB
    uint start = molRef.MolStart(molIndexB[0]) + largeBB[0];
    centerB = coordCurrRef.Get(start);
    //Set V1 to a random vector and calculate two vector perpendicular to V1
    SetBasis(cavB, prng.RandomUnitVect());
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


inline uint IntraMoleculeExchange3::Prep(const real subDraw,
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
        //COM is fixed, rotation around sphere
        newMolA[n].SetSeed(centerB, cavity, true, true, false);
        //perform trial move in cavity in sourceBox for oldMolA
        //COM is fixed, rotation around sphere
        oldMolA[n].SetSeed(centerA, cavity, true, true, false);
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


inline uint IntraMoleculeExchange3::Transform()
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
    molRef.kinds[kindIndexB[n]].BuildGrowOld(oldMolB[n], molIndexB[n]);
  }

  //Insert kindL to cavity of  center A using CD-CBMC
  for(uint n = 0; n < numInCavB; n++) {
    molRef.kinds[kindIndexB[n]].BuildGrowNew(newMolB[n], molIndexB[n]);
    ShiftMol(n, false);
    cellList.AddMol(molIndexB[n], sourceBox, coordCurrRef);
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

inline void IntraMoleculeExchange3::CalcEn()
{
  W_recip = 1.0;
  recipDiffA = 0.0, recipDiffB = 0.0;
  correctDiff = 0.0;
  //No need to calculate the correction term for kindS since it is
  // inserted rigid body. We just need it for kindL
  if(!overlap) {
    for(uint n = 0; n < numInCavB; n++) {
      correctDiff += calcEwald->SwapCorrection(newMolB[n]);
      correctDiff -= calcEwald->SwapCorrection(oldMolB[n]);
    }
    recipDiffA = calcEwald->SwapRecip(newMolA, oldMolA);
    recipDiffB = calcEwald->SwapRecip(newMolB, oldMolB);

    W_recip = exp(-1.0 * ffRef.beta * (recipDiffA + recipDiffB +
                                       correctDiff));
  }
}

inline real IntraMoleculeExchange3::GetCoeff() const
{
  real ratioF =  num::Factorial(numSCavA - 1) * num::Factorial(numSCavB) /
                   (num::Factorial(numSCavA - exchangeRatio) *
                    num::Factorial(numSCavB + exchangeRatio - 1));

  return ratioF;
}

inline void IntraMoleculeExchange3::Accept(const uint rejectState,
    const uint step)
{
  IntraMoleculeExchange1::Accept(rejectState, step);
}

#endif
