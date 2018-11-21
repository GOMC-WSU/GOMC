/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CFCMC_H
#define CFCMC_H

#if ENSEMBLE==GCMC || ENSEMBLE==GEMC

#include "MoveBase.h"
#include "TrialMol.h"
#include <numeric>

using std::vector;

class CFCMC : public MoveBase
{
public:

  CFCMC(System &sys, StaticVals const& statV) :
    ffRef(statV.forcefield), molLookRef(sys.molLookupRef),
    lambdaRef(sys.lambdaRef), MoveBase(sys, statV) 
    {
      totalMolecule = comCurrRef.Count();
      eqCycle = 1;
      lambdaWindow = 20;
      histFreq = 1000;
      lambdaMax = 1.0 / (double)(lambdaWindow);
      nuTolerance = 1e-6;
      uint totKind = molRef.GetKindsCount();
      nu.resize(totKind, 0.01);
      hist.resize(totKind);
      bias.resize(totKind);
      for(uint k = 0; k < totKind; k++) {
        hist[k].resize(lambdaWindow, 0);
        bias[k].resize(lambdaWindow, 0.0);
      }
    }

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual uint Transform();
  virtual void CalcEn();
  virtual void Accept(const uint earlyReject, const uint step);
  virtual void PrintAcceptKind();

private:

  double GetCoeff() const;
  uint GetBoxPairAndMol(const double subDraw, const double movPerc);
  void ShiftMolToSourceBox();
  void ShiftMolToDestBox();
  void UpdateBias();

  MolPick molPick;
  uint totalMolecule;
  uint sourceBox, destBox;
  uint pStartCFCMC, pLenCFCMC;
  uint molIndex, kindIndex;
  uint lambdaWindow, histFreq;
  uint eqCycle;
  vector< vector<uint> > hist;

  double W_tc, W_recip;
  double correct_old, correct_new, self_old, self_new;
  double lambdaMax, nuTolerance;
  double lambdaOld, lambdaNew;
  double *lambdaRef;
  vector< vector<double> > bias;
  vector< double > nu;


  cbmc::TrialMol oldMol, newMol;
  Intermolecular tcLose, tcGain, recipLose, recipGain;
  MoleculeLookup & molLookRef;
  Forcefield const& ffRef;
};

void CFCMC::PrintAcceptKind() {
  for(uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted CFCMC-Transfer ", molRef.kinds[k].name.c_str());
    for(uint b = 0; b < BOX_TOTAL; b++) {
      if(moveSetRef.GetTrial(b, mv::CFCMC, k) > 0)
        printf("%10.5f ", (100.0 * moveSetRef.GetAccept(b, mv::CFCMC, k)));
      else
        printf("%10.5f ", 0.0);
    }
    std::cout << std::endl;
  }
}

inline uint CFCMC::GetBoxPairAndMol(const double subDraw, const double movPerc)
{
  // Need to call a function to pick a molecule that is not fixed but cannot be
  // swap between boxes. (beta != 1, beta !=2)
  uint state = prng.PickMolAndBoxPair2(molIndex, kindIndex, sourceBox, destBox,
                                       subDraw, movPerc);
#if ENSEMBLE == GCMC
  if(state == mv::fail_state::NO_MOL_OF_KIND_IN_BOX && sourceBox == mv::BOX1) {
    std::cout << "Error: There are no molecules of kind " <<
              molRef.kinds[kindIndex].name << " left in reservoir.\n";
    exit(EXIT_FAILURE);
  }
#endif

  if (state == mv::fail_state::NO_FAIL) {
    pStartCFCMC = pLenCFCMC = 0;
    molRef.GetRangeStartLength(pStartCFCMC, pLenCFCMC, molIndex);
  }
  return state;
}

inline uint CFCMC::Prep(const double subDraw, const double movPerc)
{
  overlap = false;
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if (state == mv::fail_state::NO_FAIL) {
    lambdaOld = 1.0;
    lambdaNew = 1.0 - lambdaMax;
    //Start with full interaction in sourceBox, zero interaction in destBox
    lambdaRef[sourceBox * totalMolecule + totalMolecule] = lambdaOld;
    lambdaRef[destBox * totalMolecule + totalMolecule] = 1.0 - lambdaOld;
    newMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMol = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMol.SetCoords(coordCurrRef, pStartCFCMC);
  }
  return state;
}


inline uint CFCMC::Transform()
{
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  molRef.kinds[kindIndex].Build(oldMol, newMol, molIndex);
  //Transfer the molecule to destBox
  newMol.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
  comCurrRef.SetNew(molIndex, destBox);
  molLookRef.ShiftMolBox(molIndex, sourceBox, destBox, kindIndex);
  cellList.AddMol(molIndex, destBox, coordCurrRef);
  //overlap = newMol.HasOverlap();
  return mv::fail_state::NO_FAIL;
}

inline void CFCMC::CalcEn()
{
  W_tc = 1.0;
  W_recip = 1.0;
  correct_old = 0.0;
  correct_new = 0.0;
  self_old = 0.0;
  self_new = 0.0;

  if(ffRef.useLRC) {
    tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
    tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
    W_tc = exp(-1.0 * ffRef.beta * (tcGain.energy + tcLose.energy));
  }

  if (newMol.GetWeight() != 0.0 && !overlap) {
    correct_new = calcEwald->SwapCorrection(newMol);
    correct_old = calcEwald->SwapCorrection(oldMol);
    self_new = calcEwald->SwapSelf(newMol);
    self_old = calcEwald->SwapSelf(oldMol);
    recipGain.energy =
      calcEwald->SwapDestRecip(newMol, destBox, sourceBox, molIndex);
    recipLose.energy =
      calcEwald->SwapSourceRecip(oldMol, sourceBox, molIndex);
    //need to contribute the self and correction energy
    W_recip = exp(-1.0 * ffRef.beta * (recipGain.energy + recipLose.energy +
                                       correct_new - correct_old +
                                       self_new - self_old));
  }

}

inline double CFCMC::GetCoeff() const
{
#if ENSEMBLE == GEMC
  return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) /
         (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
         boxDimRef.volume[destBox] * boxDimRef.volInv[sourceBox];
#elif ENSEMBLE == GCMC
  if (sourceBox == mv::BOX0) { //Delete case
    if(ffRef.isFugacity) {
      return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) *
             boxDimRef.volInv[sourceBox] /
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) *
             boxDimRef.volInv[sourceBox] *
             exp(-BETA * molRef.kinds[kindIndex].chemPot);
    }
  } else { //Insertion case
    if(ffRef.isFugacity) {
      return boxDimRef.volume[destBox] /
             (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
             (BETA * molRef.kinds[kindIndex].chemPot);
    } else {
      return boxDimRef.volume[destBox] /
             (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
             exp(BETA * molRef.kinds[kindIndex].chemPot);
    }
  }
#endif
}

inline void CFCMC::Accept(const uint rejectState, const uint step)
{
  bool result;
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    double molTransCoeff = GetCoeff();
    double Wo = oldMol.GetWeight();
    double Wn = newMol.GetWeight();
    double Wrat = Wn / Wo * W_tc * W_recip;

    //safety to make sure move will be rejected in overlap case
    if(!overlap) {
      result = prng() < molTransCoeff * Wrat;
    } else
      result = false;

    if(result) {
      //Add tail corrections
      sysPotRef.boxEnergy[sourceBox].tc += tcLose.energy;
      sysPotRef.boxEnergy[destBox].tc += tcGain.energy;
      //Add rest of energy.
      sysPotRef.boxEnergy[sourceBox] -= oldMol.GetEnergy();
      sysPotRef.boxEnergy[destBox] += newMol.GetEnergy();
      //Add Reciprocal energy
      sysPotRef.boxEnergy[sourceBox].recip += recipLose.energy;
      sysPotRef.boxEnergy[destBox].recip += recipGain.energy;
      //Add correction energy
      sysPotRef.boxEnergy[sourceBox].correction -= correct_old;
      sysPotRef.boxEnergy[destBox].correction += correct_new;
      //Add self energy
      sysPotRef.boxEnergy[sourceBox].self -= self_old;
      sysPotRef.boxEnergy[destBox].self += self_new;

      //Calculate the change of force 
      calcEnRef.MoleculeForceSub(atomForceRef, molForceRef, molIndex,sourceBox);
      calcEnRef.MoleculeForceAdd(newMol.GetCoords(), atomForceRef, molForceRef,
                                 molIndex, destBox);

      //Set coordinates, new COM; shift index to new box's list
      newMol.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
      comCurrRef.SetNew(molIndex, destBox);
      molLookRef.ShiftMolBox(molIndex, sourceBox, destBox,
                             kindIndex);
      cellList.AddMol(molIndex, destBox, coordCurrRef);


      //Zero out box energies to prevent small number
      //errors in double.
      if (molLookRef.NumInBox(sourceBox) == 0) {
        sysPotRef.boxEnergy[sourceBox].Zero();
        sysPotRef.boxVirial[sourceBox].Zero();
      } else if (molLookRef.NumInBox(sourceBox) == 1) {
        sysPotRef.boxEnergy[sourceBox].inter = 0;
        sysPotRef.boxVirial[sourceBox].inter = 0;
        sysPotRef.boxEnergy[sourceBox].real = 0;
        sysPotRef.boxVirial[sourceBox].real = 0;
      }

      for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
        calcEwald->UpdateRecip(b);
      }

      //Retotal
      sysPotRef.Total();
    } else {
      cellList.AddMol(molIndex, sourceBox, coordCurrRef);
      //when weight is 0, MolDestSwap() will not be executed, thus cos/sin
      //molRef will not be changed. Also since no memcpy, doing restore
      //results in memory overwrite
      if (newMol.GetWeight() != 0.0 && !overlap) {
        calcEwald->RestoreMol(molIndex);
      }
    }

  } else //we didn't even try because we knew it would fail
    result = false;

  moveSetRef.Update(mv::CFCMC, result, step, destBox, kindIndex);
}


void CFCMC::ShiftMolToSourceBox()
{
  cellList.RemoveMol(molIndex, destBox, coordCurrRef);
  //Set coordinates, new COM; shift index to new box's list
  oldMol.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
  comCurrRef.SetNew(molIndex, sourceBox);
  molLookRef.ShiftMolBox(molIndex, destBox, sourceBox, kindIndex);
  cellList.AddMol(molIndex, sourceBox, coordCurrRef);
}

void CFCMC::ShiftMolToDestBox()
{ 
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  //Set coordinates, new COM; shift index to new box's list
  newMol.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
  comCurrRef.SetNew(molIndex, destBox);
  molLookRef.ShiftMolBox(molIndex, sourceBox, destBox, kindIndex);
  cellList.AddMol(molIndex, destBox, coordCurrRef);
}

void CFCMC::UpdateBias()
{
  if(nu[kindIndex] <= nuTolerance)
    return;

  uint idx = 0;
  //Find the index that lead to insert or delete in box0
  if(sourceBox == mv::BOX0) {
    idx = (uint)(lambdaOld * lambdaWindow);
  } else {
    idx = (uint)((1.0 - lambdaOld) * lambdaWindow);
  }

  hist[kindIndex][idx] += 1;
  bias[kindIndex][idx] -= nu[kindIndex];

  uint trial = std::accumulate(hist[kindIndex].begin(),
                               hist[kindIndex].end(), 0);

  if((trial + 1) % histFreq == 0) {
    uint maxVisited = *max_element(hist[kindIndex].begin(),
                                   hist[kindIndex].end());
    uint minVisited = *min_element(hist[kindIndex].begin(),
                                   hist[kindIndex].end());
    //check to see if all the bin is visited atleast 30% of 
    // the most visited bin.                               
    if(minVisited > 0.3 * maxVisited) {
      nu[kindIndex] = 0.5 * nu[kindIndex];
      std::fill_n(hist[kindIndex].begin(), lambdaWindow, 0);
    }                               
  }
}

#endif

#endif
