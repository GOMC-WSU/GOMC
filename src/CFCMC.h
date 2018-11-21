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
      relaxMolecule.resize(BOX_TOTAL);
      totalMolecule = comCurrRef.Count();
      eqCycle = 1;
      relaxRadiusSq = 10 * 10;
      lambdaWindow = 10;
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
  void FindRelaxingMolecules(uint box);
  void InflatingMolecule();
  bool AcceptInflating();
  void CalcEnCFCMC(bool calcNewEn);
  void RelaxingMolecules();

  MolPick molPick;
  uint totalMolecule;
  uint sourceBox, destBox;
  uint pStartCFCMC, pLenCFCMC;
  uint molIndex, kindIndex;
  uint lambdaWindow, histFreq;
  uint eqCycle;
  vector< vector<uint> > hist;
  vector< vector<uint> > relaxMolecule;

  double W_tc, W_recip;
  double correct_old, correct_new, self_old, self_new;
  double lambdaMax, nuTolerance;
  double lambdaOld, lambdaNew;
  double relaxRadiusSq;
  double *lambdaRef;
  vector< vector<double> > bias;
  vector< double > nu;

  bool accepted;


  cbmc::TrialMol oldMol, newMol;
  Intermolecular tcLose, tcGain, recipLose, recipGain;
  Energy oldEnergy[BOX_TOTAL], newEnergy[BOX_TOTAL];
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
    //Transform the picked molecule to random location in dest box
    XYZ oldCOM = comCurrRef.Get(molIndex);
    XYZ pickCOM;
    prng.FillWithRandom(pickCOM, boxDimRef, destBox);
    XYZArray temp(pLenCFCMC);
    coordCurrRef.CopyRange(temp, pStartCFCMC, 0, pLenCFCMC);
    boxDimRef.UnwrapPBC(temp, sourceBox, oldCOM);
    pickCOM -= oldCOM;
    temp.AddRange(0, pLenCFCMC, pickCOM);
    boxDimRef.WrapPBC(temp, destBox);
    newMol.SetCoords(temp, 0);
  }
  return state;
}


inline uint CFCMC::Transform()
{
  //Find the relaxing molecules in source box before shifting the molecule
  FindRelaxingMolecules(sourceBox);
  //Transfer the molecule to destBox
  ShiftMolToDestBox();
  //Find the relaxing molecules in dest box after shifting the molecule
  FindRelaxingMolecules(destBox);
  InflatingMolecule();
  return mv::fail_state::NO_FAIL;
}

inline void CFCMC::CalcEn() {
  return;
}

inline void CFCMC::CalcEnCFCMC(bool calcNewEn)
{
  W_tc = 1.0;
  W_recip = 1.0;
  correct_old = 0.0;
  correct_new = 0.0;
  self_old = 0.0;
  self_new = 0.0;
  tcLose.Zero();
  tcGain.Zero();


  uint idxNew = (uint)(lambdaNew * lambdaWindow);
  if(idxNew == 0 || idxNew == 1) {
    if(ffRef.useLRC) {
      tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
      tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
      W_tc = exp(-1.0 * ffRef.beta * (tcGain.energy + tcLose.energy));
    }
  }

  if(!calcNewEn) {
    //calculate inter and intra energy before changing lambda
    ShiftMolToSourceBox();
    calcEnRef.SingleMoleculeInter(oldEnergy[sourceBox], atomForceNew,
                                  molForceNew, molIndex, sourceBox);
    ShiftMolToDestBox();
    calcEnRef.SingleMoleculeInter(oldEnergy[destBox], atomForceNew,
                                  molForceNew, molIndex, destBox);

  } else {
    //calculate inter and intra energy after changing lambda
    ShiftMolToSourceBox();
    calcEnRef.SingleMoleculeInter(newEnergy[sourceBox], atomForceNew,
                                  molForceNew, molIndex, sourceBox);
    ShiftMolToDestBox();
    calcEnRef.SingleMoleculeInter(newEnergy[destBox], atomForceNew,
                                  molForceNew, molIndex, destBox);
  }
}

inline double CFCMC::GetCoeff() const
{
  uint idxNew = (uint)(lambdaNew * lambdaWindow);
  uint idxOld = (uint)(lambdaOld * lambdaWindow);
  double biasCoef = exp(bias[kindIndex][idxNew] - bias[kindIndex][idxOld]);

  if(lambdaNew < 1.0 and lambdaNew > 0.0) {
    return biasCoef;
  } else {
    #if ENSEMBLE == GEMC
      return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) /
            (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
            boxDimRef.volume[destBox] * boxDimRef.volInv[sourceBox] * biasCoef;
    #elif ENSEMBLE == GCMC
      if (sourceBox == mv::BOX0) { //Delete case
        if(ffRef.isFugacity) {
          return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) *
                boxDimRef.volInv[sourceBox] /
                (BETA * molRef.kinds[kindIndex].chemPot) * biasCoef;
        } else {
          return (double)(molLookRef.NumKindInBox(kindIndex, sourceBox)) *
                boxDimRef.volInv[sourceBox] *
                exp(-BETA * molRef.kinds[kindIndex].chemPot) * biasCoef;
        }
      } else { //Insertion case
        if(ffRef.isFugacity) {
          return boxDimRef.volume[destBox] /
                (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
                (BETA * molRef.kinds[kindIndex].chemPot) * biasCoef;
        } else {
          return boxDimRef.volume[destBox] /
                (double)(molLookRef.NumKindInBox(kindIndex, destBox) + 1) *
                exp(BETA * molRef.kinds[kindIndex].chemPot) * biasCoef;
        }
      }
    #endif
  }
}

inline void CFCMC::Accept(const uint rejectState, const uint step)
{
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    moveSetRef.Update(mv::CFCMC, accepted, step, destBox, kindIndex);
    if(!accepted) {   
      //Start with full interaction in sourceBox, zero interaction in destBox
      lambdaRef[sourceBox * totalMolecule + totalMolecule] = 1.0;
      lambdaRef[destBox * totalMolecule + totalMolecule] = 0.0;
    }
  } else{
    moveSetRef.Update(mv::CFCMC, false, step, destBox, kindIndex);
  }
}


inline void CFCMC::ShiftMolToSourceBox()
{
  cellList.RemoveMol(molIndex, destBox, coordCurrRef);
  //Set coordinates, new COM; shift index to new box's list
  oldMol.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
  comCurrRef.SetNew(molIndex, sourceBox);
  molLookRef.ShiftMolBox(molIndex, destBox, sourceBox, kindIndex);
  cellList.AddMol(molIndex, sourceBox, coordCurrRef);
}

inline void CFCMC::ShiftMolToDestBox()
{ 
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  //Set coordinates, new COM; shift index to new box's list
  newMol.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
  comCurrRef.SetNew(molIndex, destBox);
  molLookRef.ShiftMolBox(molIndex, sourceBox, destBox, kindIndex);
  cellList.AddMol(molIndex, destBox, coordCurrRef);
}

inline void CFCMC::UpdateBias()
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

inline void CFCMC::FindRelaxingMolecules(uint box)
{
  //Molecule has to be transfered before calling this function
  relaxMolecule[box].clear();
  XYZ center = comCurrRef.Get(molIndex);
  double minLengthSq;
  XYZ diff;

  MoleculeLookup::box_iterator n = molLookRef.BoxBegin(box);
  MoleculeLookup::box_iterator end = molLookRef.BoxEnd(box);
  while (n != end) {
    if(*n != molIndex) {
      diff = comCurrRef.Get(*n) - center;
      minLengthSq = boxDimRef.MinImage(diff, box).LengthSq();
      if( minLengthSq < relaxRadiusSq){
        //if molecule can move and its not the growing molecule
        if(!molLookRef.IsFix(*n)) {
          relaxMolecule[box].push_back(*n);
        }
      }
    }
    n++;
  }
}

inline void CFCMC::InflatingMolecule()
{
  do {
    UpdateBias();
    //Calculate the old energy
    CalcEnCFCMC(false);
    //Update the interaction in sourceBox and destBox
    lambdaRef[sourceBox * totalMolecule + totalMolecule] = lambdaNew;
    lambdaRef[destBox * totalMolecule + totalMolecule] = 1.0 - lambdaNew;
    //Calculate the new energy
    CalcEnCFCMC(true);
    accepted = AcceptInflating();
    if(accepted) {
      lambdaOld = lambdaNew;
      RelaxingMolecules();
    }
    //pick new lambda
    lambdaNew = lambdaOld + (prng.randInt(1) ? lambdaMax : -lambdaMax);
  } while(lambdaOld > 0.0 && lambdaOld < 1.0);
}

inline bool CFCMC::AcceptInflating()
{
  double molTransCoeff = GetCoeff();
  double W1 = exp(-BETA * (newEnergy[sourceBox].Total() -
                          oldEnergy[sourceBox].Total()));
  double W2 = exp(-BETA * (newEnergy[destBox].Total() -
                          oldEnergy[destBox].Total()));
  double Wrat = W1 * W2 * W_tc * W_recip;

  bool result = prng() < molTransCoeff * Wrat;

  if(result) {
    //Add tail corrections
    sysPotRef.boxEnergy[sourceBox].tc += tcLose.energy;
    sysPotRef.boxEnergy[destBox].tc += tcGain.energy;
    //Add rest of energy.
    sysPotRef.boxEnergy[sourceBox] -= oldEnergy[sourceBox];
    sysPotRef.boxEnergy[sourceBox] += newEnergy[sourceBox];
    sysPotRef.boxEnergy[destBox] -= oldEnergy[destBox];
    sysPotRef.boxEnergy[destBox] += newEnergy[destBox];

    //Retotal
    sysPotRef.Total();
  }

  return result;
}


inline void CFCMC::RelaxingMolecules()
{
  return;
}


#endif

#endif
