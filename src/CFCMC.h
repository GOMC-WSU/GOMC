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
      hist.resize(BOX_TOTAL);
      bias.resize(BOX_TOTAL);
      for(uint b = 0; b < BOX_TOTAL; b++) {
	hist[b].resize(totKind);
	bias[b].resize(totKind);
      }
      for(uint b = 0; b < BOX_TOTAL; b++) {
	for(uint k = 0; k < totKind; k++) {
	  hist[b][k].resize(lambdaWindow + 1, 0);
	  bias[b][k].resize(lambdaWindow + 1, 0.0);
	}
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
  void AcceptRelaxing(uint box);
  void CalcEnCFCMC(bool calcNewEn);
  void CalcEnRelaxing(uint box);
  void TransformRelaxing(uint box);
  void RelaxingMolecules();
  uint GetLambdaIdx(double lambda) const 
    {return ((uint)(lambda * lambdaWindow));}
  

  MolPick molPick;
  uint totalMolecule;
  uint sourceBox, destBox;
  uint pStartCFCMC, pLenCFCMC;
  uint molIndex, kindIndex;
  uint lambdaWindow, histFreq;
  uint eqCycle;
  vector< vector< vector<uint> > > hist;
  vector< vector<uint> > relaxMolecule;

  double W_tc, W_recip;
  double correct_old, correct_new, self_old, self_new;
  double lambdaMax, nuTolerance;
  double lambdaOld, lambdaNew;
  double relaxRadiusSq;
  double *lambdaRef;
  vector< vector< vector<double> > > bias;
  vector< double > nu;

  //variable needs for relaxing
  uint b, m, mk, pStart, pLen;
  XYZ newCOM;
  XYZArray newMolPos;
  Intermolecular inter_LJ, inter_Real, recip;


  cbmc::TrialMol oldMolCFCMC, newMolCFCMC;
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
  for(uint k = 0; k < molRef.GetKindsCount(); k++) {
    std::cout << "hist " << molRef.kinds[k].name.c_str() << ": ";
    for(uint i = 0; i <= lambdaWindow; i++) {
      std::cout <<  hist[0][k][i] << " ";
    }
    std::cout << std::endl;
  }

  for(uint k = 0; k < molRef.GetKindsCount(); k++) {
    std::cout << "Bias " << molRef.kinds[k].name.c_str() << ": ";
    for(uint i = 0; i <= lambdaWindow; i++) {
      std::cout <<  bias[0][k][i] << " ";
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
    newMolCFCMC = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMolCFCMC = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMolCFCMC.SetCoords(coordCurrRef, pStartCFCMC);
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
    newMolCFCMC.SetCoords(temp, 0);
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


  uint idxNew = GetLambdaIdx(lambdaNew);
  //if lambda is zero, it means molecule successfully  transfered.
  if(idxNew == 0 && calcNewEn) {
    if(ffRef.useLRC) {
      ShiftMolToSourceBox();
      tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
      tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
      W_tc = exp(-1.0 * ffRef.beta * (tcGain.energy + tcLose.energy));
      ShiftMolToDestBox();
    }
  }

  ShiftMolToSourceBox();
  if(calcNewEn) {
    //calculate inter and intra energy before changing lambda
    calcEnRef.SingleMoleculeInter(newEnergy[sourceBox], atomForceNew,
                                  molForceNew, molIndex, sourceBox);
  } else {
    //calculate inter and intra energy after changing lambda
    calcEnRef.SingleMoleculeInter(oldEnergy[sourceBox], atomForceNew,
                                  molForceNew, molIndex, sourceBox);
  }


  ShiftMolToDestBox();
  if(calcNewEn) {
    //calculate inter and intra energy before changing lambda
    calcEnRef.SingleMoleculeInter(newEnergy[destBox], atomForceNew,
                                  molForceNew, molIndex, destBox);
  } else {
    //calculate inter and intra energy after changing lambda
    calcEnRef.SingleMoleculeInter(oldEnergy[destBox], atomForceNew,
                                  molForceNew, molIndex, destBox);
  }

}

inline double CFCMC::GetCoeff() const
{
  uint idxSNew = GetLambdaIdx(lambdaNew);
  uint idxSOld = GetLambdaIdx(lambdaOld);
  uint idxDNew = GetLambdaIdx(1.0 - lambdaNew);
  uint idxDOld = GetLambdaIdx(1.0 - lambdaOld);
  double biasCoef = exp(bias[sourceBox][kindIndex][idxSNew] -
			bias[sourceBox][kindIndex][idxSOld]);
  biasCoef *= exp(bias[destBox][kindIndex][idxDNew] -
		  bias[destBox][kindIndex][idxDOld]);

  //if lambda source is > 0, its not a full molecule
  if(idxSNew > 0) {
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
  bool result = false;
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    uint idxNew = GetLambdaIdx(lambdaNew);
    //If lambdaNew is zero, it means molecule transfered to dest.
    result = (idxNew == 0);
    if(!result) {   
      //Set full interaction in sourceBox, zero interaction in destBox
      lambdaRef[sourceBox * totalMolecule + molIndex] = 1.0;
      lambdaRef[destBox * totalMolecule + molIndex] = 0.0;
      //Shift the molecule back
      ShiftMolToSourceBox();
    }
  } 
  moveSetRef.Update(mv::CFCMC, result, step, destBox, kindIndex);
}


inline void CFCMC::ShiftMolToSourceBox()
{
  cellList.RemoveMol(molIndex, destBox, coordCurrRef);
  //Set coordinates, new COM; shift index to new box's list
  oldMolCFCMC.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
  comCurrRef.SetNew(molIndex, sourceBox);
  molLookRef.ShiftMolBox(molIndex, destBox, sourceBox, kindIndex);
  cellList.AddMol(molIndex, sourceBox, coordCurrRef);
}

inline void CFCMC::ShiftMolToDestBox()
{ 
  cellList.RemoveMol(molIndex, sourceBox, coordCurrRef);
  //Set coordinates, new COM; shift index to new box's list
  newMolCFCMC.GetCoords().CopyRange(coordCurrRef, 0, pStartCFCMC, pLenCFCMC);
  comCurrRef.SetNew(molIndex, destBox);
  molLookRef.ShiftMolBox(molIndex, sourceBox, destBox, kindIndex);
  cellList.AddMol(molIndex, destBox, coordCurrRef);
}

inline void CFCMC::UpdateBias()
{
  if(nu[kindIndex] <= nuTolerance)
    return;

  //Find the index for source and dest box
  uint idxS = GetLambdaIdx(lambdaOld);
  uint idxD = GetLambdaIdx(1.0 - lambdaOld);
#if ENSEMBLE == GCMC
  if(sourceBox == mv::BOX0) {
    hist[sourceBox][kindIndex][idxS] += 1;
    bias[sourceBox][kindIndex][idxS] -= nu[kindIndex];
  } else {
    hist[destBox][kindIndex][idxD] += 1;
    bias[destBox][kindIndex][idxD] -= nu[kindIndex];
  }
#else
  hist[sourceBox][kindIndex][idxS] += 1;
  bias[sourceBox][kindIndex][idxS] -= nu[kindIndex];
  hist[destBox][kindIndex][idxD] += 1;
  bias[destBox][kindIndex][idxD] -= nu[kindIndex];
#endif
  uint box[2];
  box[0] = sourceBox;
  box[1] = destBox;

  for(uint b = 0; b < 2; b++) {
    uint trial = std::accumulate(hist[box[b]][kindIndex].begin(),
				 hist[box[b]][kindIndex].end(), 0);
    if((trial + 1) % histFreq == 0) {
      uint maxVisited = *max_element(hist[box[b]][kindIndex].begin(),
				     hist[box[b]][kindIndex].end());
      uint minVisited = *min_element(hist[box[b]][kindIndex].begin(),
				     hist[box[b]][kindIndex].end());
      //check to see if all the bin is visited atleast 30% of 
      // the most visited bin.                               
      if(minVisited > 0.3 * maxVisited) {
	nu[kindIndex] *= 0.5;
	std::fill_n(hist[box[b]][kindIndex].begin(), lambdaWindow + 1, 0);
      }                               
    }
  }
}

inline void CFCMC::FindRelaxingMolecules(uint box)
{
  //Molecule has to be transfered before calling this function
  relaxMolecule[box].clear();
  if(box >= BOXES_WITH_U_NB) {
    return;
  }

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
  //Start with full interaction in sourceBox, zero interaction in destBox
  lambdaOld = 1.0;
  //pick new lambda
  lambdaNew = lambdaOld + (prng.randInt(1) ? lambdaMax : -lambdaMax);
  UpdateBias();
  while(lambdaNew >= 0.0 && lambdaNew <= 1.0) {
    //Set the interaction in source and destBox 
    lambdaRef[sourceBox * totalMolecule + molIndex] = lambdaOld;
    lambdaRef[destBox * totalMolecule + molIndex] = 1.0 - lambdaOld;
    //Calculate the old energy
    CalcEnCFCMC(false);
    //Update the interaction in sourceBox and destBox
    lambdaRef[sourceBox * totalMolecule + molIndex] = lambdaNew;
    lambdaRef[destBox * totalMolecule + molIndex] = 1.0 - lambdaNew;
    //Calculate the new energy
    CalcEnCFCMC(true);
    bool accepted = AcceptInflating();
    if(accepted) {
      lambdaOld = lambdaNew; 
      RelaxingMolecules();
    }
    UpdateBias();
    uint idxOld = GetLambdaIdx(lambdaOld);
    if(idxOld == 0 || idxOld == lambdaWindow) {
      break;
    }
    //pick new lambda
    lambdaNew = lambdaOld + (prng.randInt(1) ? lambdaMax : -lambdaMax);
  }
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
    swap(atomForceRef, atomForceNew);
    swap(molForceRef, molForceNew);

    //Retotal
    sysPotRef.Total();
  }

  return result;
}


inline void CFCMC::RelaxingMolecules()
{
  long sourceSteps = eqCycle * relaxMolecule[sourceBox].size();
  long destSteps = eqCycle * relaxMolecule[destBox].size();
  ShiftMolToSourceBox();
  if(sourceBox < BOXES_WITH_U_NB) {
    for(uint s = 0; s < sourceSteps; s++) {
      TransformRelaxing(sourceBox);
      CalcEnRelaxing(sourceBox);
      AcceptRelaxing(sourceBox);
    }
  }
  ShiftMolToDestBox();
  if(destBox < BOXES_WITH_U_NB) {
    for(uint s = 0; s < destSteps; s++) {
      TransformRelaxing(destBox);
      CalcEnRelaxing(destBox);
      AcceptRelaxing(destBox);
    }
  }
}

inline void CFCMC::TransformRelaxing(uint b)
{
  uint pickedMol = prng.randIntExc(relaxMolecule[b].size());
  m = relaxMolecule[b][pickedMol];
  mk = molRef.GetMolKind(m);
  pStart = 0, pLen = 0;
  molRef.GetRangeStartLength(pStart, pLen, m);
  newMolPos.Uninit();
  newMolPos.Init(pLen);
  newCOM = comCurrRef.Get(m);
  bool trans = prng.randInt(1);
  trans |= (pLen == 1); // single molecule, translation only

  if(trans) {
    coordCurrRef.TranslateRand(newMolPos, newCOM, pStart, pLen,
                               m, b, moveSetRef.Scale(b, mv::DISPLACE, mk));
  } else {
    coordCurrRef.RotateRand(newMolPos, pStart, pLen, m, b,
                            moveSetRef.Scale(b, mv::ROTATE, mk));
  }
}

inline void CFCMC::CalcEnRelaxing(uint b)
{
  cellList.RemoveMol(m, b, coordCurrRef);
  overlap = false;
  //calculate LJ interaction and real term of electrostatic interaction
  atomForceRef.CopyRange(atomForceNew, 0, 0, atomForceRef.Count());
  molForceRef.CopyRange(molForceNew, 0, 0, molForceRef.Count());
  overlap = calcEnRef.MoleculeInter(inter_LJ, inter_Real, newMolPos, 
                                    atomForceNew, molForceNew, m, b);
  if(!overlap) {
    //calculate reciprocate term of electrostatic interaction
    recip.energy = calcEwald->MolReciprocal(newMolPos, m, b);
  }
}

inline void CFCMC::AcceptRelaxing(uint b)
{
  bool res = false;
  double pr = prng();
  res = pr < exp(-BETA * (inter_LJ.energy + inter_Real.energy +
                          recip.energy));
  bool result = res && !overlap;

  if (result) {
    //Set new energy.
    // setting energy and virial of LJ interaction
    sysPotRef.boxEnergy[b].inter += inter_LJ.energy;
    // setting energy and virial of coulomb interaction
    sysPotRef.boxEnergy[b].real += inter_Real.energy;
    // setting energy and virial of recip term
    sysPotRef.boxEnergy[b].recip += recip.energy;;

    //Copy coords
    newMolPos.CopyRange(coordCurrRef, 0, pStart, pLen);
    comCurrRef.Set(m, newCOM);
    swap(atomForceRef, atomForceNew);
    swap(molForceRef, molForceNew);
    calcEwald->UpdateRecip(b);

    sysPotRef.Total();
  }
  // It means that Recip energy is calculated and move not accepted
  if(!result && !overlap) {
    calcEwald->RestoreMol(m);
  }
  cellList.AddMol(m, b, coordCurrRef);
}


#endif

#endif
