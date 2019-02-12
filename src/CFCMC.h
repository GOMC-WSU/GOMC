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
      relaxSteps = statV.cfcmcVal.relaxSteps;
      lambdaWindow = statV.cfcmcVal.window;
      histFreq = statV.cfcmcVal.updateBiasFreq;
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
  bool AcceptInflating();
  void AcceptRelaxing(uint box);
  void CalcEnCFCMC(double lambdaOldS, double lambdaNewS);
  void CalcEnRelaxing(uint box);
  void TransformRelaxing(uint box);
  void RelaxingMolecules();
  
  uint totalMolecule;
  uint sourceBox, destBox;
  uint pStartCFCMC, pLenCFCMC;
  uint molIndex, kindIndex;
  uint lambdaWindow, histFreq;
  uint lambdaIdxOld, lambdaIdxNew;
  uint relaxSteps;
  bool overlapCFCMC;
  vector< vector< vector<long int> > > hist;

  double W_tc, W_recip;
  double correctDiffSource, correctDiffDest, selfDiffSource, selfDiffDest;
  double recipDiffSource, recipDiffDest;
  double lambdaMax, nuTolerance;
  double *lambdaRef;
  double molInSourceBox, molInDestBox;  
  vector< vector< vector<double> > > bias;
  vector< double > nu;

  //variable needs for relaxing
  uint b, m, mk, pStart, pLen;
  XYZ newCOM;
  XYZArray newMolPos;
  Intermolecular inter_LJ, inter_Real, recip;
  double sDraw, mPerc;


  cbmc::TrialMol oldMolCFCMC, newMolCFCMC;
  Intermolecular tcLose, tcGain;
  Energy oldEnergy[BOX_TOTAL], newEnergy[BOX_TOTAL];
  MoleculeLookup & molLookRef;
  Forcefield const& ffRef;
  SystemPotential sysPotNew;
};

void CFCMC::PrintAcceptKind() {
  for(uint k = 0; k < molRef.GetKindsCount(); k++) {
    printf("%-30s %-5s ", "% Accepted CFCMC-Transfer ",
	   molRef.kinds[k].name.c_str());
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
  overlapCFCMC = false;
  sDraw = subDraw;
  mPerc = movPerc;
  uint state = GetBoxPairAndMol(subDraw, movPerc);
  if(state == mv::fail_state::NO_FAIL) {
    newMolCFCMC = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, destBox);
    oldMolCFCMC = cbmc::TrialMol(molRef.kinds[kindIndex], boxDimRef, sourceBox);
    oldMolCFCMC.SetCoords(coordCurrRef, pStartCFCMC);
   //Unwrap the old coordinate for using in new coordinate after wraping
    XYZArray mol(pLenCFCMC); 
    coordCurrRef.CopyRange(mol, pStartCFCMC, 0, pLenCFCMC); 
    boxDimRef.UnwrapPBC(mol, sourceBox, comCurrRef.Get(molIndex));
    boxDimRef.WrapPBC(mol, destBox); 
    //Later it will shift to random COM
    newMolCFCMC.SetCoords(mol, 0);
    // store number of molecule in box before shifting molecule
    molInSourceBox = (double)molLookRef.NumKindInBox(kindIndex, sourceBox);
    molInDestBox = (double)molLookRef.NumKindInBox(kindIndex, destBox);
  }
  return state;
}


inline uint CFCMC::Transform()
{
  //Start with full interaction in sourceBox, zero interaction in destBox
  lambdaIdxOld = lambdaWindow;
  lambdaIdxNew = lambdaIdxOld - 1;
  double lambdaOld = (double)(lambdaIdxOld) * lambdaMax;
  double lambdaNew = (double)(lambdaIdxNew) * lambdaMax;
  //Update the interaction in destBox
  lambdaRef[destBox * totalMolecule + molIndex] = 1.0 - lambdaNew;
  //Start growing the fractional molecule in destBox
  molRef.kinds[kindIndex].BuildIDNew(newMolCFCMC, molIndex);
  overlapCFCMC = newMolCFCMC.HasOverlap();
  //Add bonded energy because we dont considered in DCRotate.cpp 
  newMolCFCMC.AddEnergy(calcEnRef.MoleculeIntra(newMolCFCMC, molIndex));
  ShiftMolToDestBox();
  UpdateBias();

  do{
    //Set the interaction in source and destBox 
    lambdaRef[sourceBox * totalMolecule + molIndex] = lambdaOld;
    lambdaRef[destBox * totalMolecule + molIndex] = 1.0 - lambdaOld;
    if(lambdaIdxNew == 0) {
      //removing the fractional molecule in last steps using CBMC in sourceBox
      molRef.kinds[kindIndex].BuildIDOld(oldMolCFCMC, molIndex);
      //Add bonded energy because we dont considered in DCRotate.cpp 
      oldMolCFCMC.AddEnergy(calcEnRef.MoleculeIntra(oldMolCFCMC, molIndex));
    } else if(lambdaIdxNew == lambdaWindow) {
      //removing the inserted fractional molecule using CBMC in destBox
      cellList.RemoveMol(molIndex, destBox, coordCurrRef);
      molRef.kinds[kindIndex].BuildIDOld(newMolCFCMC, molIndex);
      //Add bonded energy because we dont considered in DCRotate.cpp 
      newMolCFCMC.AddEnergy(calcEnRef.MoleculeIntra(newMolCFCMC, molIndex));
      cellList.AddMol(molIndex, destBox, coordCurrRef);
    }
    //Calculate the old and new energy in source and destBox(if we dont do CBMC)
    CalcEnCFCMC(lambdaOld, lambdaNew);
    //Accept or reject the inflation
    bool acceptedInflate = AcceptInflating();
    if(acceptedInflate) {
      lambdaIdxOld = lambdaIdxNew; 
      //Update the interaction in sourceBox and destBox
      lambdaRef[sourceBox * totalMolecule + molIndex] = lambdaNew;
      lambdaRef[destBox * totalMolecule + molIndex] = 1.0 - lambdaNew;
    }
    RelaxingMolecules();
    UpdateBias();
    //pick new lambda in the neighborhood
    lambdaIdxNew = lambdaIdxOld + (prng.randInt(1) ? 1 : -1);
    lambdaOld = (double)(lambdaIdxOld) * lambdaMax;
    lambdaNew = (double)(lambdaIdxNew) * lambdaMax;
  } while(lambdaIdxOld > 0 && lambdaIdxOld < lambdaWindow);

  return mv::fail_state::NO_FAIL;
}

inline void CFCMC::CalcEn() {
  return;
}

inline void CFCMC::CalcEnCFCMC(double lambdaOldS, double lambdaNewS)
{
  W_tc = 1.0;
  W_recip = 1.0;
  correctDiffSource = 0.0;
  correctDiffDest = 0.0;
  selfDiffSource = 0.0;
  selfDiffDest = 0.0;
  tcLose.Zero();
  tcGain.Zero();

  if(overlapCFCMC) {
    //Do not calculate the energy difference if we have overlap
    return;
  }

  //Calculating long range correction
  if(ffRef.useLRC) { 
    if(sourceBox == mv::BOX0) {
      //Deletion move
      if(lambdaIdxOld == lambdaWindow) {
	//tansition from lambda 1.0 to lower value	
	ShiftMolToSourceBox();
	tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
	tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
	W_tc = exp(-1.0 * ffRef.beta * (tcGain.energy + tcLose.energy));
	ShiftMolToDestBox();
      } else if(lambdaIdxNew == lambdaWindow) {
	//tansition from lambda lower value ro 1.0	
	tcLose = calcEnRef.MoleculeTailChange(destBox, kindIndex, false);
	tcGain = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, true);
	W_tc = exp(-1.0 * ffRef.beta * (tcGain.energy + tcLose.energy));
      }
    } else {
      //Insertion move
      if(lambdaIdxNew == 0) {
	//transition from lambda lower to 1.0 in destBox
	ShiftMolToSourceBox();
	tcLose = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, false);
	tcGain = calcEnRef.MoleculeTailChange(destBox, kindIndex, true);
	W_tc = exp(-1.0 * ffRef.beta * (tcGain.energy + tcLose.energy));
	ShiftMolToDestBox();
      }
    }
  }

  //No need to calculate energy when performing CBMC
  ShiftMolToSourceBox();
  oldEnergy[sourceBox].Zero();
  newEnergy[sourceBox].Zero();
  if(lambdaIdxNew != 0) {
    //calculate inter energy for lambda new and old in source Box
    calcEnRef.SingleMoleculeInter(oldEnergy[sourceBox], newEnergy[sourceBox],
				  atomForceNew, molForceNew, lambdaOldS,
				  lambdaNewS, molIndex, sourceBox);
  }


  ShiftMolToDestBox();
  oldEnergy[destBox].Zero();
  newEnergy[destBox].Zero();
  if(lambdaIdxOld != lambdaWindow && lambdaIdxNew != lambdaWindow) {
    //calculate inter energy for lambda new and old in dest Box
    calcEnRef.SingleMoleculeInter(oldEnergy[destBox], newEnergy[destBox],
				  atomForceNew, molForceNew, 1.0 - lambdaOldS,
				  1.0 - lambdaNewS, molIndex, destBox);
  }

  //Calculate self and correction difference for lambdaNew and lambdaOld
  //For electrostatic we use lamda**5 
  double coefDiffS = pow(lambdaNewS, 5) - pow(lambdaOldS, 5);
  double coefDiffD = pow(1.0 - lambdaNewS, 5) - pow(1.0 - lambdaOldS, 5);
  correctDiffSource = coefDiffS * calcEwald->SwapCorrection(oldMolCFCMC);
  correctDiffDest = coefDiffD * calcEwald->SwapCorrection(newMolCFCMC);
  selfDiffSource = coefDiffS * calcEwald->SwapSelf(oldMolCFCMC);
  selfDiffDest = coefDiffD * calcEwald->SwapSelf(newMolCFCMC);
  //calculate Recprocal Difference in source and dest box
  recipDiffSource = calcEwald->CFCMCRecip(oldMolCFCMC.GetCoords(), lambdaOldS,
					  lambdaNewS, molIndex, sourceBox);
  recipDiffDest = calcEwald->CFCMCRecip(newMolCFCMC.GetCoords(), 1.0-lambdaOldS,
					1.0 - lambdaNewS, molIndex, destBox);

  //need to contribute the self and correction energy
  W_recip = exp(-1.0 * ffRef.beta * (recipDiffSource + recipDiffDest +
				     correctDiffSource + correctDiffDest +
				     selfDiffSource + selfDiffDest));
}

inline double CFCMC::GetCoeff() const
{
  double coef = 1.0;
  uint idxSNew = lambdaIdxNew;
  uint idxSOld = lambdaIdxOld;
  uint idxDNew = lambdaWindow - lambdaIdxNew;
  uint idxDOld = lambdaWindow - lambdaIdxOld;;
  double biasCoef = exp(bias[sourceBox][kindIndex][idxSNew] -
			bias[sourceBox][kindIndex][idxSOld]);
  biasCoef *= exp(bias[destBox][kindIndex][idxDNew] -
		  bias[destBox][kindIndex][idxDOld]);
  //biasCoef = 1.0;

  //if lambda source is > 0, its not a full molecule
#if ENSEMBLE == GEMC
  if(idxSOld == lambdaWindow) {
    coef *= molInSourceBox / (molInDestBox + 1.0);
    coef *= boxDimRef.volume[destBox] * boxDimRef.volInv[sourceBox];
    coef *= 0.5;
  }
  if(idxSNew == lambdaWindow) {  
    coef *= (molInDestBox + 1.0) / molInSourceBox;
    coef *= boxDimRef.volInv[destBox] * boxDimRef.volume[sourceBox];
    coef *= 2.0;
  } else if(idxSNew == 0) {
    coef *= 2.0;
  }

  return coef * biasCoef;
  
#elif ENSEMBLE == GCMC
  if(sourceBox == mv::BOX0) {
    //deletion 
    if(idxSOld == lambdaWindow) {
      coef *= molInSourceBox * boxDimRef.volInv[sourceBox];
      coef *= exp(-BETA * molRef.kinds[kindIndex].chemPot);
      coef *= 0.5;
    }
    if(idxSNew == lambdaWindow) {  
      coef *= boxDimRef.volume[sourceBox] / molInSourceBox;
      coef *= exp(BETA * molRef.kinds[kindIndex].chemPot);
      coef *= 2.0;
    } else if(idxSNew == 0) {
      coef *= 2.0;
    }
    return coef * biasCoef;   
  } else {
    //insertion
    if(idxDOld == 0) {
      coef *= 0.5;
    } 
    if(idxDNew == 0) {
      coef *= 2.0;
    } else if(idxDNew == lambdaWindow) {  
      coef *= boxDimRef.volume[destBox] / (molInDestBox + 1.0);
      coef *= exp(BETA * molRef.kinds[kindIndex].chemPot);
      coef *= 2.0;
    }
    return coef * biasCoef; 
  } 
#endif

}

inline void CFCMC::Accept(const uint rejectState, const uint step)
{
  bool result = false;
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    //If lambdaIdxOld is zero, it means molecule transfered to destBox.
    result = (lambdaIdxOld == 0);
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
  uint idxS = lambdaIdxOld;
  uint idxD = lambdaWindow - lambdaIdxOld;
#if ENSEMBLE == GCMC
  if(sourceBox == mv::BOX0) {
    hist[sourceBox][kindIndex][idxS] += 1;
    //bias[sourceBox][kindIndex][idxS] -= nu[kindIndex];
  } else {
    hist[destBox][kindIndex][idxD] += 1;
    //bias[destBox][kindIndex][idxD] -= nu[kindIndex];
  }
#else
  hist[sourceBox][kindIndex][idxS] += 1;
  //bias[sourceBox][kindIndex][idxS] -= nu[kindIndex];
  hist[destBox][kindIndex][idxD] += 1;
  //bias[destBox][kindIndex][idxD] -= nu[kindIndex];
#endif
  uint box[2];
  box[0] = sourceBox;
  box[1] = destBox;

  printf("lambdaIdxOld: %d, lambdaIdxNew: %d\n", lambdaIdxOld, lambdaIdxNew);

  for(uint b = 0; b < 2; b++) {
    //uint trial = std::accumulate(hist[box[b]][kindIndex].begin(), hist[box[b]][kindIndex].end(), 0);
    uint trial = moveSetRef.GetTrial(box[b], mv::CFCMC, kindIndex);
    if((trial + 1) % histFreq == 0) {
      uint maxVisited = *max_element(hist[box[b]][kindIndex].begin(),
				     hist[box[b]][kindIndex].end());
      uint minVisited = *min_element(hist[box[b]][kindIndex].begin(),
				     hist[box[b]][kindIndex].end());

      for(uint i = 0; i <= lambdaWindow; i++) {
	if(hist[box[b]][kindIndex][i] == 0) {
	  bias[box[b]][kindIndex][i] = 1.0;
	} else {
	  bias[box[b]][kindIndex][i] =
	    log(1.0 / (double)(hist[box[b]][kindIndex][i]));
	}
      }
      
      //check to see if all the bin is visited atleast 95% of 
      // the most visited bin.                               
      if(minVisited > 0.95 * maxVisited) {
	nu[kindIndex] *= 0.5;
	std::fill_n(hist[box[b]][kindIndex].begin(), lambdaWindow + 1, 0);
      }                               
    }
  }
}


inline bool CFCMC::AcceptInflating()
{
  double molTransCoeff = GetCoeff();
  double W1 = exp(-BETA * (newEnergy[sourceBox].Total() -
                          oldEnergy[sourceBox].Total()));
  double W2 = exp(-BETA * (newEnergy[destBox].Total() -
                          oldEnergy[destBox].Total()));
  //override the weight and energy if we used CBMC.
  if(lambdaIdxNew == 0) {
    W1 = 1.0 / oldMolCFCMC.GetWeight();
    oldEnergy[sourceBox] = oldMolCFCMC.GetEnergy();
    oldMolCFCMC.Reset();
  } else if(lambdaIdxNew == lambdaWindow) {
    W2 = 1.0 / newMolCFCMC.GetWeight();
    oldEnergy[destBox] = newMolCFCMC.GetEnergy();
    newMolCFCMC.Reset();
  }
  if(lambdaIdxOld == lambdaWindow) {
    W2 = newMolCFCMC.GetWeight();
    newEnergy[destBox] = newMolCFCMC.GetEnergy();
    newMolCFCMC.Reset();
  }
  double Wrat = W1 * W2 * W_tc * W_recip;

  bool result = prng() < molTransCoeff * Wrat;
  //Reject the move if we had overlaped
  result = result && !overlapCFCMC;

  //printf("DestBox: %d, Result: %d, lambdaIDOld: %d, lambdaIDNew: %d, W1: %f, W2: %f,Wtot: %f \n", destBox, result,lambdaIdxOld, lambdaIdxNew, W1, W2, Wrat);

  if(result) {
    //Add tail corrections
    sysPotRef.boxEnergy[sourceBox].tc += tcLose.energy;
    sysPotRef.boxEnergy[destBox].tc += tcGain.energy;
    //Add rest of energy.
    sysPotRef.boxEnergy[sourceBox] -= oldEnergy[sourceBox];
    sysPotRef.boxEnergy[sourceBox] += newEnergy[sourceBox];
    sysPotRef.boxEnergy[destBox] -= oldEnergy[destBox];
    sysPotRef.boxEnergy[destBox] += newEnergy[destBox];
    //Add correction energy
    sysPotRef.boxEnergy[sourceBox].correction += correctDiffSource;
    sysPotRef.boxEnergy[destBox].correction += correctDiffDest;
    //Add self energy
    sysPotRef.boxEnergy[sourceBox].self += selfDiffSource;
    sysPotRef.boxEnergy[destBox].self += selfDiffDest;
    //Add Reciprocal energy
    sysPotRef.boxEnergy[sourceBox].recip += recipDiffSource;
    sysPotRef.boxEnergy[destBox].recip += recipDiffDest;


    swap(atomForceRef, atomForceNew);
    swap(molForceRef, molForceNew);

    calcEwald->UpdateRecip(sourceBox);
    calcEwald->UpdateRecip(destBox);

    //Retotal
    sysPotRef.Total();
  }

  return result;
}


inline void CFCMC::RelaxingMolecules()
{
  ShiftMolToSourceBox();
  if(sourceBox < BOXES_WITH_U_NB) {
    for(uint s = 0; s < relaxSteps; s++) {
      TransformRelaxing(sourceBox);
      CalcEnRelaxing(sourceBox);
      AcceptRelaxing(sourceBox);
    }
    oldMolCFCMC.SetCoords(coordCurrRef, pStartCFCMC);
  }

  ShiftMolToDestBox();
  if(destBox < BOXES_WITH_U_NB) {
    for(uint s = 0; s < relaxSteps; s++) {
      TransformRelaxing(destBox);
      CalcEnRelaxing(destBox);
      AcceptRelaxing(destBox);
    }
    newMolCFCMC.SetCoords(coordCurrRef, pStartCFCMC);
  }
}

inline void CFCMC::TransformRelaxing(uint b)
{
  //Randomely pick a molecule in Box
  prng.PickMol(m, mk, b, sDraw, mPerc);
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
