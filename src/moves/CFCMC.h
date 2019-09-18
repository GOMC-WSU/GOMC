/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
      lambdaRef(sys.lambdaRef), MoveBase(sys, statV), MP(sys, statV)
    {
      if(statV.cfcmcVal.enable) {
        MPEnable = statV.cfcmcVal.MPEnable;
        steps = 0;
        histFreq = statV.GetPerAdjust(); 
        relaxSteps = statV.cfcmcVal.relaxSteps;
        lambdaWindow = statV.cfcmcVal.lambdaVDW.size() - 1;
        lambdaCoulomb = statV.cfcmcVal.lambdaCoulomb;
        lambdaVDW = statV.cfcmcVal.lambdaVDW; 
        flatness = statV.cfcmcVal.histFlatness;
        nuTolerance = 1e-6;
        uint totKind = molRef.GetKindsCount();
        nu.resize(BOX_TOTAL);
        hist.resize(BOX_TOTAL);
        bias.resize(BOX_TOTAL);
        kCount.resize(BOX_TOTAL);
        firstPrint.resize(BOX_TOTAL);
        for(uint b = 0; b < BOX_TOTAL; b++) {
          nu[b].resize(totKind, 0.01);
          hist[b].resize(totKind);
          bias[b].resize(totKind);
          kCount[b].resize(totKind);
          firstPrint[b].resize(totKind, true);
        }
        for(uint b = 0; b < BOX_TOTAL; b++) {
          for(uint k = 0; k < totKind; k++) {
            hist[b][k].resize(lambdaWindow + 1, 0);
            bias[b][k].resize(lambdaWindow + 1, 0.0);
          }
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
  void CalcEnCFCMC(uint lambdaIdxOldS, uint lambdaIdxNewS);
  void CalcEnRelaxing(uint box);
  uint TransformRelaxing(uint box);
  void RelaxingMolecules();
  
  Lambda & lambdaRef;
  uint sourceBox, destBox;
  uint pStartCFCMC, pLenCFCMC;
  uint molIndex, kindIndex;
  uint lambdaIdxOld, lambdaIdxNew;
  uint box[2];
  uint relaxSteps, lambdaWindow, histFreq;
  bool overlapCFCMC;
  vector< vector < bool > > firstPrint;
  vector< vector< vector<long int> > > hist;
  vector< vector< uint > > kCount;

  double W_tc, W_recip;
  double correctDiffSource, correctDiffDest, selfDiffSource, selfDiffDest;
  double recipDiffSource, recipDiffDest, tcDiffSource, tcDiffDest;
  double nuTolerance, flatness, molInSourceBox, molInDestBox; 
  vector< vector< vector<double> > > bias;
  vector< double > lambdaCoulomb, lambdaVDW;
  vector< vector<double> > nu;

  //variable needs for relaxing
  MultiParticle MP;
  bool MPEnable;
  uint b, m, mk, pStart, pLen;
  long steps;
  XYZ newCOM;
  XYZArray newMolPos;
  Intermolecular inter_LJ, inter_Real, recip;


  cbmc::TrialMol oldMolCFCMC, newMolCFCMC;
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

  if(state == mv::fail_state::NO_FAIL) {
    pStartCFCMC = pLenCFCMC = 0;
    molRef.GetRangeStartLength(pStartCFCMC, pLenCFCMC, molIndex);
    box[0] = sourceBox;
    box[1] = destBox;
  }
  return state;
}

inline uint CFCMC::Prep(const double subDraw, const double movPerc)
{
  overlapCFCMC = false;
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
    for(uint b = 0; b < BOX_TOTAL; b++) {
      for (uint k = 0; k < molRef.GetKindsCount(); ++k) {
        kCount[b][k] = molLookRef.NumKindInBox(k, b);
        if(b == sourceBox && k == kindIndex) {
          //consider the fraction molecule as different molecule kind
          --kCount[b][k];
        } 
      }
    }
  }
  return state;
}


inline uint CFCMC::Transform()
{
  //Start with full interaction in sourceBox, zero interaction in destBox
  //SInce we have the lambda for growing molecule, in sourceBox, lambdaWindow
  // correspond to full interaction and (lambdaWindow - X) is for destBox
  lambdaIdxOld = lambdaWindow;
  lambdaIdxNew = lambdaIdxOld - 1;
  //Update the interaction in destBox
  lambdaRef.Set(lambdaVDW[lambdaWindow - lambdaIdxNew],
                lambdaCoulomb[lambdaWindow - lambdaIdxNew], molIndex,
                kindIndex, destBox);
  //Start growing the fractional molecule in destBox
  molRef.kinds[kindIndex].BuildIDNew(newMolCFCMC, molIndex);
  overlapCFCMC = newMolCFCMC.HasOverlap();
  //Add bonded energy because we dont considered in DCRotate.cpp 
  newMolCFCMC.AddEnergy(calcEnRef.MoleculeIntra(newMolCFCMC, molIndex));
  ShiftMolToDestBox();

  do{
    //Set the interaction in source and destBox
    lambdaRef.Set(lambdaVDW[lambdaIdxOld], lambdaCoulomb[lambdaIdxOld],
                  molIndex, kindIndex, sourceBox);
    lambdaRef.Set(lambdaVDW[lambdaWindow - lambdaIdxOld],
                 lambdaCoulomb[lambdaWindow - lambdaIdxOld], molIndex,
                 kindIndex, destBox);
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
    CalcEnCFCMC(lambdaIdxOld, lambdaIdxNew);
    //Accept or reject the inflation
    bool acceptedInflate = AcceptInflating();
    if(acceptedInflate) {
      lambdaIdxOld = lambdaIdxNew; 
      //Update the interaction in sourceBox and destBox
      lambdaRef.Set(lambdaVDW[lambdaIdxNew], lambdaCoulomb[lambdaIdxNew],
                    molIndex, kindIndex, sourceBox);
      lambdaRef.Set(lambdaVDW[lambdaWindow - lambdaIdxNew],
                  lambdaCoulomb[lambdaWindow - lambdaIdxNew], molIndex,
                  kindIndex, destBox);
    }
    RelaxingMolecules();
    //Dont update Bias if move resulted in overLap
    if(!overlapCFCMC) {
      UpdateBias();
    }
    //pick new lambda in the neighborhood
    lambdaIdxNew = lambdaIdxOld + (prng.randInt(1) ? 1 : -1);
  } while(lambdaIdxOld > 0 && lambdaIdxOld < lambdaWindow);

  return mv::fail_state::NO_FAIL;
}

inline void CFCMC::CalcEn() {
  return;
}

inline void CFCMC::CalcEnCFCMC(uint lambdaIdxOldS, uint lambdaIdxNewS)
{
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

  if(overlapCFCMC) {
    //Do not calculate the energy difference if we have overlap
    return;
  }
    double lambdaOld_VDW_S = lambdaVDW[lambdaIdxOldS];
    double lambdaNew_VDW_S = lambdaVDW[lambdaIdxNewS];
    double lambdaOld_VDW_D = lambdaVDW[lambdaWindow - lambdaIdxOldS];
    double lambdaNew_VDW_D = lambdaVDW[lambdaWindow - lambdaIdxNewS];
    double lambdaOld_Coulomb_S = lambdaCoulomb[lambdaIdxOldS];
    double lambdaNew_Coulomb_S = lambdaCoulomb[lambdaIdxNewS];
    double lambdaOld_Coulomb_D = lambdaCoulomb[lambdaWindow - lambdaIdxOldS];
    double lambdaNew_Coulomb_D = lambdaCoulomb[lambdaWindow - lambdaIdxNewS];

  //Calculating long range correction
  if(ffRef.useLRC) { 
    //Calculate LRC difference for lambdaNew and lambdaOld
    tcDiffSource = calcEnRef.MoleculeTailChange(sourceBox, kindIndex, 
                                                kCount[sourceBox], lambdaOld_VDW_S, 
                                                lambdaNew_VDW_S);
    tcDiffDest = calcEnRef.MoleculeTailChange(destBox, kindIndex,
					                                    kCount[destBox], lambdaOld_VDW_D,
                                              lambdaNew_VDW_D);
    W_tc = exp(-1.0 * ffRef.beta * (tcDiffSource + tcDiffDest));
  }

  //No need to calculate energy when performing CBMC
  ShiftMolToSourceBox();
  if(lambdaIdxNew != 0) {
    //calculate inter energy for lambda new and old in source Box
    calcEnRef.SingleMoleculeInter(oldEnergy[sourceBox], newEnergy[sourceBox],
                                  lambdaOld_VDW_S,lambdaNew_VDW_S,
                                  lambdaOld_Coulomb_S, lambdaNew_Coulomb_S,
                                  molIndex, sourceBox);
  }


  ShiftMolToDestBox();
  if(lambdaIdxOld != lambdaWindow && lambdaIdxNew != lambdaWindow) {
    //calculate inter energy for lambda new and old in dest Box
    calcEnRef.SingleMoleculeInter(oldEnergy[destBox], newEnergy[destBox],
                                  lambdaOld_VDW_D, lambdaNew_VDW_D,
                                  lambdaOld_Coulomb_D, lambdaNew_Coulomb_D, 
                                  molIndex, destBox);
  }


  //Calculate self and correction difference for lambdaNew and lambdaOld
  //For electrostatic we use linear scaling
  double coefDiffS = lambdaNew_Coulomb_S - lambdaOld_Coulomb_S;
  double coefDiffD = lambdaNew_Coulomb_D - lambdaOld_Coulomb_D;
  correctDiffSource = coefDiffS * calcEwald->SwapCorrection(oldMolCFCMC);
  correctDiffDest = coefDiffD * calcEwald->SwapCorrection(newMolCFCMC);
  selfDiffSource = coefDiffS * calcEwald->SwapSelf(oldMolCFCMC);
  selfDiffDest = coefDiffD * calcEwald->SwapSelf(newMolCFCMC);
  //calculate Recprocal Difference in source and dest box
  recipDiffSource = calcEwald->CFCMCRecip(oldMolCFCMC.GetCoords(),
                                          lambdaOld_Coulomb_S,
                                          lambdaNew_Coulomb_S, 
                                          molIndex, sourceBox);
  recipDiffDest = calcEwald->CFCMCRecip(newMolCFCMC.GetCoords(), 
                                        lambdaOld_Coulomb_D, 
                                        lambdaNew_Coulomb_D, 
                                        molIndex, destBox);

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
  double biasCoef = 1.0;
  biasCoef *= exp(bias[sourceBox][kindIndex][idxSNew] - 
		  bias[sourceBox][kindIndex][idxSOld]);
  biasCoef *= exp(bias[destBox][kindIndex][idxDNew] -
		  bias[destBox][kindIndex][idxDOld]);

  //transition from lambda = 1 to any value means deletion
  //transition to lambda = 1 from any value means insertion
#if ENSEMBLE == GEMC
  if(idxSOld == lambdaWindow) {
    coef *= molInSourceBox * boxDimRef.volInv[sourceBox];
    coef *= 0.5;
  }
  if(idxSNew == lambdaWindow) { 
    coef *= boxDimRef.volume[sourceBox] / molInSourceBox;
    coef *= 2.0;
  } else if(idxSNew == 0) {
    //same as (idxDNew == lambdaWindow)
    coef *= boxDimRef.volume[destBox] / (molInDestBox + 1.0);
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
  steps = step;
  //If we didn't skip the move calculation
  if(rejectState == mv::fail_state::NO_FAIL) {
    //If lambdaIdxOld is zero, it means molecule transfered to destBox.
    result = (lambdaIdxOld == 0);
    if(result) {
      //Set full interaction in destBox, zero interaction in sourceBox
      lambdaRef.UnSet(destBox, sourceBox);
    } else {   
      //Set full interaction in sourceBox, zero interaction in destBox
      lambdaRef.UnSet(sourceBox, destBox);
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
  //Find the index for source and dest box
  uint idxS = lambdaIdxOld;
  uint idxD = lambdaWindow - lambdaIdxOld;
  vector < uint > idx;
  idx.push_back(idxS);
  idx.push_back(idxD);

  for(uint b = 0; b < 2; b++) { 
    hist[box[b]][kindIndex][idx[b]] += 1;

    //In Bias, if lambda is 1.0, it is also 0.0 for continueity
    if((idxS == lambdaWindow) || (idxS == 0)) {   
      hist[box[b]][kindIndex][idx[1 - b]] += 1;
    }

    //Stop the modifying bias if we converged
    if(nu[box[b]][kindIndex] <= nuTolerance) {
      if(firstPrint[box[b]][kindIndex]) {
        printf("STOPED MODIFYING BIAS FOR %s. \n",
        molRef.kinds[kindIndex].name.c_str());
        firstPrint[box[b]][kindIndex] = false;
      }
      //Check after equilibration if all the bin visited atleast X% of 
      // the most visited bin
      long trial = std::accumulate(hist[box[b]][kindIndex].begin(), hist[box[b]][kindIndex].end(), 0);
      if((trial + 1) % (histFreq * 10) == 0) {
        long int maxVisited = *max_element(hist[box[b]][kindIndex].begin(),
				   hist[box[b]][kindIndex].end());
        long int minVisited = *min_element(hist[box[b]][kindIndex].begin(),
              hist[box[b]][kindIndex].end()); 
        //check to see if all the bin visited atleast X% of the most visited bin
        if(minVisited < flatness * maxVisited) {
          //nu[box[b]][kindIndex] = 0.01;
          //std::fill_n(hist[box[b]][kindIndex].begin(), lambdaWindow + 1, 0);
          printf("Warning [%d][%4s]: Minimum visited lambda state (%ld) is not %4.2f of Maximum visited lambda state (%ld)! \n",
                box[b], molRef.kinds[kindIndex].name.c_str(), minVisited,
                flatness, maxVisited);
        } 
      }
      continue;
    } 

    //Update the bias for both box
    bias[box[b]][kindIndex][idx[b]] -= nu[box[b]][kindIndex];
    //In Bias, if lambda is 1.0, it is also 0.0 for continueity
    if((idxS == lambdaWindow) || (idxS == 0)) {
      bias[box[b]][kindIndex][idx[1 - b]] -= nu[box[b]][kindIndex];
    }

#if ENSEMBLE == GCMC
    //We dont consider biasing in reservoir
    bias[1][kindIndex][idxD] = bias[1][kindIndex][idxS] = 0.0;
    hist[1][kindIndex][idxD] = hist[1][kindIndex][idxS] = 0;
#endif

    //long trial = moveSetRef.GetTrial(box[b], mv::CFCMC, kindIndex);
    long trial = std::accumulate(hist[box[b]][kindIndex].begin(), hist[box[b]][kindIndex].end(), 0);
    //Check flatness 
    if((trial + 1) % histFreq == 0) {
      uint maxVisited = *max_element(hist[box[b]][kindIndex].begin(),
				   hist[box[b]][kindIndex].end());
      uint minVisited = *min_element(hist[box[b]][kindIndex].begin(),
            hist[box[b]][kindIndex].end()); 
      //check to see if all the bin visited atleast X% of the most visited bin
      if(minVisited > flatness * maxVisited) {
        nu[box[b]][kindIndex] *= 0.5;
        std::fill_n(hist[box[b]][kindIndex].begin(), lambdaWindow + 1, 0);
        printf("Bias-Adj [%d][%4s]: %4.10f \n",
               box[b], molRef.kinds[kindIndex].name.c_str(),
               nu[box[b]][kindIndex]);
      } else {
        printf("Hist     [%d][%4s]: [", box[b],
              molRef.kinds[kindIndex].name.c_str());
        for(uint i = 0; i <= lambdaWindow; i++) {
          printf("%8ld", hist[box[b]][kindIndex][i]);
        }
        std::cout << "] \n";
        //Reset the histogram to reevaluate it
        //std::fill_n(hist[box[b]][kindIndex].begin(), lambdaWindow + 1, 0);
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
/*
  if(!overlapCFCMC) {
    printf("lambda[%5s]: %-2d -> %-2d : WS: %5.1e : WD: %5.1e : Wtot: %5.1e : Coef: %5.1e : SourceBox: %d \n", molRef.kinds[kindIndex].name.c_str(), lambdaIdxOld, lambdaIdxNew, W1, W2, Wrat, molTransCoeff, sourceBox);
  }
*/
  if(result) {
    //Add tail corrections
    sysPotRef.boxEnergy[sourceBox].tc += tcDiffSource;
    sysPotRef.boxEnergy[destBox].tc += tcDiffDest;
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

    calcEwald->UpdateRecip(sourceBox);
    calcEwald->UpdateRecip(destBox);

    //Retotal
    sysPotRef.Total();
    //set single move accept to true for multiparticle
    moveSetRef.SetSingleMoveAccepted();
  } 
  overlapCFCMC = false;

  return result;
}


inline void CFCMC::RelaxingMolecules()
{
  ShiftMolToSourceBox();
  if(sourceBox < BOXES_WITH_U_NB) {
    for(uint s = 0; s < relaxSteps; s++) {
      uint state = TransformRelaxing(sourceBox);
      if(state == mv::fail_state::NO_FAIL) {
        CalcEnRelaxing(sourceBox);
        AcceptRelaxing(sourceBox);
      }
    }
    oldMolCFCMC.SetCoords(coordCurrRef, pStartCFCMC);
  }

  ShiftMolToDestBox();
  if(destBox < BOXES_WITH_U_NB) {
    for(uint s = 0; s < relaxSteps; s++) {
      uint state = TransformRelaxing(destBox);
      if(state == mv::fail_state::NO_FAIL) {
        CalcEnRelaxing(destBox);
        AcceptRelaxing(destBox);
      }
    }
    newMolCFCMC.SetCoords(coordCurrRef, pStartCFCMC);
  }
}

inline uint CFCMC::TransformRelaxing(uint b)
{
  uint state = mv::fail_state::NO_FAIL;
  if(MPEnable) {
    MP.PrepCFCMC(b);
    MP.Transform();
  } else {
    //Randomely pick a molecule in Box
    uint state = prng.PickMol(m, mk, b);
    if(state == mv::fail_state::NO_FAIL) {
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
  }
  return state;
}

inline void CFCMC::CalcEnRelaxing(uint b)
{
  if(MPEnable) {
    MP.CalcEn();
  } else {
    cellList.RemoveMol(m, b, coordCurrRef);
    overlap = false;
    //calculate LJ interaction and real term of electrostatic interaction
    overlap = calcEnRef.MoleculeInter(inter_LJ, inter_Real, newMolPos, m, b);
    if(!overlap) {
      //calculate reciprocate term of electrostatic interaction
      recip.energy = calcEwald->MolReciprocal(newMolPos, m, b);
    }
  }
}

inline void CFCMC::AcceptRelaxing(uint b)
{
  if(MPEnable) {
    MP.Accept(mv::fail_state::NO_FAIL, steps);
  } else {
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
      calcEwald->UpdateRecip(b);

      sysPotRef.Total();
    }
    cellList.AddMol(m, b, coordCurrRef);
  }
}


#endif

#endif
