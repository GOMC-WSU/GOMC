#include "EnsemblePreprocessor.h"
#include "TrialMol.h"

// An abstract class
class ElectrostaticBase
{
	// Data members of class
public:
	// Pure Virtual Function
	//virtual void show() = 0;

  // calculate self term for a box
  virtual double BoxSelf(uint box) const = 0;

  // calculate correction term for a molecule
  virtual double MolCorrection(uint molIndex, uint box) const = 0;

  // calculate correction term after swap move, with lambda = 1
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol) const = 0;

  // calculate correction term after swap move, with system lambda
  virtual double SwapCorrection(const cbmc::TrialMol &trialMol,
                                const uint molIndex) const = 0;

  // calculate self term after swap move
  virtual double SwapSelf(const cbmc::TrialMol &trialMol) const = 0;

  // It's called in free energy calculation to calculate the change in
  // self energy in all lambda states
  virtual void ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const = 0;

  // It's called in free energy calculation to calculate the change in
  // correction energy in all lambda states
  virtual void ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                                const std::vector<double> &lambda_Coul,
                                const uint iState, const uint molIndex,
                                const uint box) const = 0;


  std::vector<int> particleKind;
  std::vector<int> particleMol;

  // starting index of molecule
  std::vector<int> startMol;
  // length of each molecule
  std::vector<int> lengthMol;
  // starting atom index of each box
  int boxStart[BOX_TOTAL];
  // ending atom index of each box
  int boxEnd[BOX_TOTAL];
  // atom charges
  std::vector<double> particleCharge;
  // which atoms don't have charge
  std::vector<bool> particleHasNoCharge;
/* Other members */
};
