/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCDATA_H
#define DCDATA_H
#include <algorithm>
#include <vector>

#include "BasicTypes.h"
#include "CBMC.h"
#include "Setup.h"
#include "System.h"
#include "XYZArray.h"

class Forcefield;

namespace cbmc {
// Class to avoid reallocating arrays for CBMC
// Could be refactored into an object pool. This would be easier if I had
// bothered to write accessors
class DCData {
public:
  explicit DCData(System &sys, const Forcefield &forcefield, const Setup &set);
  ~DCData();

  const CalculateEnergy &calc;

  const Ewald *calcEwald;

  const Forcefield &ff;
  const BoxDimensions &axes;
  PRNG &prng;

  const uint nAngleTrials;
  const uint nDihTrials;
  const uint nLJTrialsFirst;
  const uint nLJTrialsNth;
  uint totalTrials;

  // used for both angles and dihedrals
  double *angles;
  double *angleWeights;
  double *angleEnergy;

  XYZArray &positions; // candidate positions for inclusion (alias for
                       // multiPositions[0])
  double *inter;       // intermolecule energies, reused for new and old
  double *real;        // short range coulomb interaction
  double *ljWeights;
  double *bonded;
  double *oneFour;
  double *nonbonded;     // calculated nonbonded 1_N LJ and coulomb energies
  double *nonbonded_1_4; // calculated nonbonded 1_4 LJ and coulomb energies
  double *nonbonded_1_3; // calculated nonbonded 1_3 LJ and coulomb energies

  double *interT;     // For DCRotateCOM, we have combined first and Nth trial
  double *realT;      // For DCRotateCOM, we have combined first and Nth trial
  double *ljWeightsT; // For DCRotateCOM, we have combined first and Nth trial
  bool *overlap;      // For detecting overlap for each LJ trial
  bool *overlapT;     // For detecting overlap for each LJ trial. Used in
                      // DCRotateCOM

  XYZArray multiPositions[MAX_BONDS];
};

inline DCData::DCData(System &sys, const Forcefield &forcefield,
                      const Setup &set)
    :

      calc(sys.calcEnergy), ff(forcefield), axes(sys.boxDimRef), prng(sys.prng),
      nAngleTrials(set.config.sys.cbmcTrials.bonded.ang),
      nDihTrials(set.config.sys.cbmcTrials.bonded.dih),
      nLJTrialsFirst(set.config.sys.cbmcTrials.nonbonded.first),
      nLJTrialsNth(set.config.sys.cbmcTrials.nonbonded.nth),
      positions(*multiPositions) {
  calcEwald = sys.GetEwald();
  uint maxLJTrials = nLJTrialsFirst;
  if (nLJTrialsNth > nLJTrialsFirst)
    maxLJTrials = nLJTrialsNth;

  totalTrials = nLJTrialsFirst * nLJTrialsNth;
  if (totalTrials == 0)
    totalTrials = maxLJTrials;

  for (uint i = 0; i < MAX_BONDS; ++i) {
    multiPositions[i] = XYZArray(maxLJTrials);
  }
  inter = new double[maxLJTrials];
  real = new double[maxLJTrials];
  bonded = new double[maxLJTrials];
  oneFour = new double[maxLJTrials];
  nonbonded = new double[maxLJTrials];
  ljWeights = new double[maxLJTrials];
  overlap = new bool[maxLJTrials];

  interT = new double[totalTrials];
  realT = new double[totalTrials];
  ljWeightsT = new double[totalTrials];
  overlapT = new bool[totalTrials];

  uint trialMax = std::max(nAngleTrials, nDihTrials);
  angleEnergy = new double[trialMax];
  angleWeights = new double[trialMax];
  angles = new double[trialMax];
  nonbonded_1_3 = new double[trialMax];
  nonbonded_1_4 = new double[trialMax];
}

inline DCData::~DCData() {
  delete[] inter;
  delete[] real;
  delete[] bonded;
  delete[] oneFour;
  delete[] nonbonded;
  delete[] nonbonded_1_4;
  delete[] nonbonded_1_3;
  delete[] ljWeights;
  delete[] angles;
  delete[] angleWeights;
  delete[] angleEnergy;
  delete[] interT;
  delete[] realT;
  delete[] ljWeightsT;
  delete[] overlap;
  delete[] overlapT;
}

} // namespace cbmc

#endif /*DCDATA_H*/
