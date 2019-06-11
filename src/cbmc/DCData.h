/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCDATA_H
#define DCDATA_H
#include "BasicTypes.h"
#include "XYZArray.h"
#include "Setup.h"
#include "System.h"
#include "CBMC.h"
#include <vector>
#include <algorithm>

class Forcefield;


namespace cbmc
{
//Class to avoid reallocating arrays for CBMC
//Could be refactored into an object pool. This would be easier if I had bothered to write accessors
class DCData
{
public:
  explicit  DCData(System& sys, const Forcefield& forcefield,
                   const Setup& set);
  ~DCData();

  const CalculateEnergy& calc;

  const Ewald  *calcEwald;

  const Forcefield& ff;
  const BoxDimensions& axes;
  PRNG& prng;

  const uint nAngleTrials;
  const uint nDihTrials;
  const uint nLJTrialsFirst;
  const uint nLJTrialsNth;
  uint totalTrials;

  //used for both angles and dihedrals
  real* angles;
  real* angleWeights;
  real* angleEnergy;

  XYZArray& positions;     //candidate positions for inclusion (alias for multiPositions[0])
  real* inter;          //intermolecule energies, reused for new and old
  real* real_en;           //short range coulomb interaction
  real* ljWeights;
  real* bonded;
  real* oneFour;
  real* nonbonded;      //calculated nonbonded 1_N LJ and coulomb energie
  real* nonbonded_1_4;  //calculated nonbonded 1_4 LJ and coulomb energie
  real* nonbonded_1_3;  //calculated nonbonded 1_3 LJ and coulomb energie

  real* interT;     //For DCRotateCOM, we have combined first and Nth trial
  real* realT;      //For DCRotateCOM, we have combined first and Nth trial
  real* ljWeightsT; //For DCRotateCOM, we have combined first and Nth trial
  bool* overlap;      //For detecting overlap for each LJ trial
  bool* overlapT;     //For detecting overlap for each LJ trial. Used in DCRotateCOM

  XYZArray multiPositions[MAX_BONDS];
};

inline DCData::DCData(System& sys, const Forcefield& forcefield, const Setup& set):

  calc(sys.calcEnergy), ff(forcefield),
  prng(sys.prng), axes(sys.boxDimRef),
  nAngleTrials(set.config.sys.cbmcTrials.bonded.ang),
  nDihTrials(set.config.sys.cbmcTrials.bonded.dih),
  nLJTrialsFirst(set.config.sys.cbmcTrials.nonbonded.first),
  nLJTrialsNth(set.config.sys.cbmcTrials.nonbonded.nth),
  positions(*multiPositions)
{
  calcEwald = sys.GetEwald();
  uint maxLJTrials = nLJTrialsFirst;
  if ( nLJTrialsNth > nLJTrialsFirst )
    maxLJTrials = nLJTrialsNth;

  totalTrials = nLJTrialsFirst * nLJTrialsNth;
  if(totalTrials == 0)
    totalTrials = maxLJTrials;

  for(uint i = 0; i < MAX_BONDS; ++i) {
    multiPositions[i] = XYZArray(maxLJTrials);
  }
  inter = new real[maxLJTrials];
  real_en = new real[maxLJTrials];
  bonded = new real[maxLJTrials];
  oneFour = new real[maxLJTrials];
  nonbonded = new real[maxLJTrials];
  ljWeights = new real[maxLJTrials];
  overlap = new bool[maxLJTrials];

  interT = new real[totalTrials];
  realT = new real[totalTrials];
  ljWeightsT = new real[totalTrials];
  overlapT = new bool[totalTrials];

  uint trialMax = std::max(nAngleTrials, nDihTrials);
  angleEnergy = new real[trialMax];
  angleWeights = new real[trialMax];
  angles = new real[trialMax];
  nonbonded_1_3 = new real[trialMax];
  nonbonded_1_4 = new real[trialMax];
}

inline DCData::~DCData()
{
  delete[] inter;
  delete[] real_en;
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

}

#endif
