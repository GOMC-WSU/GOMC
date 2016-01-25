/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCDATA_H
#define DCDATA_H
#include "../../lib/BasicTypes.h"
#include "../XYZArray.h"
#include "../Setup.h"
#include "../System.h"
#include "../CBMC.h"
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
      const Forcefield& ff;
      const BoxDimensions& axes;
      PRNG& prng;

      const uint nAngleTrials;
      const uint nDihTrials;
      const uint nLJTrialsFirst;
      const uint nLJTrialsNth;

      //used for both angles and dihedrals
      double* angles;
      double* angleWeights;
      double* angleEnergy;

      XYZArray& positions;     //candidate positions for inclusion (alias for multiPositions[0])
      double* inter;          //intermolecule energies, reused for new and old
      double* bonded;
      double* nonbonded;      //calculated nonbonded 1_N energies
      double* nonbonded_1_4;  //calculated nonbonded 1_4 energies and 1_3 in
                              //case of Martini forcefield 
	  double* self;			//energies for Ewald calculations
	  double* correction;
	  double* real;
      double* ljWeights;

      XYZArray multiPositions[MAX_BONDS];
   };

inline DCData::DCData(System& sys, const Forcefield& forcefield, const Setup& set):
  calc(sys.calcEnergy), ff(forcefield), prng(sys.prng), axes(sys.boxDimRef),
  nAngleTrials(set.config.sys.cbmcTrials.bonded.ang),
  nDihTrials(set.config.sys.cbmcTrials.bonded.dih),
  nLJTrialsFirst(set.config.sys.cbmcTrials.nonbonded.first),
  nLJTrialsNth(set.config.sys.cbmcTrials.nonbonded.nth), 
  positions(*multiPositions)
{
   uint maxLJTrials = nLJTrialsFirst;
   if ( nLJTrialsNth > nLJTrialsFirst )
     maxLJTrials = nLJTrialsNth;
   for(uint i = 0; i < MAX_BONDS; ++i)
   {
      multiPositions[i] = XYZArray(maxLJTrials);
   }
   inter = new double[maxLJTrials];
   bonded = new double[maxLJTrials];
   nonbonded = new double[maxLJTrials];
   nonbonded_1_4 = new double[maxLJTrials];
   ljWeights = new double[maxLJTrials];
   self = new double[maxLJTrials];
   correction = new double[maxLJTrials];
   real = new double[maxLJTrials];

   uint trialMax = std::max(nAngleTrials, nDihTrials);
   angleEnergy = new double[trialMax];
   angleWeights = new double[trialMax];
   angles = new double[trialMax];
}

inline DCData::~DCData()
{
   delete[] inter;
   delete[] bonded;
   delete[] nonbonded;
   delete[] nonbonded_1_4;
   delete[] self;
   delete[] correction;
   delete[] real;
   delete[] ljWeights;
   delete[] angles;
   delete[] angleWeights;
   delete[] angleEnergy;
}

}

#endif

