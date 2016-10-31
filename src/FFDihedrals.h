/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.8
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_DIHEDRALS_H
#define FF_DIHEDRALS_H

#include <math.h>       //cos, pow

#include "FFConst.h"    //GetRot
#include "SubdividedArray.h" //Subdivisions in master array for dih sets.

#include "../lib/NumLib.h"    //Sq
#include "../lib/BasicTypes.h" //For "uint"

namespace ff_setup { class Dihedral; }

//FFDihKind
//Stores parameters for dihedral kinds
//dihedrals may require several sets of parameters to calculate,
//these are stored in consecutive array elements, each kind's
//parameter sets are delinated by subdiv
//Dihedrals may be calculated by several different styles
class FFDihedrals
{
  public:
   //calculate the energy of dih kind at angle phi
   double Calc(const uint kind, const double phi) const;
   //Initialize with data from parameter files
   void Init(ff_setup::Dihedral const& dih);

   FFDihedrals(void) : Kchi(NULL), delta(NULL), n(NULL) {}
   ~FFDihedrals(void);

 private:

   //dih kind params
   SubdividedArray subdiv;
   double * Kchi, * delta;
   uint *n;
};

inline double FFDihedrals::Calc(const uint kind, const double phi) const
{
   double sum = 0.0;
   for(uint i = subdiv.Begin(kind); i != subdiv.End(kind); ++i)
   {
      sum += Kchi[i] * (1 + cos(n[i] * phi - delta[i]));
   }
   return sum;
}

#endif /*FF_DIHEDRALS_H*/
