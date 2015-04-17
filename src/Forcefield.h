/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FORCEFIELD_H
#define FORCEFIELD_H

//Member classes
#include "FFParticle.h"
#include "FFBonds.h"
#include "FFAngles.h"
#include "FFDihedrals.h"
#include "PRNG.h"

namespace config_setup { class FFValues; }
class FFSetup;
class Setup;

class Forcefield
{
 public:   
      
   //Initialize contained FFxxxx structs from setup data
   void Init(const Setup& set);


   FFParticle particles;      //!<For LJ/Mie energy between unbonded atoms
   FFBonds bonds;             //!<For bond stretching energy
   FFAngles angles;           //!<For 3-atom bending energy
   FFDihedrals dihedrals;     //!<For 4-atom torsional rotation energy
   bool useLRC;               //!<Use long-range tail corrections if true
   double T_in_K;             //!<System temp in Kelvin
   double rCut;               //!<Cutoff radius for LJ/Mie potential (angstroms)
   double rCutSq;             //!<Cutoff radius for LJ/Mie potential squared (a^2)
   double rCutOver2;          //!<Cutoff radius for LJ/Mie potential over 2 (a)
   //XXX 1-4 pairs are not yet implemented
   double scl_14;             //!<Scaling factor for 1-4 pairs' ewald interactions
   double beta;               //!<Thermodynamic beta = 1/(T) K^-1)


 private:
   //Initialize primitive member variables from setup data
   void InitBasicVals(config_setup::FFValues const& val);

};

#endif /*FORCEFIELD_H*/

