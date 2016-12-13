/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.9
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FORCEFIELD_H
#define FORCEFIELD_H

//Member classes
#include "FFParticle.h"
#include "FFShift.h"
#include "FFSwitch.h"
#include "FFSwitchMartini.h"
#include "FFBonds.h"
#include "FFAngles.h"
#include "FFDihedrals.h"
#include "PRNG.h"

namespace config_setup
{
//class FFValues;
struct FFKind;
//class Temperature;
struct SystemVals;
}
class FFSetup;
class Setup;
class FFPrintout;

class Forcefield
{
public:
  friend class FFPrintout;

  Forcefield();
  ~Forcefield();
  //Initialize contained FFxxxx structs from setup data
  void Init(const Setup& set);


  FFParticle * particles;    //!<For LJ/Mie energy between unbonded atoms
  // for LJ, shift and switch type
  FFBonds bonds;             //!<For bond stretching energy
  FFAngles * angles;         //!<For 3-atom bending energy
  FFDihedrals dihedrals;     //!<For 4-atom torsional rotation energy
  bool useLRC;               //!<Use long-range tail corrections if true
  double T_in_K;             //!<System temp in Kelvin
  double rCut;               //!<Cutoff radius for LJ/Mie potential (angstroms)
  double rCutSq;             //!<Cutoff radius for LJ/Mie potential squared (a^2)
  double rCutOver2;          //!<Cutoff radius for LJ/Mie potential over 2 (a)

  double rOn;                // for switch tup of LJ (angstroms)
  double ronSq;              // for switch tup of LJ (a^2)

  //XXX 1-4 pairs are not yet implemented
  double scl_14;           //!<Scaling factor for 1-4 pairs' ewald interactions
  double beta;             //!<Thermodynamic beta = 1/(T) K^-1)

  bool OneThree, OneFour, OneN;  // to include 1-3, 1-4 and more interaction
  bool electrostatic, ewald;     //to consider columb interaction
  double alpha, recip_rcut;      // ewald terms
  uint vdwKind;      // to define VdW type, standard, shift or switch
  uint exckind;      // to define  exclude kind, 1-2, 1-3, 1-4

private:
  //Initialize primitive member variables from setup data

  void InitBasicVals(config_setup::SystemVals const& val,
                     config_setup::FFKind const& ffKind);

};

#endif /*FORCEFIELD_H*/
