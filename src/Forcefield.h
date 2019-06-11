/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
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
class FFParticle;

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
  FFBonds bonds;                  //!<For bond stretching energy
  FFAngles * angles;              //!<For 3-atom bending energy
  FFDihedrals dihedrals;          //!<For 4-atom torsional rotation energy
  bool useLRC;                    //!<Use long-range tail corrections if true
  real T_in_K;                  //!<System temp in Kelvin
  real beta;                    //!<Thermodynamic beta = 1/(T) K^-1)
  real rCut, rCutSq;            //!<Cutoff radius for LJ/Mie potential (angstroms)
  real rCutLow, rCutLowSq;      //!<Cutoff min for Electrostatic (angstroms)
  real rCutCoulomb[BOX_TOTAL];  //!<Cutoff Coulomb interaction(angstroms)
  real rCutCoulombSq[BOX_TOTAL]; //!<Cutoff Coulomb interaction(angstroms)
  real alpha[BOX_TOTAL];        //Ewald sum terms
  real alphaSq[BOX_TOTAL];      //Ewald sum terms
  real recip_rcut[BOX_TOTAL];   //Ewald sum terms
  real recip_rcut_Sq[BOX_TOTAL]; //Ewald sum terms
  real tolerance;               //Ewald sum terms
  real rswitch;                 //Switch distance
  real dielectric;              //dielectric for martini
  real scaling_14;              //!<Scaling factor for 1-4 pairs' ewald interactions

  bool OneThree, OneFour, OneN;   //To include 1-3, 1-4 and more interaction
  bool electrostatic, ewald;      //To consider columb interaction
  bool vdwGeometricSigma;         //For sigma combining rule
  bool isMartini;
  uint vdwKind;                   //To define VdW type, standard, shift or switch
  uint exckind;                   //To define  exclude kind, 1-2, 1-3, 1-4
#if ENSEMBLE == GCMC
  bool isFugacity;                //To check if we are using fugacity instead of chemical potential
#endif

private:
  //Initialize primitive member variables from setup data

  void InitBasicVals(config_setup::SystemVals const& val,
                     config_setup::FFKind const& ffKind);

};

#endif /*FORCEFIELD_H*/
