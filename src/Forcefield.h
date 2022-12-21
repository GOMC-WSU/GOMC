/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FORCEFIELD_H
#define FORCEFIELD_H

// Member classes
#include "FFAngles.h"
#include "FFBonds.h"
#include "FFDihedrals.h"
#include "FFParticle.h"

namespace config_setup {
// class FFValues;
struct FFKind;
// class Temperature;
struct SystemVals;
} // namespace config_setup
class FFSetup;
class Setup;
class FFPrintout;
struct FFParticle;

typedef double BOX_SIZE_DOUBLE_ARRAY[BOX_TOTAL];

class Forcefield {
public:
  friend class FFPrintout;

  Forcefield();
  ~Forcefield();
  // Initialize contained FFxxxx structs from setup data
  void Init(const Setup &set);
  uint GetWolfKind(void);
  uint GetCoulKind(void); 
  void SetWolfKind(uint);
  void SetCoulKind(uint); 
  void SetWolfAlpha(double, uint);
  BOX_SIZE_DOUBLE_ARRAY& GetWolfAlpha(void);
  BOX_SIZE_DOUBLE_ARRAY& GetWolfFactor1(void);
  BOX_SIZE_DOUBLE_ARRAY& GetWolfFactor2(void);
  BOX_SIZE_DOUBLE_ARRAY& GetWolfFactor3(void);
  BOX_SIZE_DOUBLE_ARRAY& GetRCutCoulomb(void);
  BOX_SIZE_DOUBLE_ARRAY& GetRCutCoulombSq(void);

  FFParticle *particles; //!< For LJ/Mie energy between unbonded atoms
  // for LJ, shift and switch type
  FFBonds bonds;             //!< For bond stretching energy
  FFAngles *angles;          //!< For 3-atom bending energy
  FFDihedrals dihedrals;     //!< For 4-atom torsional rotation energy
  bool useLRC;               //!< Use long-range tail corrections if true
  bool useIPC;               //!< Use impulse pressure corrections if true
  double T_in_K;             //!< System temp in Kelvin
  double beta;               //!< Thermodynamic beta = 1/(T) K^-1)
  double rCut, rCutSq;       //!< Cutoff radius for LJ/Mie potential (angstroms)
  double rCutLow, rCutLowSq; //!< Cutoff min for Electrostatic (angstroms)
  double rCutCoulomb[BOX_TOTAL];   //!< Cutoff Coulomb interaction(angstroms)
  double rCutCoulombSq[BOX_TOTAL]; //!< Cutoff Coulomb interaction(angstroms)
  double alpha[BOX_TOTAL];         // Ewald sum terms
  double alphaSq[BOX_TOTAL];       // Ewald sum terms
  double recip_rcut[BOX_TOTAL];    // Ewald sum terms
  double recip_rcut_Sq[BOX_TOTAL]; // Ewald sum terms
  double tolerance;                // Ewald sum terms
  double rswitch;                  // Switch distance
  double dielectric;               // dielectric for martini
  double scaling_14; //!< Scaling factor for 1-4 pairs' ewald interactions
  double sc_alpha;   // Free energy parameter
  double sc_sigma, sc_sigma_6; // Free energy parameter

  bool OneThree, OneFour, OneN; // To include 1-3, 1-4 and more interaction
  bool electrostatic, ewald;    // To consider columb interaction
  bool vdwGeometricSigma;       // For sigma combining rule
  bool isMartini;
  bool exp6;
  bool freeEnergy, sc_coul; // Free energy parameter
  uint vdwKind;             // To define VdW type, standard, shift or switch
  uint exckind;             // To define  exclude kind, 1-2, 1-3, 1-4
  uint sc_power;            // Free energy parameter
#if ENSEMBLE == GCMC
  bool isFugacity; // To check if we are using fugacity instead of chemical
                   // potential
#endif
  // Initialize primitive member variables from setup data
  bool wolf;
  double wolfAlpha[BOX_TOTAL]; //alpha term for Wolf Electrostatic and constant factors
  double wolfFactor1[BOX_TOTAL]; //alpha term for Wolf Electrostatic and constant factors
  double wolfFactor2[BOX_TOTAL];  //alpha term for Wolf Electrostatic and constant factors
  double wolfFactor3[BOX_TOTAL]; //alpha term for Wolf Electrostatic and constant factors
  uint coulKind, wolfKind;        //To define Coul type (if Wolf), dampened shift potential or dampened shift force


private:

  void InitBasicVals(config_setup::SystemVals const &val,
                     config_setup::FFKind const &ffKind);
};

#endif /*FORCEFIELD_H*/
