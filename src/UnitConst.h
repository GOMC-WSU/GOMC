/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef UNIT_CONST_H
#define UNIT_CONST_H

//Constants for unit conversion
namespace unit
{

// Molecules: molecules/A3 --> mol/cm3
//
// 1cm3=1e24 A3; 6.0221413e23 molecules/mol
//
// const 1e24/6.02...e23 A3/molecules * mol/cm3
//
static const real MOLECULES_PER_A3_TO_MOL_PER_CM3 = 1.660539277;

//
// K * molecule / A3 --> bar
//
// P(K*molecule/A3) * 1e30 A3/1 m3 * 1 mol/6.02...e23 molecule *
//   1.987204 cal/mol/K * 4.184 J /1 cal * ...
//   1 J/m3 (J/m3 = N*m/m^3 = N/m2 = Pa) ...
//   1 bar / 1e5 Pa
//   = 1e30*1.987204*4.184/6.0221413e23/1e5
//
static const real K_MOLECULE_PER_A3_TO_BAR = 138.0649;

// convert K to mN/m (dyne)
static const real K_TO_MN_PER_M = 1.380649;

//
// bar --> A3 / (K * molecule)
//
// P(bar) * 1e5 Pa / 1 bar *  1 J/m3 / 1 Pa
//   1 cal / 4.184 J * 1 mol*K / 1.987204 K *
//   6.02...e23 molecule / 1 mol * 1m3 / 1e30 A3 -->
//   A3 / (K * molecule)
//
static const real BAR_TO_K_MOLECULE_PER_A3 = 0.007242971;

//
// K --> kJ/mol
//
// delta U or P*(delta V) (K) * 0.0019... kcal/K/mol * 1000 cal/1 kcal *
//     4.184 J/1 cal * 1 kJ/1000 J
//  = 0.001987204*4.184
//
static const real K_TO_KJ_PER_MOL = 0.008314462;

//
// dimensionality of the virial (3D = 3)
//
static const uint DIMENSIONALITY = 3;

}

#endif /*UNIT_CONST_H*/
