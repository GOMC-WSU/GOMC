/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOVES_CONST_H
#define MOVES_CONST_H

#include "EnsemblePreprocessor.h"
#include "BasicTypes.h"

#include <vector>
#include <string>


namespace mv
{
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
//////////////////////////////////////////////////////////
// ENSEMBLE KINDS
const uint GEMC_NVT = 0;
const uint GEMC_NPT = 1;
#endif

//GENERAL MOVE CONSTANTS

//////////////////////////////////////////////////////////
const uint DISPLACE = 0;
const uint ROTATE = 1;
const uint MULTIPARTICLE = 2;
#if ENSEMBLE == NVT
const uint INTRA_SWAP = 3;
const uint REGROWTH = 4;
const uint INTRA_MEMC = 5;
const uint CRANKSHAFT = 6;
const uint MOVE_KINDS_TOTAL = 7;
#elif ENSEMBLE == GCMC
const uint INTRA_SWAP = 3;
const uint REGROWTH = 4;
const uint INTRA_MEMC = 5;
const uint CRANKSHAFT = 6;
const uint MEMC = 7;
const uint MOL_TRANSFER = 8;
const uint MOVE_KINDS_TOTAL = 9;
#elif ENSEMBLE == GEMC
const uint VOL_TRANSFER = 3;
const uint INTRA_SWAP = 4;
const uint REGROWTH = 5;
const uint INTRA_MEMC = 6;
const uint CRANKSHAFT = 7;
const uint MEMC = 8;
const uint MOL_TRANSFER = 9;
const uint MOVE_KINDS_TOTAL = 10;
#elif ENSEMBLE == NPT
const uint VOL_TRANSFER = 3;
const uint INTRA_SWAP = 4;
const uint REGROWTH = 5;
const uint INTRA_MEMC = 6;
const uint CRANKSHAFT = 7;
const uint MOVE_KINDS_TOTAL = 8;
#endif

const uint BOX0 = 0;
const uint BOX1 = 1;

const uint PICK_ONE_MOLECULE = 0;
const uint PICK_ALL_MOLECULE = 1;

const uint TRANS_DISPLACE = 0;
const uint TRANS_ROTATE = 1;
const uint TRANS_SCALE_COM = 2;
//Vlugt CBMC scheme, simple gen.
//
//Weight comes from the intramolecular potential (LJ(inter, intra)+TC)
//While intra only is used by stochastic generation and pick of final
//trial pos to use.
//
//This is only valid, generally for linear molecules.
const uint TRANS_GEN_TRIAL_LINEAR = 3;
//Pick angles that match Boltzmann acc. for n_ch_bend
//Pass them to dihedrals, pick one via Boltzman.
//
//Disadvantage:
// n_ch_lj is limited to a single generated conformation.
// LJ particle is then locked into place as we must always accept
// pick it when generating trial conf. in the old box.
// hence if we had more trial conformations over a free new box particle
// it'd introduce a bias.
//
// i.e. We must pick ONE and only one of each "thing" per conformation
// phase, hence in the last phase we always have to pick the old conf.
// in the old box, hence you can't do biased sel. of locations based
// on the LJ pot.
//
//IN OLD BOX:
// if first angle, usual actual angle, dihedral, lj in old.(i=j=k=1)
// if first dihedral, dihedral, lj is old (j=k=1)
// if first lj ... just LJ is old., i.e. in old box
// Now each phase incorporates the weight of the last, so the pick
const uint GEN_TRIAL_DECOUPLED = 4;
//Pick angles that match boltz. acc., then pick dihedrals based on
// wi * Wb / Wt
//We now have if i = j = k = 1 is the only case where the old conf.
//is used.  Otherwise since each step depends on the last
const uint GEN_TRIAL_COUPLED = 5;
const uint GEN_TRIAL_COUPLED_DECOUPLED = 6;
const uint TRANS_KINDS_TOTAL = 7;

const uint IT_SINGLE_MOL = 0;
const uint IT_ALL_MOL = 1;
const uint IT_KINDS_TOTAL = 2;

//////////////////////////////////////////////////////////

//NVT : 1. Disp (box 0)         2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. IntraSwap (box 0)    5. Regrowth (box 0)   6. IntraMEMC (box 0) 
//      7. CrankShaft (box 0)
//
//GCMC: 1. Disp (box 0)         2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. IntraSwap (box 0)    5. Regrowth (box 0)   6. IntraMEMC (box 0)
//      7. CrankShaft (box 0)   8. MEMC (box 0)       9. Deletion (box 0)
//      10. Insertion (box 0)
//
//GEMC: 1. Disp (box 0)          2. Disp (box 1)
//      3. MultiParticle (box 0) 4. MultiParticle (box 1)
//      5. Rotate (box 0)        6. Rotate (box 1)
//      7. Vol. (b0->b1)         8. Vol. (b1->b0)
//      9. IntraSwap (box 0)    10. IntraSwap (box 1)
//     11. Regrowth (box 0)     12. Regrowth (box 1)
//     13. IntraMEMC (box 0)    14. IntraMEMC (box 1)
//     15. CrankShaft (box 0)   16. CrankShaft (box 1)
//     17. MEMC (box 0)         18. MEMC (box 1)
//     19. Mol Trans (b0->b1),  20. Mol Trans (b1->b0)
//
//NPT : 1. Disp (box 0)         2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. Vol. (box 0)         5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)    8. CrankShaft (box 0)
//AUTO REJECTION OR ACCEPTANCE FLAGS

//early exit flags.
namespace auto_accept
{
const uint ONLY_IN_BOX_ROT_OR_DISP = 0;
}

namespace fail_state
{
const uint NO_FAIL = 1;
const uint ROTATE_ON_SINGLE_ATOM = 2;
const uint NO_MOL_OF_KIND_IN_BOX = 3;
const uint INNER_CUTOFF_NEW_TRIAL_POS = 4;
const uint VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF = 5;
}

}

#endif /*MOVES_CONST_H*/