/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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


namespace fail_state
{
const uint NO_FAIL = 1;
const uint ROTATE_ON_SINGLE_ATOM = 2;
const uint NO_MOL_OF_KIND_IN_BOX = 3;
}

}

#endif /*MOVES_CONST_H*/