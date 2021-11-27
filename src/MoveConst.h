/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
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
const uint MULTIPARTICLE_BM = 3;
const uint INTRA_SWAP = 4;
const uint REGROWTH = 5;
const uint INTRA_MEMC = 6;
const uint CRANKSHAFT = 7;
const uint INTRA_TARGETED_SWAP = 8;
#if ENSEMBLE == NVT
const uint MOVE_KINDS_TOTAL = 9;
#elif ENSEMBLE == GCMC
const uint MEMC = 9;
const uint MOL_TRANSFER = 10;
const uint NE_MTMC = 11;
const uint TARGETED_SWAP = 12;
const uint MOVE_KINDS_TOTAL = 13;
#elif ENSEMBLE == GEMC
const uint MEMC = 9;
const uint MOL_TRANSFER = 10;
const uint NE_MTMC = 11;
const uint TARGETED_SWAP = 12;
const uint VOL_TRANSFER = 13;
const uint MOVE_KINDS_TOTAL = 14;
#elif ENSEMBLE == NPT
const uint VOL_TRANSFER = 9;
const uint MOVE_KINDS_TOTAL = 10;
#endif

const uint BOX0 = 0;
const uint BOX1 = 1;

//////////////////////////////////////////////////////////

//NVT : 1. Disp (box 0)             2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. MultiParticle_BM (box 0) 5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)        8. CrankShaft (box 0) 9. IntraTargetedSwap (box 0)
//
//GCMC: 1. Disp (box 0)             2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. MultiParticle_BM (box 0) 5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)        8. CrankShaft (box 0) 9. IntraTargetedSwap (box 0) 
//     10. MEMC (box 0)            11. Deletion (box 0)  12. Insertion (box 0)
//     13. NE_MTMC (box 0)         14. NE_MTMC (box 1)   15. TargetedSwap (box 0)
//     16. TargetedSwap (box 1)
//
//GEMC: 1. Disp (box 0)               2. Disp (box 1)
//      3. MultiParticle (box 0)      4. MultiParticle (box 1)
//      5. MultiParticle_BM(box 0)    6. MultiParticle_BM(box 1)
//      7. Rotate (box 0)             8. Rotate (box 1)
//      9. IntraSwap (box 0)         10. IntraSwap (box 1)
//     11. Regrowth (box 0)          12. Regrowth (box 1)
//     13. IntraMEMC (box 0)         14. IntraMEMC (box 1)
//     15. CrankShaft (box 0)        16. CrankShaft (box 1)
//     17. IntraTargetedSwap (box 0) 18. IntraTargetedSwap (box 1)
//     19. MEMC (b0->b1)             20. MEMC (b1->b0)
//     21. Mol Trans (b0->b1)        22. Mol Trans (b1->b0)
//     23. NE_MTMC (b0->b1)          24. NE_MTMC (b1->b0)
//     25. TargetedSwap (b0->b1)     26. TargetedSwap (b1->b0)
//     27. Vol. (b0->b1)             28. Vol. (b1->b0)
//
//NPT : 1. Disp (box 0)             2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. MultiParticle_BM (box 0) 5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)        8. CrankShaft (box 0) 9. IntraTargetedSwap (box 0)
//     10. Vol. (box 0)         
//AUTO REJECTION OR ACCEPTANCE FLAGS


namespace fail_state
{
const uint NO_FAIL = 1;
const uint ROTATE_ON_SINGLE_ATOM = 2;
const uint NO_MOL_OF_KIND_IN_BOX = 3;
}


#ifndef NDEBUG
std::string printMoveName(uint moveConst)
{
  std::string moveName;
  switch (moveConst) {
    case 0:
      moveName = "DISPLACE";
      break;
    case 1:
      moveName = "ROTATE";
      break;
    case 2:
      moveName = "MULTIPARTICLE";
      break;
    case 3:
      moveName = "BROWNIAN MULTIPARTICLE";
      break;
    case 4:
      moveName = "INTRA_SWAP";
      break;
    case 5:
      moveName = "REGROWTH";
      break;
    case 6:
      moveName = "INTRA_MEMC";
      break;
    case 7:
      moveName = "CRANKSHAFT";
      break;
    case 8:
      moveName = "INTRA_TARGETED_SWAP";
      break;
#if ENSEMBLE == NPT
    case 9:
      moveName = "VOL_TRANSFER";
      break;
#elif ENSEMBLE == GCMC || ENSEMBLE == GEMC
    case 9:
      moveName = "MEMC";
      break;
    case 10:
      moveName = "MOL_TRANSFER";
      break;
    case 11:
      moveName = "NE_MTMC";
      break;
    case 12:
      moveName = "TARGETED_SWAP";
      break;
#elif ENSEMBLE == GEMC
    case 13:
      moveName = "VOL_TRANSFER";
      break;
#endif
    default:
      moveName = "Update printMove() function!!! Undefined";
  }

  return moveName;
}
#endif /*NDEBUG*/
} //end namespace mv
#endif /*MOVES_CONST_H*/
