/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef MOVES_CONST_H
#define MOVES_CONST_H

#include <string>
#include <vector>

#include "BasicTypes.h"
#include "EnsemblePreprocessor.h"

namespace mv {
#if ENSEMBLE == GEMC || ENSEMBLE == NPT
//////////////////////////////////////////////////////////
// ENSEMBLE KINDS
const uint GEMC_NVT = 0;
const uint GEMC_NPT = 1;
#endif

// GENERAL MOVE CONSTANTS

//////////////////////////////////////////////////////////
const uint DISPLACE = 0;
const uint ROTATE = 1;
const uint MULTIPARTICLE = 2;
const uint MULTIPARTICLE_BM = 3;
<<<<<<< HEAD
#if ENSEMBLE == NVT
const uint INTRA_SWAP = 4;
const uint REGROWTH = 5;
const uint INTRA_MEMC = 6;
const uint CRANKSHAFT = 7;
const uint MOVE_KINDS_TOTAL = 8;
#elif ENSEMBLE == GCMC
const uint INTRA_SWAP = 4;
const uint REGROWTH = 5;
const uint INTRA_MEMC = 6;
const uint CRANKSHAFT = 7;
const uint MEMC = 8;
const uint MOL_TRANSFER = 9;
const uint CFCMC = 10;
const uint MOVE_KINDS_TOTAL = 11;
#elif ENSEMBLE == GEMC
const uint VOL_TRANSFER = 4;
const uint INTRA_SWAP = 5;
const uint REGROWTH = 6;
const uint INTRA_MEMC = 7;
const uint CRANKSHAFT = 8;
const uint MEMC = 9;
const uint MOL_TRANSFER = 10;
const uint CFCMC = 11;
const uint MOVE_KINDS_TOTAL = 12;
#elif ENSEMBLE == NPT
const uint VOL_TRANSFER = 4;
const uint INTRA_SWAP = 5;
const uint REGROWTH = 6;
const uint INTRA_MEMC = 7;
const uint CRANKSHAFT = 8;
const uint MOVE_KINDS_TOTAL = 9;
=======
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
>>>>>>> 0d5a1882dac7ee07e3118a3faf3dd3cfe681cb16
#endif

const uint BOX0 = 0;
const uint BOX1 = 1;

//////////////////////////////////////////////////////////

<<<<<<< HEAD
//NVT : 1. Disp (box 0)           2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. MultiParticle_BM(box0) 5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)      8. CrankShaft (box 0)
//
//GCMC: 1. Disp (box 0)         2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. MultiParticle_BM(box0)
//      5. IntraSwap (box 0)    6. Regrowth (box 0)   7. IntraMEMC (box 0)
//      8. CrankShaft (box 0)   9. MEMC (box 0)       10. Deletion (box 0)
//      11. Insertion (box 0)   12. CFCMC (box 0)     13. CFCMC (box 1)
//
//GEMC: 1. Disp (box 0)           2. Disp (box 1)
//      3. MultiParticle (box 0)  4. MultiParticle (box 1)
//      5. MultiParticle_BM(box0) 6. MultiParticle_BM (box1)
//      7. Rotate (box 0)        8. Rotate (box 1)
//      9. Vol. (b0->b1)        10. Vol. (b1->b0)
//     11. IntraSwap (box 0)    12. IntraSwap (box 1)
//     13. Regrowth (box 0)     14. Regrowth (box 1)
//     15. IntraMEMC (box 0)    16. IntraMEMC (box 1)
//     17. CrankShaft (box 0)   18. CrankShaft (box 1)
//     19. MEMC (box 0)         20. MEMC (box 1)
//     21. Mol Trans (b0->b1)   22. Mol Trans (b1->b0)
//     23. CFCMC (b0->b1)       24. CFCMC (b1->b0)
//
//NPT : 1. Disp (box 0)         2. Rotate (box 0)     3. MultiParticle (box 0)
//      4. MultiParticle_BM (box0)
//      5. Vol. (box 0)         6. IntraSwap (box 0)  7. Regrowth (box 0)
//      8. IntraMEMC (box 0)    9. CrankShaft (box 0)
//AUTO REJECTION OR ACCEPTANCE FLAGS
=======
// NVT : 1. Disp (box 0)             2. Rotate (box 0)     3. MultiParticle (box
// 0)
//      4. MultiParticle_BM (box 0) 5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)        8. CrankShaft (box 0) 9. IntraTargetedSwap
//      (box 0)
//
// GCMC: 1. Disp (box 0)             2. Rotate (box 0)     3. MultiParticle (box
// 0)
//      4. MultiParticle_BM (box 0) 5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)        8. CrankShaft (box 0) 9. IntraTargetedSwap
//      (box 0)
//     10. MEMC (box 0)            11. Deletion (box 0)  12. Insertion (box 0)
//     13. NE_MTMC (box 0)         14. NE_MTMC (box 1)   15. TargetedSwap (box
//     0)
//     16. TargetedSwap (box 1)
//
// GEMC: 1. Disp (box 0)               2. Disp (box 1)
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
// NPT : 1. Disp (box 0)             2. Rotate (box 0)     3. MultiParticle (box
// 0)
//      4. MultiParticle_BM (box 0) 5. IntraSwap (box 0)  6. Regrowth (box 0)
//      7. IntraMEMC (box 0)        8. CrankShaft (box 0) 9. IntraTargetedSwap
//      (box 0)
//     10. Vol. (box 0)
// AUTO REJECTION OR ACCEPTANCE FLAGS
>>>>>>> 0d5a1882dac7ee07e3118a3faf3dd3cfe681cb16

namespace fail_state {
const uint NO_FAIL = 1;
const uint ROTATE_ON_SINGLE_ATOM = 2;
const uint NO_MOL_OF_KIND_IN_BOX = 3;
} // namespace fail_state

} // end namespace mv
#endif /*MOVES_CONST_H*/
