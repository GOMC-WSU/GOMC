#ifndef MOVES_CONST_H
#define MOVES_CONST_H

#include "EnsemblePreprocessor.h"
#include "../lib/BasicTypes.h"

#include <vector>
#include <string>


namespace mv
{
#if ENSEMBLE == GEMC
   //////////////////////////////////////////////////////////
   // ENSEMBLE KINDS
   const uint GEMC_NVT = 0;
   const uint GEMC_NPT = 1;
#endif

   //GENERAL MOVE CONSTANTS

   //////////////////////////////////////////////////////////
   const uint DISPLACE = 0;
   const uint ROTATE = 1;
#if ENSEMBLE == NVT
   const uint INTRA_SWAP = 2;
   const uint MOVE_KINDS_TOTAL = 3;
#elif ENSEMBLE == GCMC
   const uint MOL_TRANSFER=2;
   const uint INTRA_SWAP = 3;
   const uint MOVE_KINDS_TOTAL=4;
#elif ENSEMBLE == GEMC
   const uint VOL_TRANSFER=2;
   const uint MOL_TRANSFER=3;
   const uint INTRA_SWAP = 4;
   const uint MOVE_KINDS_TOTAL=5;
#endif

   const uint BOX0 = 0;
   const uint BOX1 = 1;

   const uint PICK_ONE_MOLECULE=0;
   const uint PICK_ALL_MOLECULE=1;

   const uint TRANS_DISPLACE=0;
   const uint TRANS_ROTATE=1;
   const uint TRANS_SCALE_COM=2;
   //Vlugt CBMC scheme, simple gen.
   //
   //Weight comes from the intramolecular potential (LJ(inter, intra)+TC)
   //While intra only is used by stochastic generation and pick of final
   //trial pos to use.
   //
   //This is only valid, generally for linear molecules.
   const uint TRANS_GEN_TRIAL_LINEAR=3; 
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
   const uint GEN_TRIAL_DECOUPLED=4;
   //Pick angles that match boltz. acc., then pick dihedrals based on
   // wi * Wb / Wt
   //We now have if i = j = k = 1 is the only case where the old conf.
   //is used.  Otherwise since each step depends on the last
   const uint GEN_TRIAL_COUPLED=5;
   const uint GEN_TRIAL_COUPLED_DECOUPLED=6;
   const uint TRANS_KINDS_TOTAL=7;

   const uint IT_SINGLE_MOL=0;
   const uint IT_ALL_MOL=1;
   const uint IT_KINDS_TOTAL=2;

   //////////////////////////////////////////////////////////
   
   //NVT : 1. Disp (box 0) 2. Rotate (box 0) 3. IntraSwap (box 0)
   //GCMC: 1. Disp (box 0) 2. Rotate (box 0) 3. Deletion (box 0)
   //      4. Insertion (box 0) 5. IntraSwap (box 0)
   //GEMC: 1. Disp (box 0) 2. Disp (box 1) 3. Rotate (box 0) 4. Rotate (box 1)
   //      5. Vol. (b0->b1) 6. Vol. (b1->b0) 7. Mol Trans (b0->b1), lin.
   //      8. Mol Trans (b1->b0), lin. 9. IntraSwap (box 0)
   //     10. IntraSwap (box 1)

#if ENSEMBLE == NVT
   const uint COUNT = 3;
   const uint SCALEABLE = 2;
#elif ENSEMBLE == GCMC
   const uint COUNT = 5;
   const uint SCALEABLE = 2;
#elif ENSEMBLE == GEMC
   const uint COUNT = 10;
   const uint SCALEABLE = 6;
#endif


   //AUTO REJECTION OR ACCEPTANCE FLAGS

   //early exit flags.
   namespace auto_accept 
   { 
      const uint ONLY_IN_BOX_ROT_OR_DISP=0; 
   }

   namespace fail_state
   { 
      const uint NO_FAIL = 1;
      const uint ROTATE_ON_SINGLE_ATOM = 2;
      const uint NO_MOL_OF_KIND_IN_BOX = 3;
      const uint INNER_CUTOFF_NEW_TRIAL_POS=4;
      const uint VOL_TRANS_WOULD_SHRINK_BOX_BELOW_CUTOFF=5;
   }

   inline void GetMoveMajIndex(uint & maj, uint & subDiv, const uint sub)
   { 
#if ENSEMBLE == NVT || ENSEMBLE == GEMC
      maj = sub/BOX_TOTAL;
      subDiv = BOX_TOTAL;
#elif ENSEMBLE == GCMC
      if (maj > MOL_TRANSFER)
	 maj--;
      //Adjust for single box moves.
      maj = sub;
      subDiv = BOX_TOTAL;
      if (maj < MOL_TRANSFER)
	 subDiv = 1;
#endif
   }

   inline uint GetMoveSubIndex(const uint maj, const uint b = 0)
   { 
#if ENSEMBLE == GEMC || ENSEMBLE == NVT
      return maj*BOX_TOTAL + b; 
#else

       if(maj == mv::INTRA_SWAP)
	 return maj + b + 1; //in GCMC we care just about box 0. [3+0+1]
       else
	 return maj+b; //in GCMC we have deletion and insertion,
                       //we care about two box. [2+0 or 2+1]
#endif
   }

   //Names of above moves as strings for output.
   std::vector<std::string> MoveNames();
   const std::vector<std::string> MOVE_NAME(MoveNames());
   std::vector<std::string> ScaleMoveNames();
   const std::vector<std::string> SCALE_MOVE_NAME(ScaleMoveNames());
   //Used enums -- immutable and take no space

}

#endif /*MOVES_CONST_H*/
