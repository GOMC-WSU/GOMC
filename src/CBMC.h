#ifndef CBMC_H
#define CBMC_H

/*    CBMC.h
*     Base Class for CBMC algorithms
*     Also includes Factory Functions for same
*
*/

#include "../lib/BasicTypes.h"

class MolPick;
class Forcefield;
class MoleculeKind;
class Setup;
class System;

namespace cbmc
{
   class TrialMol;

   class CBMC
   {
   public:
      //Builds a new molecule using a CBMC algorithm, oldMol and newMol
      //will be modified to contain the energies of the old and new sites
      virtual void Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex) = 0;

      virtual ~CBMC() {}
   };

   //Max allowed bonds to any atom
   static const uint MAX_BONDS = 6;
   //Factory function, determines, prepares and returns appropriate CBMC
   CBMC* MakeCBMC(System& sys, const Forcefield& ff,
      const MoleculeKind& kind, const Setup& set);
}


#endif
