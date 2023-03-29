/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef CBMC_H
#define CBMC_H

/*    CBMC.h
 *     Base Class for CBMC algorithms
 *     Also includes Factory Functions for same
 *
 */

#include "BasicTypes.h"
#include "Forcefield.h"

struct MolPick;
class Forcefield;
class MoleculeKind;
class Setup;
class System;

namespace cbmc {
class TrialMol;

class CBMC {
public:
  // Builds a new molecule using a CBMC algorithm, oldMol and newMol
  // will be modified to contain the energies of the old and new sites
  virtual void Build(TrialMol &oldMol, TrialMol &newMol, uint molIndex) = 0;

  // Regrowing the molecule using a CBMC algorithm, oldMol and newMol
  virtual void Regrowth(TrialMol &oldMol, TrialMol &newMol, uint molIndex) = 0;

  // Rotate the atoms between two nodes arounde  the vector that connects two
  // nodes using crank shaft algorithm, oldMol and newMol
  virtual void CrankShaft(TrialMol &oldMol, TrialMol &newMol,
                          uint molIndex) = 0;

  // Rigid insertion of molecule and perform position and rotational trial
  virtual void BuildIDNew(TrialMol &newMol, uint molIndex) = 0;
  virtual void BuildIDOld(TrialMol &oldMol, uint molIndex) = 0;

  // Build the molecule using CD-CBMC
  virtual void BuildNew(TrialMol &newMol, uint molIndex) = 0;
  virtual void BuildOld(TrialMol &oldMol, uint molIndex) = 0;

  // Grow the molecule from predefined atom (node)
  virtual void BuildGrowNew(TrialMol &newMol, uint molIndex) = 0;
  virtual void BuildGrowOld(TrialMol &oldMol, uint molIndex) = 0;

  // grow the molecule with starting seed within subVolume
  virtual void BuildGrowInCav(TrialMol &oldMol, TrialMol &newMol,
                              uint molIndex) = 0;

  virtual ~CBMC() {}
};

// Max allowed bonds to any atom
static const uint MAX_BONDS = 6;
// Factory function, determines, prepares and returns appropriate CBMC
CBMC *MakeCBMC(System &sys, const Forcefield &ff, const MoleculeKind &kind,
               const Setup &set);
} // namespace cbmc

#endif
