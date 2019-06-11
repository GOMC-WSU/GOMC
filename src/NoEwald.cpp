/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "NoEwald.h"


NoEwald::NoEwald(StaticVals & stat, System & sys) :
  Ewald(stat, sys) {}

void NoEwald::Init() {}

void NoEwald::AllocMem()
{
  return;
}

void NoEwald::RecipInit(uint box, BoxDimensions const& boxAxes)
{
  return;
}

//calculate reciprocate term for a box
void NoEwald::BoxReciprocalSetup(uint box, XYZArray const& molCoords)
{
  return;
}


//calculate reciprocate term for a box
real NoEwald::BoxReciprocal(uint box) const
{
  return 0.0;
}


//calculate reciprocate force term for a box
Virial NoEwald::ForceReciprocal(Virial& virial, uint box) const
{
  return virial;
}

//calculate correction force term for a box
Virial NoEwald::ForceCorrection(Virial& virial, uint box) const
{
  return virial;
}


//calculate reciprocate term for displacement and rotation move
real NoEwald::MolReciprocal(XYZArray const& molCoords,
                              const uint molIndex, const uint box)
{
  return 0.0;
}


//calculate self term for a box
real NoEwald::BoxSelf(BoxDimensions const& boxAxes, uint box) const
{
  return 0.0;
}


//calculate correction term for a molecule
real NoEwald::MolCorrection(uint molIndex, uint box) const
{
  return 0.0;
}


//calculate reciprocate term in destination box for swap move
real NoEwald::SwapDestRecip(const cbmc::TrialMol &newMol,
                              const uint box,
                              const int molIndex)
{
  return 0.0;
}


//calculate reciprocate term in source box for swap move
real NoEwald::SwapSourceRecip(const cbmc::TrialMol &oldMol,
                                const uint box, const int molIndex)
{
  return 0.0;
}


//calculate reciprocate term for inserting some molecules (kindA) in destination
// box and removing a molecule (kindB) from destination box
real NoEwald::SwapRecip(const std::vector<cbmc::TrialMol> &newMol,
                          const std::vector<cbmc::TrialMol> &oldMol)
{
  return 0.0;
}


//calculate self term after swap move
real NoEwald::SwapSelf(const cbmc::TrialMol& trialMol) const
{
  return 0.0;
}

//calculate correction term after swap move
real NoEwald::SwapCorrection(const cbmc::TrialMol& trialMol) const
{
  return 0.0;
}

//back up reciptocate value to Ref (will be called during initialization)
void NoEwald::SetRecipRef(uint box)
{
  return;
}

//update reciprocate values
void NoEwald::UpdateRecip(uint box)
{
  return;
}

//update the hx,y,z hsqr and prefact
void NoEwald::UpdateRecipVec(uint box)
{
  return;
}

//restore cosMol and sinMol
void NoEwald::RestoreMol(int molIndex)
{
  return;
}

//restore the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void NoEwald::exgMolCache()
{
  return;
}

//backup the whole cosMolRef & sinMolRef into cosMolBoxRecip & sinMolBoxRecip
void NoEwald::backupMolCache()
{
  return;
}

void NoEwald::UpdateVectorsAndRecipTerms()
{
  return;
}
