/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.50
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
double NoEwald::BoxReciprocal(uint box) const
{
  return 0.0;
}

//calculate reciprocate force term for a box with molCoords
void NoEwald::BoxForceReciprocal(XYZArray const& molCoords,
                                 XYZArray& atomForceRec, XYZArray& molForceRec,
                                 uint box)
{
  return;
}

//calculate reciprocate force term for a box
Virial NoEwald::VirialReciprocal(Virial& virial, uint box) const
{
  return virial;
}


//calculate reciprocate term for displacement and rotation move
double NoEwald::MolReciprocal(XYZArray const& molCoords,
                              const uint molIndex, const uint box)
{
  return 0.0;
}


//calculate reciprocate term for lambdaNew and Old with same coordinates
double NoEwald::CFCMCRecip(XYZArray const& molCoords, const double lambdaOld,
                           const double lambdaNew, const uint molIndex,
                           const uint box)
{
  return 0.0;
}


//calculate self term for a box
double NoEwald::BoxSelf(BoxDimensions const& boxAxes, uint box) const
{
  return 0.0;
}


//calculate correction term for a molecule
double NoEwald::MolCorrection(uint molIndex, uint box) const
{
  return 0.0;
}

//It's called in free energy calculation to calculate the change in
// self energy in all lambda states
void NoEwald::ChangeSelf(Energy *energyDiff, Energy &dUdL_Coul,
                         const std::vector<double> &lambda_Coul,
                         const uint iState, const uint molIndex,
                         const uint box) const
{
  return;
}

//It's called in free energy calculation to calculate the change in
// correction energy in all lambda states
void NoEwald::ChangeCorrection(Energy *energyDiff, Energy &dUdL_Coul,
                               const std::vector<double> &lambda_Coul,
                               const uint iState, const uint molIndex,
                               const uint box) const
{
  return;
}

//It's called in free energy calculation to calculate the change in
// reciprocal energy in all lambda states
void NoEwald::ChangeRecip(Energy *energyDiff, Energy &dUdL_Coul,
                          const std::vector<double> &lambda_Coul,
                          const uint iState, const uint molIndex,
                          const uint box) const
{
  return;
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
