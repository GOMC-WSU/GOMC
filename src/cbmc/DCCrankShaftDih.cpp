/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCCrankShaftDih.h"
#include "TrialMol.h"
#include "DCData.h"
#include "PRNG.h"
#include "MolSetup.h"
#include "Forcefield.h"
#include "NumLib.h"
#include "CalculateEnergy.h"
#include "XYZArray.h"
#include <numeric>
#include <cassert>

namespace cbmc
{
  struct FindA1 {
  FindA1(uint x) : x(x) {};
  bool operator()(const mol_setup::Bond& b) {
    return (b.a1 == x);
  }
  uint x;
};

DCCrankShaftDih::DCCrankShaftDih(DCData* data, const mol_setup::MolKind& kind,
                                uint a0, uint a1, uint a2, uint a3) : 
                                data(data), a0(a0), a1(a1), a2(a2), a3(a3)
{
  using namespace mol_setup;
  using namespace std;
  vector<bool> visited(kind.atoms.size(), false);
  //Find all the atoms that bonds with atoms a1
  vector<Bond> bonds = AtomBonds(kind, a1);
  //Remove the a0-a1 bond
  bonds.erase(remove_if(bonds.begin(), bonds.end(), FindA1(a0)), bonds.end());

  //Loop through atoms that are bonded to a1
  for(uint b = 0; b < bonds.size(); b++) {
    //Store the atom index if it doesnot exist
    if(!visited[bonds[b].a0]) {
      atoms.push_back(bonds[b].a0);
      visited[bonds[b].a0] = true;
    }

    vector<Bond> temp = AtomBonds(kind, bonds[b].a1);
    //Remove a2-a3 bonds
    temp.erase(remove_if(temp.begin(), temp.end(), FindA1(a3)), temp.end());
    for(uint i = 0; i < temp.size(); i++) {
      if(!visited[temp[i].a0]) {
        bonds.push_back(temp[i]);
      }
    }
  }

  numAtom = atoms.size();

  multiPosRotions = new XYZArray[numAtom];
  for(uint i = 0; i < numAtom; ++i) {
    multiPosRotions[i] = XYZArray(data->nLJTrialsNth);
  }

  if(data->nLJTrialsNth < 1) {
      std::cout << "Error: CBMC secondary atom trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

void DCCrankShaftDih::PrepareOld(TrialMol& oldMol, uint molIndex)
{
  for(uint i = 0; i < numAtom; i++) {
    //Unwrap the coordinates with respect to a0.
    XYZ coord = data->axes.UnwrapPBC(oldMol.AtomPosition(atoms[i]), oldMol.GetBox(),
                            oldMol.AtomPosition(a0));
    multiPosRotions[i].Set(0, coord);
    
  }
}

void DCCrankShaftDih::PrepareNew(TrialMol& newMol, uint molIndex)
{
  for(uint i = 0; i < numAtom; i++) {
    //Unwrap the coordinates with respect to a0.
    XYZ coord = data->axes.UnwrapPBC(newMol.AtomPosition(atoms[i]), newMol.GetBox(),
                            newMol.AtomPosition(a0));
    multiPosRotions[i].Set(0, coord);
    
  }
}

    
void DCCrankShaftDih::BuildOld(TrialMol& oldMol, uint molIndex)
{
  PRNG& prng = data->prng;
  XYZArray& positions = data->positions;
  uint nLJTrials = data->nLJTrialsFirst;
  double* inter = data->inter;
  double* real = data->real;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);

  if(oldMol.COMFix()) {
    nLJTrials = 1;
  } else {
    prng.FillWithRandom(positions, nLJTrials, data->axes, oldMol.GetBox());
  }
  positions.Set(0, data->axes.WrapPBC(oldMol.AtomPosition(atom), oldMol.GetBox()));
  data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                           oldMol.GetBox(), nLJTrials);

  for (uint trial = 0; trial < nLJTrials; ++trial) {
    stepWeight += exp(-1 * data->ff.beta *
                      (inter[trial] + real[trial]));
  }
  oldMol.MultWeight(stepWeight / nLJTrials);
  oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0],
                          0.0, 0.0, 0.0));
  oldMol.ConfirmOldAtom(atom);
}

void DCCrankShaftDih::BuildNew(TrialMol& newMol, uint molIndex)
{
  PRNG& prng = data->prng;
  XYZArray& positions = data->positions;
  uint nLJTrials = data->nLJTrialsFirst;
  double* inter = data->inter;
  double* real = data->real;
  double* ljWeights = data->ljWeights;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);

  if(newMol.COMFix()) {
    nLJTrials = 1;
    positions.Set(0, data->axes.WrapPBC(newMol.GetCavityCenter(), newMol.GetBox()));
  } else {
    prng.FillWithRandom(positions, nLJTrials, data->axes, newMol.GetBox());
  }
  data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                           newMol.GetBox(), nLJTrials);

  double stepWeight = 0;
  for (uint trial = 0; trial < nLJTrials; ++trial) {
    ljWeights[trial] = exp(-1 * data->ff.beta *
                           (inter[trial] + real[trial]));
    stepWeight += ljWeights[trial];
  }
  uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  newMol.MultWeight(stepWeight / nLJTrials);
  newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner],
                          0.0, 0.0, 0.0));
  newMol.AddAtom(atom, positions[winner]);
}
}
