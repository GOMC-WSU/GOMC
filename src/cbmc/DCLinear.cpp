/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include <cassert>
#include "DCLinear.h"
#include "DCSingle.h"
#include "DCOnSphere.h"

using namespace cbmc;

DCLinear::DCLinear(System& sys, const Forcefield& ff,
                   const MoleculeKind& kind, const Setup& set) :
  data(sys, ff, set)
{
  mol_setup::MolMap::const_iterator it = set.mol.kindMap.find(kind.name);
  assert(it != set.mol.kindMap.end());
  const mol_setup::MolKind setupKind = it->second;
  uint size = kind.NumAtoms();
  atomSize = size;

  if(atomSize < 3)
  {
    forward.push_back(new DCSingle(&data, 0));
    backward.push_back(new DCSingle(&data, size - 1));

    if(size < 2)
      return;
  
    forward.push_back(new DCOnSphere(&data, setupKind, 1, 0));
    backward.push_back(new DCOnSphere(&data, setupKind, size - 2, size - 1));
  }
  else
  {
    graph = new DCGraph(sys, ff, kind, set);
  }
}

DCLinear::~DCLinear()
{
  if(atomSize < 3){
    for(uint i = 0; i < forward.size(); ++i) {
      delete forward[i];
      delete backward[i];
    }
  }
  else {
    delete graph;
  }
}

void DCLinear::Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  if(atomSize < 3)
  {
    std::vector<DCComponent*>& comps =data.prng.randInt(1) ? forward : backward;
    for(uint i = 0; i < comps.size(); ++i)
    {
      comps[i]->PrepareNew(newMol, molIndex);
      comps[i]->BuildNew(newMol, molIndex);
    }

    for(uint i = 0; i < comps.size(); ++i)
    {
      comps[i]->PrepareOld(oldMol, molIndex);
      comps[i]->BuildOld(oldMol, molIndex);
    }
  }
  else
  {
    graph->Build(oldMol, newMol, molIndex);
  }
}

void DCLinear::Regrowth(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  //perform Intra-Swap move within the same box
  if(atomSize < 2) {
    return Build(oldMol, newMol, molIndex);
  }
  else if(atomSize < 3) {
    //we only have two atoms in molecule: atom 0, 1
    uint fix = data.prng.randInt(1);
    uint grow = 1 - fix;
    //If fix == 0, forward (build atom 1), else backward (build atom 0)
    std::vector<DCComponent*>& comps = fix ? backward : forward;
    
    //copy the coordinate of the fix atom
    newMol.AddAtom(fix, oldMol.AtomPosition(fix));
    oldMol.ConfirmOldAtom(fix);
    //build the second atom
    comps[1]->PrepareNew(newMol, molIndex);
    comps[1]->BuildNew(newMol, molIndex);
    comps[1]->PrepareOld(oldMol, molIndex);
    comps[1]->BuildOld(oldMol, molIndex);
  }
  else {
    graph->Regrowth(oldMol, newMol, molIndex);
  }
}
