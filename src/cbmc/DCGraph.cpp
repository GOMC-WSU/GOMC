/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCGraph.h"
#include "DCFreeHedron.h"
#include "DCLinkedHedron.h"
#include "MolSetup.h"
#include "Setup.h"
#include "MoleculeKind.h"
#include <cassert>
#include <map>

namespace cbmc
{
DCGraph::DCGraph(System& sys, const Forcefield& ff,
                 const MoleculeKind& kind, const Setup& set)
  : data(sys, ff, set)
{
  using namespace mol_setup;
  MolMap::const_iterator it = set.mol.kindMap.find(kind.name);
  assert(it != set.mol.kindMap.end());
  const MolKind setupKind = it->second;

  std::vector<uint> atomToNode(setupKind.atoms.size(), 0);
  std::vector<uint> bondCount(setupKind.atoms.size(), 0);
  for (uint b = 0; b < setupKind.bonds.size(); ++b) {
    const Bond& bond = setupKind.bonds[b];
    ++bondCount[bond.a0];
    ++bondCount[bond.a1];
  }
  for (uint atom = 0; atom < setupKind.atoms.size(); ++atom) {
    if (bondCount[atom] < 2) {
      atomToNode[atom] = -1;
    } else {
      std::vector<Bond> bonds = AtomBonds(setupKind, atom);
      atomToNode[atom] = nodes.size();
      nodes.push_back(Node());
      Node& node = nodes.back();
      node.starting = new DCFreeHedron(&data, setupKind, atom,
                                       bonds[0].a1);
      for (uint i = 0; i < bonds.size(); ++i) {
        uint partner = bonds[i].a1;
        if(bondCount[partner] == 1) {
          continue;
        }
        Edge e = Edge(partner, new DCLinkedHedron(&data, setupKind, partner, atom));
        node.edges.push_back(e);
      }
    }
  }
  //reassign destination values from atom indices to node indices
  for (uint i = 0; i < nodes.size(); ++i) {
    for (uint j = 0; j < nodes[i].edges.size(); ++j) {
      uint& dest = nodes[i].edges[j].destination;
      dest = atomToNode[dest];
      assert(dest != -1);
    }
  }
}

void DCGraph::Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  visited[current] = true;
  fringe = nodes[current].edges;
  DCComponent* comp = nodes[current].starting;
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  //Advance along edges, building as we go
  while (!fringe.empty()) {
    uint pick = data.prng.randIntExc(fringe.size());
    comp = fringe[pick].component;
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    //travel to new node, remove traversed edge
    current = fringe[pick].destination;
    fringe[pick] = fringe.back();
    fringe.pop_back();
    visited[current] = true;
    //add edges to unvisited nodes
    for(uint i = 0; i < nodes[current].edges.size(); ++i) {
      Edge& e = nodes[current].edges[i];
      if(!visited[e.destination]) {
        fringe.push_back(e);
      }
    }
  }
}

DCGraph::~DCGraph()
{
  for(uint v = 0; v < nodes.size(); ++v) {
    Node& node = nodes[v];
    delete node.starting;
    for(uint e = 0; e < node.edges.size(); ++ e) {
      delete node.edges[e].component;
    }
  }
}

}
