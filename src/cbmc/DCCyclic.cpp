/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCCyclic.h"
#include "DCFreeHedron.h"
#include "DCLinkedHedron.h"
#include "DCFreeHedronSeed.h"
#include "DCFreeCycle.h"
#include "DCLinkedCycle.h"
#include "DCCloseCycle.h"
#include "DCFreeCycleSeed.h"
#include "DCRotateCOM.h"
#include "DCCrankShaftDih.h"
#include "DCCrankShaftAng.h"
#include "FloydWarshallCycle.h"
#include <cassert>
#include <map>

namespace cbmc
{
DCGraphDCCyclic::DCCyclic(System& sys, const Forcefield& ff,
                         const MoleculeKind& kind, const Setup& set)
  : data(sys, ff, set)
{
  using namespace mol_setup;
  MolMap::const_iterator it = set.mol.kindMap.find(kind.name);
  assert(it != set.mol.kindMap.end());
  const MolKind setupKind = it->second;

  idExchange = new DCRotateCOM(&data, setupKind);

  std::vector<uint> atomToNode(setupKind.atoms.size(), 0);
  std::vector<uint> bondCount(setupKind.atoms.size(), 0);
  FloydWarshallCycle fwc(setupKind.atoms.size());
  //Count the number of bonds for each atom
  for (uint b = 0; b < setupKind.bonds.size(); ++b) {
    const Bond& bond = setupKind.bonds[b];
    ++bondCount[bond.a0];
    ++bondCount[bond.a1];
    fwc.AddEdge(bond.a0, bond.a1);
  }
  cyclicAtom = fwc.GetAllUniqueCyclesAndCommonCycles();

  //Find the node (number of bound > 1)
  //Construct the starting node (DCFreeHedron or DCFreeCycle)
  //Construct the Linking node (DCLinkHedron or DCLinkedCycle)
  for (uint atom = 0; atom < setupKind.atoms.size(); ++atom) {
    if (bondCount[atom] < 2) {
      atomToNode[atom] = -1;
    } else {
      //Get the information of other Atoms that are bonded to the atom
      std::vector<Bond> bonds = AtomBonds(setupKind, atom);
      atomToNode[atom] = nodes.size();
      //Add atom to the node list and initialize it with DCFreeHedron, atom and
      // the first partner of the atom
      nodes.push_back(Node());
      Node& node = nodes.back();
      //Atoms bonded to atom will be build from focus (atom) in random loc.
      node.starting = new DCFreeHedron(&data, setupKind, atom,
                                       bonds[0].a1);
      //Atoms bonded to atom will be build from focus (atom) in specified loc.
      node.restarting = new DCFreeHedronSeed(&data, setupKind, atom,
                                             bonds[0].a1);
      //set the atom index of the node
      node.atomIndex = atom;
      //Loop through all the bonds
      for (uint i = 0; i < bonds.size(); ++i) {
        uint partner = bonds[i].a1;
        //Store partner index for each node
        node.partnerIndex.push_back(partner);
        if(bondCount[partner] == 1) {
          continue;
        }
        //Add partner to the edge list of node and initialize it with partner
        //and the atom in DCLinkedHedron or DCLinkedCycle or DCCloseCycle
        //Atoms will be build from prev(atom) to focus (partner)
        Edge e = Edge(partner, new DCLinkedHedron(&data, setupKind, partner,
                      atom));
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

  //InitCrankShaft(setupKind);
}