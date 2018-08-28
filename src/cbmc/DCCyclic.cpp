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
#include "DCFreeCycleSeed.h"
#include "DCRotateCOM.h"
#include "DCCrankShaftDih.h"
#include "DCCrankShaftAng.h"
#include "FloydWarshallCycle.h"
#include <cassert>
#include <map>

namespace cbmc
{
DCCyclic::DCCyclic(System& sys, const Forcefield& ff,
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
  isRing.resize(setupKind.atoms.size(), false);
  ringIdx.resize(setupKind.atoms.size(), -1);
  FloydWarshallCycle fwc(setupKind.atoms.size());
  //Count the number of bonds for each atom
  for (uint b = 0; b < setupKind.bonds.size(); ++b) {
    const Bond& bond = setupKind.bonds[b];
    ++bondCount[bond.a0];
    ++bondCount[bond.a1];
    fwc.AddEdge(bond.a0, bond.a1);
  }
  cyclicAtoms = fwc.GetAllUniqueCyclesAndCommonCycles();

  //Find the node (number of bound > 1)
  //Construct the starting node (DCFreeHedron or DCFreeCycle)
  //Construct the Linking node (DCLinkHedron or DCLinkedCycle or DCCloseCycle)
  for (uint atom = 0; atom < setupKind.atoms.size(); ++atom) {
    if (bondCount[atom] < 2) {
      atomToNode[atom] = -1;
      isRing[atom] = false;
      ringIdx[atom] = -1;
    } else {
      //Get the information of other Atoms that are bonded to the atom
      std::vector<Bond> bonds = AtomBonds(setupKind, atom);
      atomToNode[atom] = nodes.size();
      //Add atom to the node list and initialize it with DCFreeHedron, atom and
      // the first partner of the atom
      nodes.push_back(Node());
      Node& node = nodes.back();
      for (uint i = 0; i < cyclicAtoms.size(); i++) {
        if (std::find(cyclicAtoms[i].begin(), cyclicAtoms[i].end(), atom)
            != cyclicAtoms[i].end()) {
            isRing[atom] = true;
            ringIdx[atom] = i;
            break;
        }
      }
      //Check if the node belong to a ring or not
      if(isRing[atom]) {
        //Atoms bonded to atom will be build from focus (atom) in random loc.
        node.starting = new DCFreeCycle(&data, setupKind, cyclicAtoms[ringIdx[atom]],
                                        atom, bonds[0].a1);
        //Atoms bonded to atom will be build from focus (atom) in specified loc.
        node.restarting = new DCFreeCycleSeed(&data, setupKind, cyclicAtoms[ringIdx[atom]],
                                              atom, bonds[0].a1);
      } else {
        //Atoms bonded to atom will be build from focus (atom) in random loc.
        node.starting = new DCFreeHedron(&data, setupKind, atom,
                                         bonds[0].a1);
        //Atoms bonded to atom will be build from focus (atom) in specified loc.
        node.restarting = new DCFreeHedronSeed(&data, setupKind, atom,
                                               bonds[0].a1);
      }
      
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
        
        //Check to see if the partner belongs to a ring or not
        bool ring = false;
        uint ringIndex = -1;
        for (uint i = 0; i < cyclicAtoms.size(); i++) {
            if (std::find(cyclicAtoms[i].begin(), cyclicAtoms[i].end(), partner)
                != cyclicAtoms[i].end()) {
                ring = true;
                ringIndex = i;
                break;
            }
        }

        if (ring) {
            //Add partner to the edge list of node and initialize it with partner
            //and the atom in DCLinkedHedron or DCLinkedCycle or DCCloseCycle
            //Atoms will be build from prev(atom) to focus(partner)
            Edge e = Edge(partner, new DCLinkedCycle(&data, setupKind, cyclicAtoms[ringIndex],
                                                    partner,atom)); 
            node.edges.push_back(e);
        } else {
            //Add partner to the edge list of node and initialize it with partner
            //and the atom in DCLinkedHedron or DCLinkedCycle or DCCloseCycle
            //Atoms will be build from prev(atom) to focus(partner)
            Edge e = Edge(partner, new DCLinkedHedron(&data, setupKind, partner,atom));  
            node.edges.push_back(e);
        }
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

  InitCrankShaft(setupKind);
}

void DCCyclic::InitCrankShaft(const mol_setup::MolKind& kind)
{
  using namespace mol_setup;
  using namespace std;
  //Start with the atoms that form angles.
  vector<Node> tempNodes = nodes;
  vector<bool> visited(kind.atoms.size(), false);

  while(!tempNodes.empty()) {
    //start from last node, find the atom index of the node
    uint a0 = tempNodes.back().atomIndex;
    //Find the angle that end with a0
    vector<Angle> angles = AtomEndAngles(kind, a0);
    while(!angles.empty()) {
      //find the last atomindex in the angle
      uint a1 = angles.back().a1;
      uint a2 = angles.back().a2;

      if(!(visited[a0] && visited[a1] && visited[a2])) {
        bool fixAngle = false;
        //Find all the angle that forms x-a0-a1
        vector<Angle> angle = AtomMidEndAngles(kind, a0, a1);
        //Find all the angle that forms a1-a2-x
        vector<Angle> tempAng = AtomMidEndAngles(kind, a2, a1);
        //merge all the angle
        angle.insert(angle.end(), tempAng.begin(), tempAng.end());
        //Check to see if any of these angles are fixed or not.
        for(uint a = 0; a < angle.size(); a++) {
          if(data.ff.angles->AngleFixed(angle[a].kind)) {
            fixAngle = true;
          }
        }
        //Check to see if atoms that are bonded to a1 belongs to same ring or not
        bool sameRing = false;
        if(isRing[a1]) {
            //FInd the atoms that are bonded to a1
            vector<Bond> bonds = AtomBonds(kind, a1);
            for(uint b = 0; b < bonds.size(); b++) {
                uint partner = bonds[b].a1;
                if((partner == a0) || (partner == a2)) {
                    continue;
                }
                if(isRing[partner]) {
                    sameRing |= (ringIdx[a1] == ringIdx[partner]);
                }
            }
        }

        //If there was no fix angles and atom a1 and any atom bonded to a1 are not
        // in the same ring, we create DCCrankShaftAngle
        if(!fixAngle && !sameRing) {
          crankshaft.push_back(new DCCrankShaftAng(&data, kind, a0, a1, a2));
        }
        visited[a0] = true;
        visited[a1] = true;
        visited[a2] = true;
      }
      angles.pop_back();	
    }
    tempNodes.pop_back();
  }

  hasCrankShaft = (crankshaft.size() != 0);
}

void DCCyclic::CrankShaft(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  if(!hasCrankShaft) {
    //No crank shaft move means all angles are fixed.
    //Instead we perform IntraSwap move
    Build(oldMol, newMol, molIndex);
  } else {
    //Pick a random node pair
    uint pick = data.prng.randIntExc(crankshaft.size());
    //Call DCCrankShaftAng and rotate a1 node around a0-a2 shaft
    crankshaft[pick]->PrepareNew(newMol, molIndex);
    crankshaft[pick]->BuildNew(newMol, molIndex);
    crankshaft[pick]->PrepareOld(oldMol, molIndex);
    crankshaft[pick]->BuildOld(oldMol, molIndex);
  }
}

void DCCyclic::Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  //Randomely pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  //Visiting the node
  visited[current] = true;
  DCComponent* comp = nodes[current].starting;
  //Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  //Advance along edges, building as we go
  BuildEdges(oldMol, newMol, molIndex, current);
}

void DCCyclic::BuildEdges(TrialMol& oldMol, TrialMol& newMol, uint molIndex,
                         const uint cur)
{
  uint current = cur;
  //Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  //Advance along edges, building as we go
  while (!fringe.empty()) {
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent* comp = fringe[pick].connect;
    //Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    //travel to new node, remove traversed edge
    //Current node is the edge that we picked
    current = fringe[pick].destination;
    //Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    //Visiting the node
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

void DCCyclic::BuildIDNew(TrialMol& newMol, uint molIndex)
{
  idExchange->PrepareNew(newMol, molIndex);
  idExchange->BuildNew(newMol, molIndex);
}

void DCCyclic::BuildIDOld(TrialMol& oldMol, uint molIndex)
{
  idExchange->PrepareOld(oldMol, molIndex);
  idExchange->BuildOld(oldMol, molIndex);
}

DCCyclic::~DCCyclic()
{
  delete idExchange;
  for(uint v = 0; v < nodes.size(); ++v) {
    Node& node = nodes[v];
    delete node.starting;
    delete node.restarting;
    for(uint e = 0; e < node.edges.size(); ++ e) {
      delete node.edges[e].connect;
    }
  }

  for(uint i = 0; i < crankshaft.size(); i++) {
      delete crankshaft[i];
  }
}


}