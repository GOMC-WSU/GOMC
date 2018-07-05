/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCGraph.h"
#include "DCFreeHedron.h"
#include "DCLinkedHedron.h"
#include "DCFreeHedronSeed.h"
#include "DCRotateCOM.h"
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

  idExchange = new DCRotateCOM(&data, setupKind);

  std::vector<uint> atomToNode(setupKind.atoms.size(), 0);
  std::vector<uint> bondCount(setupKind.atoms.size(), 0);
  //Count the number of bonds for each atom
  for (uint b = 0; b < setupKind.bonds.size(); ++b) {
    const Bond& bond = setupKind.bonds[b];
    ++bondCount[bond.a0];
    ++bondCount[bond.a1];
  }

  //Find the node (number of bound > 1)
  //Construct the starting node (DCFreeHedron)
  //Construct the Linking node (DCLinkHedron)
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
        //and the atom in DCLinkedHedron
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

  InitCrankShaft(setupKind);
}


void DCGraph::InitCrankShaft(mol_setup::MolKind& kind)
{
  using namespace mol_setup;
  using namespace std;
  vector<Node> tempNodes = nodes;
  while(!tempNodes.empty()) {
    //start from last node, find the atom index of the node
    uint a0 = tempNodes.back().atomIndex;
    //Find the dihedrals that end with a0
    vector<Dihedral> dihs = AtomEndDihs(kind, a0);
    while(!dihs.empty()) {
      //find the last atomindex in the dihedral
      uint a3 = dihs.back().a3;
      //Check if the a3 is a node or not
      for(uint n = 0; n < tempNodes.size(); n++) {
	if(tempNodes[n].atomIndex == a3) {
	  shaftNodesDih.push_back(new DCCrankShaftDih(&data, kind, a0, a3));
	}
      }
      dihs.pop_back();	
    }
    tempNodes.pop_back();
  }
}

void DCGraph::Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
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

void DCGraph::BuildEdges(TrialMol& oldMol, TrialMol& newMol, uint molIndex,
                         const uint cur)
{
  uint current = cur;
  //Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  //Advance along edges, building as we go
  while (!fringe.empty()) {
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent* comp = fringe[pick].component;
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

void DCGraph::Regrowth(TrialMol& oldMol, TrialMol& newMol, uint molIndex)
{
  //Randomely pick a node to keep it fix and not grow it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  //Visiting the node
  visited[current] = true;
  //Copy the current node's focus coordinate
  uint seedInx = nodes[current].atomIndex;
  newMol.AddAtom(seedInx, oldMol.AtomPosition(seedInx));
  oldMol.ConfirmOldAtom(seedInx);
  //check if we want to grow all atoms from node's focus or not
  bool growAll = data.prng() < 1.0 / nodes.size();
  if(growAll) {
    DCComponent* comp = nodes[current].restarting;
    //Call DCFreeHedronSeed to build all Atoms connected to the node's focus
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);
    //Build all edges
    BuildEdges(oldMol, newMol, molIndex, current);
  } else {
    //Copy the all atoms bonded to node's focus
    for(uint b = 0; b < nodes[current].partnerIndex.size(); b++) {
      uint partner = nodes[current].partnerIndex[b];
      newMol.AddAtom(partner, oldMol.AtomPosition(partner));
      oldMol.ConfirmOldAtom(partner);
    }

    if(nodes[current].edges.size() == 1) {
      //If current is the terminal node, continue building all edges
      BuildEdges(oldMol, newMol, molIndex, current);
    } else {
      //First we pick a edge and continue copying the coordinate
      //Then continue to build the rest of the molecule from current
      //Copy the edges of the node to fringe
      fringe = nodes[current].edges;
      //randomely pick one one of the edges connected to fixNode
      uint pickFixEdg = data.prng.randIntExc(fringe.size());
      //Travel to picked edges and make it as new fixNode
      uint fixNode = fringe[pickFixEdg].destination;
      visited[fixNode] = true;
      //Copy the all atoms bonded to fixNode's focus
      for(uint b = 0; b < nodes[fixNode].partnerIndex.size(); b++) {
        uint partner = nodes[fixNode].partnerIndex[b];
        newMol.AddAtom(partner, oldMol.AtomPosition(partner));
        oldMol.ConfirmOldAtom(partner);
      }
      //Copy the edges of the new node to fringe
      fringe = nodes[fixNode].edges;
      //remove the edge that we travelled from
      for( uint f = 0; f < fringe.size(); f++) {
        if(fringe[f].destination == current)
          fringe.erase(fringe.begin() + f);
      }
      //Continue along picked edges and copy the coordinates
      while(!fringe.empty()) {
        fixNode = fringe[0].destination;
        //Copy the all atoms bonded to fixNode's focus
        for(uint b = 0; b < nodes[fixNode].partnerIndex.size(); b++) {
          uint partner = nodes[fixNode].partnerIndex[b];
          newMol.AddAtom(partner, oldMol.AtomPosition(partner));
          oldMol.ConfirmOldAtom(partner);
        }
        //Travel to new fixNode, remove traversed edge
        fringe[0] = fringe.back();
        fringe.pop_back();
        visited[fixNode] = true;
        //Add edges to unvisited nodes
        for(uint i = 0; i < nodes[fixNode].edges.size(); ++i) {
          Edge& e = nodes[fixNode].edges[i];
          if(!visited[e.destination]) {
            fringe.push_back(e);
          }
        }
      }
      //Now Start building the rest of the molecule from current
      //Copy the edges of the current node to fringe
      fringe = nodes[current].edges;
      //Remove the fixed edge from fringe
      fringe.erase(fringe.begin() + pickFixEdg);
      //Advance along edges, building as we go
      while(!fringe.empty()) {
        //Randomely pick one of the edges connected to node
        uint pick = data.prng.randIntExc(fringe.size());
        DCComponent* comp = fringe[pick].component;
        //Call DCLinkedHedron and build all Atoms connected to selected edge
        comp->PrepareNew(newMol, molIndex);
        comp->BuildNew(newMol, molIndex);
        comp->PrepareOld(oldMol, molIndex);
        comp->BuildOld(oldMol, molIndex);
        current = fringe[pick].destination;
        //Remove the edge that we visited
        fringe[pick] = fringe.back();
        fringe.pop_back();
        visited[current] = true;
        for(uint i = 0; i < nodes[current].edges.size(); ++i) {
          Edge& e = nodes[current].edges[i];
          if(!visited[e.destination]) {
            fringe.push_back(e);
          }
        }
      }
    }
  }
}

void DCGraph::BuildIDNew(TrialMol& newMol, uint molIndex)
{
  idExchange->PrepareNew(newMol, molIndex);
  idExchange->BuildNew(newMol, molIndex);
}

void DCGraph::BuildIDOld(TrialMol& oldMol, uint molIndex)
{
  idExchange->PrepareOld(oldMol, molIndex);
  idExchange->BuildOld(oldMol, molIndex);
}

void DCGraph::BuildOld(TrialMol& oldMol, uint molIndex)
{
  //Randomely pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  //Visiting the node
  visited[current] = true;
  DCComponent* comp = nodes[current].starting;
  //Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  //Advance along edges, building as we go
  //Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  //Advance along edges, building as we go
  while(!fringe.empty())
  {
    //Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent* comp = fringe[pick].component;
    //Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    //Travel to new node, remove traversed edge
    //Current node is the edge that we picked
    current = fringe[pick].destination;
    //Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    //Visiting the node
    visited[current] = true;

    //Add edges to unvisited nodes
    for(uint i = 0; i < nodes[current].edges.size(); ++i)
    {
      Edge& e = nodes[current].edges[i];
      if(!visited[e.destination])
      {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::BuildNew(TrialMol& newMol, uint molIndex)
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
  //Advance along edges, building as we go
  //Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  //Advance along edges, building as we go
  while(!fringe.empty())
  {
    //Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent* comp = fringe[pick].component;
    //Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);

    //Travel to new node, remove traversed edge
    //Current node is the edge that we picked
    current = fringe[pick].destination;
    //Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    //Visiting the node
    visited[current] = true;

    //Add edges to unvisited nodes
    for(uint i = 0; i < nodes[current].edges.size(); ++i)
    {
      Edge& e = nodes[current].edges[i];
      if(!visited[e.destination])
      {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::BuildGrowOld(TrialMol& oldMol, uint molIndex)
{
  visited.assign(nodes.size(), false);
  //Use backbone atom to start the node
  uint current = -1;
  for(uint i = 0; i < nodes.size(); i++) {
    if(nodes[i].atomIndex == oldMol.GetAtomBB(0)) {
      current = i;
      break;
    }
  }

  if(current == -1) {
    std::cout << "In MEMC-3 move, atom index "<< oldMol.GetAtomBB(0)<<
        " must be a node."
        << "This atom must be bounded to two or more atoms! \n";
    exit(1);
  }
    
  //Visiting the node
  visited[current] = true;
  DCComponent* comp = nodes[current].starting;
  //Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  //Advance along edges, building as we go
  //Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  //Advance along edges, building as we go
  while(!fringe.empty())
  {
    //Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent* comp = fringe[pick].component;
    //Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    //Travel to new node, remove traversed edge
    //Current node is the edge that we picked
    current = fringe[pick].destination;
    //Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    //Visiting the node
    visited[current] = true;

    //Add edges to unvisited nodes
    for(uint i = 0; i < nodes[current].edges.size(); ++i)
    {
      Edge& e = nodes[current].edges[i];
      if(!visited[e.destination])
      {
        fringe.push_back(e);
      }
    }
  }
}


void DCGraph::BuildGrowNew(TrialMol& newMol, uint molIndex)
{
  visited.assign(nodes.size(), false);
  //Use backbone atom to start the node
  uint current = -1;
  for(uint i = 0; i < nodes.size(); i++) {
    if(nodes[i].atomIndex == newMol.GetAtomBB(0)) {
      current = i;
      break;
    }
  }
    
  if(current == -1) {
    std::cout << "In MEMC-3 move, atom index "<< newMol.GetAtomBB(0)<<
        " must be a node."
        << "This atom must be bounded to two or more atoms! \n";
    exit(1);
  }
    
  //Visiting the node
  visited[current] = true;
  DCComponent* comp = nodes[current].starting;
  //Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  //Advance along edges, building as we go
  //Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  //Advance along edges, building as we go
  while(!fringe.empty())
  {
    //Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent* comp = fringe[pick].component;
    //Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);

    //Travel to new node, remove traversed edge
    //Current node is the edge that we picked
    current = fringe[pick].destination;
    //Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    //Visiting the node
    visited[current] = true;

    //Add edges to unvisited nodes
    for(uint i = 0; i < nodes[current].edges.size(); ++i)
    {
      Edge& e = nodes[current].edges[i];
      if(!visited[e.destination])
      {
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
    delete node.restarting;
    for(uint e = 0; e < node.edges.size(); ++ e) {
      delete node.edges[e].component;
    }
  }
  delete idExchange;
}


}
