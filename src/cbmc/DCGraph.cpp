/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "DCGraph.h"

#include <cassert>
#include <map>

#include "DCCrankShaftAng.h"
#include "DCCrankShaftDih.h"
#include "DCFreeHedron.h"
#include "DCFreeHedronSeed.h"
#include "DCLinkedHedron.h"
#include "DCRotateCOM.h"

namespace cbmc {
DCGraph::DCGraph(System &sys, const Forcefield &ff, const MoleculeKind &kind,
                 const Setup &set)
    : data(sys, ff, set) {
  using namespace mol_setup;
  MolMap::const_iterator it = set.mol.kindMap.find(kind.uniqueName);
  assert(it != set.mol.kindMap.end());
  const MolKind setupKind = it->second;

  idExchange = new DCRotateCOM(&data, setupKind);
  // init the coordinate
  coords.Init(setupKind.atoms.size());

  std::vector<uint> atomToNode(setupKind.atoms.size(), 0);
  std::vector<uint> bondCount(setupKind.atoms.size(), 0);
  // Count the number of bonds for each atom
  for (uint b = 0; b < setupKind.bonds.size(); ++b) {
    const Bond &bond = setupKind.bonds[b];
    ++bondCount[bond.a0];
    ++bondCount[bond.a1];
  }

  // Find the node (number of bound > 1)
  // Construct the starting node (DCFreeHedron)
  // Construct the Linking node (DCLinkHedron)
  for (uint atom = 0; atom < setupKind.atoms.size(); ++atom) {
    if (bondCount[atom] < 2) {
      atomToNode[atom] = -1;
    } else {
      // Get the information of other Atoms that are bonded to the atom
      std::vector<Bond> bonds = AtomBonds(setupKind, atom);
      atomToNode[atom] = nodes.size();
      // Add atom to the node list and initialize it with DCFreeHedron, atom and
      // the first partner of the atom
      nodes.push_back(Node());
      Node &node = nodes.back();
      // Atoms bonded to atom will be build from focus (atom) in random loc.
      node.starting = new DCFreeHedron(&data, setupKind, atom, bonds[0].a1);
      // Atoms bonded to atom will be build from focus (atom) in specified loc.
      node.restarting =
          new DCFreeHedronSeed(&data, setupKind, atom, bonds[0].a1);
      // set the atom index of the node
      node.atomIndex = atom;
      // Loop through all the bonds
      for (uint i = 0; i < bonds.size(); ++i) {
        uint partner = bonds[i].a1;
        // Store partner index for each node
        node.partnerIndex.push_back(partner);
        if (bondCount[partner] == 1) {
          continue;
        }
        // Add partner to the edge list of node and initialize it with partner
        // and the atom in DCLinkedHedron
        // Atoms will be build from prev(atom) to focus (partner)
        Edge e =
            Edge(partner, new DCLinkedHedron(&data, setupKind, partner, atom));
        node.edges.push_back(e);
      }
    }
  }

  // reassign destination values from atom indices to node indices
  for (uint i = 0; i < nodes.size(); ++i) {
    for (uint j = 0; j < nodes[i].edges.size(); ++j) {
      int &dest = nodes[i].edges[j].destination;
      dest = atomToNode[dest];
      assert(dest != -1);
    }
  }

  InitCrankShaft(setupKind);
}

void DCGraph::InitCrankShaft(const mol_setup::MolKind &kind) {
  using namespace mol_setup;
  std::vector<uint> bondCount(kind.atoms.size(), 0);
  std::vector<Bond> allBonds = BondsAll(kind);
  // Count the number of bonds for each atom
  for (uint b = 0; b < allBonds.size(); ++b) {
    ++bondCount[allBonds[b].a0];
    ++bondCount[allBonds[b].a1];
  }

  // Start with atoms that form dihedral
  std::vector<Dihedral> dihs = DihsAll(kind);
  for (uint d = 0; d < dihs.size(); d++) {
    // find the last atom index in the dihedral
    uint a0 = dihs[d].a0;
    uint a1 = dihs[d].a1;
    uint a2 = dihs[d].a2;
    uint a3 = dihs[d].a3;
    // ignore single bonded atoms
    if (bondCount[a0] == 1 && bondCount[a3] == 1) {
      continue;
    }

    bool fixAngle = false;
    // Find all the angle that forms x-a0-a1
    std::vector<Angle> angle = AtomMidEndAngles(kind, a0, a1);
    // Find all the angle that forms a2-a3-x
    std::vector<Angle> tempAng = AtomMidEndAngles(kind, a3, a2);
    // merge all the angle
    angle.insert(angle.end(), tempAng.begin(), tempAng.end());
    // Check to see if any of these angles are fixed or not.
    for (uint a = 0; a < angle.size(); a++) {
      if (data.ff.angles->AngleFixed(angle[a].kind)) {
        fixAngle = true;
      }
    }

    // If there was no fix angles, we create DCCrankShaftDih
    if (!fixAngle) {
      shaftNodes.push_back(new DCCrankShaftDih(&data, kind, a0, a1, a2, a3));
    }
  }

  // Continue with the atoms that form angles.
  std::vector<Angle> angles = AngsAll(kind);
  for (uint a = 0; a < angles.size(); a++) {
    // find the last atom index in the angle
    uint a0 = angles[a].a0;
    uint a1 = angles[a].a1;
    uint a2 = angles[a].a2;
    // ignore single bonded atoms
    if (bondCount[a0] == 1 && bondCount[a2] == 1) {
      continue;
    }

    bool fixAngle = false;
    // Find all the angle that forms x-a0-a1
    std::vector<Angle> angle = AtomMidEndAngles(kind, a0, a1);
    // Find all the angle that forms a1-a2-x
    std::vector<Angle> tempAng = AtomMidEndAngles(kind, a2, a1);
    // merge all the angle
    angle.insert(angle.end(), tempAng.begin(), tempAng.end());
    // Check to see if any of these angles are fixed or not.
    for (uint i = 0; i < angle.size(); i++) {
      if (data.ff.angles->AngleFixed(angle[i].kind)) {
        fixAngle = true;
      }
    }

    // If there was no fix angles, we create DCCrankShaftAngle
    if (!fixAngle) {
      shaftNodes.push_back(new DCCrankShaftAng(&data, kind, a0, a1, a2));
    }
  }

  hasCrankShaft = (shaftNodes.size() != 0);
}

void DCGraph::CrankShaft(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  if (!hasCrankShaft) {
    // No crank shaft move for molecule with less than 4 nodes.
    // Instead we perform Regrowth move within the same box
    Regrowth(oldMol, newMol, molIndex);
  } else {
    // Set coords to coordinate of actual molecule, it will be modified
    // No need to unwrap because box is same.
    oldMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
    newMol.SetCoords(coords, 0);
    // Pick a random node pair
    uint pick = data.prng.randIntExc(shaftNodes.size());
    shaftNodes[pick]->PrepareNew(newMol, molIndex);
    shaftNodes[pick]->BuildNew(newMol, molIndex);
    shaftNodes[pick]->PrepareOld(oldMol, molIndex);
    shaftNodes[pick]->BuildOld(oldMol, molIndex);
  }
}

void DCGraph::Build(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  // Randomly pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  // Visiting the node
  visited[current] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  BuildEdges(oldMol, newMol, molIndex, current);
}

void DCGraph::BuildEdges(TrialMol &oldMol, TrialMol &newMol, uint molIndex,
                         const uint cur) {
  uint current = cur;
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].component;
    // Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    // travel to new node, remove traversed edge
    // Current node is the edge that we picked
    current = fringe[pick].destination;
    // Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    // Visiting the node
    visited[current] = true;

    // add edges to unvisited nodes
    for (uint i = 0; i < nodes[current].edges.size(); ++i) {
      Edge &e = nodes[current].edges[i];
      if (!visited[e.destination]) {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::Regrowth(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  // Randomly pick a node to keep it fixed and not grow it
  int current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  // Visiting the node
  visited[current] = true;
  // Copy the current node's focus coordinate
  uint seedInx = nodes[current].atomIndex;
  newMol.AddAtom(seedInx, oldMol.AtomPosition(seedInx));
  oldMol.ConfirmOldAtom(seedInx);
  // check if we want to grow all atoms from node's focus or not
  bool growAll = data.prng() < 1.0 / nodes.size();
  if (growAll) {
    DCComponent *comp = nodes[current].restarting;
    // Call DCFreeHedronSeed to build all Atoms connected to the node's focus
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);
    // Build all edges
    BuildEdges(oldMol, newMol, molIndex, current);
  } else {
    // Copy the all atoms bonded to node's focus
    for (uint b = 0; b < nodes[current].partnerIndex.size(); b++) {
      uint partner = nodes[current].partnerIndex[b];
      newMol.AddAtom(partner, oldMol.AtomPosition(partner));
      oldMol.ConfirmOldAtom(partner);
    }
    // First we pick a edge that will be fixed and copy the coordinates
    // We continue the same until only one edge left from this node
    // If current is the terminal node, we don't enter the while loop
    // Then continue to build the rest of the molecule from current

    // Copy the edges of the node to currFringe
    currFringe = nodes[current].edges;
    while (currFringe.size() > 1) {
      // randomely pick one of the edges connected to current (fixNode)
      uint pickFixEdg = data.prng.randIntExc(currFringe.size());
      // Travel to picked edges and make it as new fixNode
      uint fixNode = currFringe[pickFixEdg].destination;
      visited[fixNode] = true;
      // Copy the all atoms bonded to fixNode's focus
      for (uint b = 0; b < nodes[fixNode].partnerIndex.size(); b++) {
        uint partner = nodes[fixNode].partnerIndex[b];
        newMol.AddAtom(partner, oldMol.AtomPosition(partner));
        oldMol.ConfirmOldAtom(partner);
      }
      // Copy the edges of the new node to fringe
      fringe = nodes[fixNode].edges;
      // remove the edge that we traveled from
      for (uint f = 0; f < fringe.size(); f++) {
        if (fringe[f].destination == current)
          fringe.erase(fringe.begin() + f);
      }
      // Continue along picked edges and copy the coordinates
      while (!fringe.empty()) {
        fixNode = fringe[0].destination;
        // Copy the all atoms bonded to fixNode's focus
        for (uint b = 0; b < nodes[fixNode].partnerIndex.size(); b++) {
          uint partner = nodes[fixNode].partnerIndex[b];
          newMol.AddAtom(partner, oldMol.AtomPosition(partner));
          oldMol.ConfirmOldAtom(partner);
        }
        // Travel to new fixNode, remove traversed edge
        fringe[0] = fringe.back();
        fringe.pop_back();
        visited[fixNode] = true;
        // Add edges to unvisited nodes
        for (uint i = 0; i < nodes[fixNode].edges.size(); ++i) {
          Edge &e = nodes[fixNode].edges[i];
          if (!visited[e.destination]) {
            fringe.push_back(e);
          }
        }
      }
      // Remove the fixed edge from currFring
      currFringe.erase(currFringe.begin() + pickFixEdg);
    }
    // Now Start building the rest of the molecule from current
    // Start with only one left edge
    // Advance along edges, building as we go
    while (!currFringe.empty()) {
      // Randomly pick one of the edges connected to node
      uint pick = data.prng.randIntExc(currFringe.size());
      DCComponent *comp = currFringe[pick].component;
      // Call DCLinkedHedron and build all Atoms connected to selected edge
      comp->PrepareNew(newMol, molIndex);
      comp->BuildNew(newMol, molIndex);
      comp->PrepareOld(oldMol, molIndex);
      comp->BuildOld(oldMol, molIndex);
      current = currFringe[pick].destination;
      // Remove the edge that we visited
      currFringe[pick] = currFringe.back();
      currFringe.pop_back();
      visited[current] = true;
      for (uint i = 0; i < nodes[current].edges.size(); ++i) {
        Edge &e = nodes[current].edges[i];
        if (!visited[e.destination]) {
          currFringe.push_back(e);
        }
      }
    }
  }
}

void DCGraph::BuildIDNew(TrialMol &newMol, uint molIndex) {
  idExchange->PrepareNew(newMol, molIndex);
  idExchange->BuildNew(newMol, molIndex);
}

void DCGraph::BuildIDOld(TrialMol &oldMol, uint molIndex) {
  idExchange->PrepareOld(oldMol, molIndex);
  idExchange->BuildOld(oldMol, molIndex);
}

void DCGraph::BuildOld(TrialMol &oldMol, uint molIndex) {
  // Randomly pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  // Visiting the node
  visited[current] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomly pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].component;
    // Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    // Travel to new node, remove traversed edge
    // Current node is the edge that we picked
    current = fringe[pick].destination;
    // Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    // Visiting the node
    visited[current] = true;

    // Add edges to unvisited nodes
    for (uint i = 0; i < nodes[current].edges.size(); ++i) {
      Edge &e = nodes[current].edges[i];
      if (!visited[e.destination]) {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::BuildNew(TrialMol &newMol, uint molIndex) {
  // Randomly pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  // Visiting the node
  visited[current] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomly pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].component;
    // Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);

    // Travel to new node, remove traversed edge
    // Current node is the edge that we picked
    current = fringe[pick].destination;
    // Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    // Visiting the node
    visited[current] = true;

    // Add edges to unvisited nodes
    for (uint i = 0; i < nodes[current].edges.size(); ++i) {
      Edge &e = nodes[current].edges[i];
      if (!visited[e.destination]) {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::BuildGrowOld(TrialMol &oldMol, uint molIndex) {
  visited.assign(nodes.size(), false);
  // Use backbone atom to start the node
  int current = -1;
  for (int i = 0; i < (int)nodes.size(); i++) {
    if (nodes[i].atomIndex == oldMol.GetAtomBB(0)) {
      current = i;
      break;
    }
  }

  if (current == -1) {
    std::cout << "Error: In MEMC-3 move, atom "
              << oldMol.GetKind().atomNames[oldMol.GetAtomBB(0)] << " in "
              << oldMol.GetKind().name << " must be a node.\n";
    std::cout << "       This atom must be bounded to two or more atoms! \n";
    exit(1);
  }

  // Visiting the node
  visited[current] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomly pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].component;
    // Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareOld(oldMol, molIndex);
    comp->BuildOld(oldMol, molIndex);

    // Travel to new node, remove traversed edge
    // Current node is the edge that we picked
    current = fringe[pick].destination;
    // Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    // Visiting the node
    visited[current] = true;

    // Add edges to unvisited nodes
    for (uint i = 0; i < nodes[current].edges.size(); ++i) {
      Edge &e = nodes[current].edges[i];
      if (!visited[e.destination]) {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::BuildGrowNew(TrialMol &newMol, uint molIndex) {
  visited.assign(nodes.size(), false);
  // Use backbone atom to start the node
  int current = -1;
  for (int i = 0; i < (int)nodes.size(); i++) {
    if (nodes[i].atomIndex == newMol.GetAtomBB(0)) {
      current = i;
      break;
    }
  }

  if (current == -1) {
    std::cout << "Error: In MEMC-3 move, atom "
              << newMol.GetKind().atomNames[newMol.GetAtomBB(0)] << " in "
              << newMol.GetKind().name << " must be a node.\n";
    std::cout << "       This atom must be bounded to two or more atoms! \n";
    exit(1);
  }

  // Visiting the node
  visited[current] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomly pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].component;
    // Call DCLinkedHedron and build all Atoms connected to selected edge
    comp->PrepareNew(newMol, molIndex);
    comp->BuildNew(newMol, molIndex);

    // Travel to new node, remove traversed edge
    // Current node is the edge that we picked
    current = fringe[pick].destination;
    // Remove the edge that we visited
    fringe[pick] = fringe.back();
    fringe.pop_back();
    // Visiting the node
    visited[current] = true;

    // Add edges to unvisited nodes
    for (uint i = 0; i < nodes[current].edges.size(); ++i) {
      Edge &e = nodes[current].edges[i];
      if (!visited[e.destination]) {
        fringe.push_back(e);
      }
    }
  }
}

void DCGraph::BuildGrowInCav(TrialMol &oldMol, TrialMol &newMol,
                             uint molIndex) {
  visited.assign(nodes.size(), false);
  // Get the seedIndex
  int sIndex;
  if (newMol.HasCav()) {
    sIndex = newMol.GetGrowingAtomIndex();
  } else if (oldMol.HasCav()) {
    sIndex = oldMol.GetGrowingAtomIndex();
  } else {
    std::cout << "Error: Calling BuildGrowInCav, but there is no cavity"
              << " defined for newMol and oldMol.\n";
    exit(EXIT_FAILURE);
  }

  // Use backbone atom to start the node
  int current = -1;
  for (int i = 0; i < (int)nodes.size(); i++) {
    if ((int)nodes[i].atomIndex == sIndex) {
      current = i;
      break;
    }
  }

  if (current == -1) {
    std::cout << "Error: In TargetedSwap or IntraTargetedSwap move, atom "
              << newMol.GetKind().atomNames[sIndex] << " in "
              << newMol.GetKind().name << " must be a node.\n";
    std::cout << "       This atom must be bounded to two or more atoms! \n";
    exit(1);
  }

  // Visiting the node
  visited[current] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  BuildEdges(oldMol, newMol, molIndex, current);
}

DCGraph::~DCGraph() {
  delete idExchange;
  for (uint v = 0; v < nodes.size(); ++v) {
    Node &node = nodes[v];
    delete node.starting;
    delete node.restarting;
    for (uint e = 0; e < node.edges.size(); ++e) {
      delete node.edges[e].component;
    }
  }

  for (uint i = 0; i < shaftNodes.size(); i++) {
    delete shaftNodes[i];
  }
}

} // namespace cbmc
