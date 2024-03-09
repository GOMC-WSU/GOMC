/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "DCCyclic.h"

#include <cassert>
#include <map>

#include "DCCrankShaftAng.h"
#include "DCFreeCycle.h"
#include "DCFreeCycleSeed.h"
#include "DCFreeHedron.h"
#include "DCFreeHedronSeed.h"
#include "DCLinkedCycle.h"
#include "DCLinkedHedron.h"
#include "DCRotateCOM.h"
#include "DCRotateOnAtom.h"
#include "FloydWarshallCycle.h"

namespace cbmc {
DCCyclic::DCCyclic(System &sys, const Forcefield &ff, const MoleculeKind &kind,
                   const Setup &set)
    : data(sys, ff, set) {
  using namespace mol_setup;
  MolMap::const_iterator it = set.mol.kindMap.find(kind.uniqueName);
  assert(it != set.mol.kindMap.end());
  const MolKind setupKind = it->second;
  totAtom = setupKind.atoms.size();

  if (totAtom < 4) {
    std::cout
        << "Error: GOMC does not support cyclic molecule with 3 atoms!\n\n";
    exit(EXIT_FAILURE);
  }

  idExchange = new DCRotateCOM(&data, setupKind);

  // init the coordinate
  coords.Init(totAtom);

  std::vector<uint> atomToNode(totAtom, 0);
  std::vector<uint> bondCount(totAtom, 0);
  isRing.resize(totAtom, false);
  ringIdx.resize(totAtom, -1);
  CircuitFinder CF(totAtom);

  // Count the number of bonds for each atom
  for (uint b = 0; b < setupKind.bonds.size(); ++b) {
    const Bond &bond = setupKind.bonds[b];
    ++bondCount[bond.a0];
    ++bondCount[bond.a1];
    CF.addEdge(bond.a0, bond.a1);
    CF.addEdge(bond.a1, bond.a0);
  }
  cyclicAtoms = CF.GetAllCommonCycles();
  // Find the ringindex that each atom belongs to
  for (uint atom = 0; atom < totAtom; ++atom) {
    if (bondCount[atom] < 2) {
      atomToNode[atom] = -1;
      isRing[atom] = false;
      ringIdx[atom] = -1;
    } else {
      for (uint i = 0; i < cyclicAtoms.size(); i++) {
        if (std::find(cyclicAtoms[i].begin(), cyclicAtoms[i].end(), atom) !=
            cyclicAtoms[i].end()) {
          isRing[atom] = true;
          ringIdx[atom] = i;
          break;
        }
      }
    }
  }

  // Find the node (number of bound > 1)
  // Construct the starting node (DCFreeHedron or DCFreeCycle)
  // Construct the Linking node (DCLinkHedron or DCLinkedCycle or DCCloseCycle)
  for (uint atom = 0; atom < totAtom; ++atom) {
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

      // Check if the node belongs to a ring or not
      if (isRing[atom]) {
        int prev = -1;
        for (uint i = 0; i < bonds.size(); i++) {
          // Use one of the atoms that is in the ring as prev
          if (isRing[bonds[i].a1]) {
            prev = bonds[i].a1;
            break;
          }
        }
        assert(prev != -1);
        // Atoms bonded to atom will be build from focus (atom) in random loc.
        node.starting = new DCFreeCycle(&data, setupKind,
                                        cyclicAtoms[ringIdx[atom]], atom, prev);
        // Atoms bonded to atom will be build from focus (atom) in specified
        // loc.
        node.restarting = new DCFreeCycleSeed(
            &data, setupKind, cyclicAtoms[ringIdx[atom]], atom, prev);
      } else {
        // Atoms bonded to atom will be build from focus (atom) in random loc.
        node.starting = new DCFreeHedron(&data, setupKind, atom, bonds[0].a1);
        // Atoms bonded to atom will be build from focus (atom) in specified
        // loc.
        node.restarting =
            new DCFreeHedronSeed(&data, setupKind, atom, bonds[0].a1);
      }

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

        if (isRing[partner]) {
          // Add partner to the edge list of node and initialize it with partner
          // and the atom in DCLinkedHedron or DCLinkedCycle or DCCloseCycle
          // Atoms will be build from prev(atom) to focus(partner)
          Edge e =
              Edge(partner, new DCLinkedCycle(&data, setupKind,
                                              cyclicAtoms[ringIdx[partner]],
                                              partner, atom));
          node.edges.push_back(e);

        } else {
          // Add partner to the edge list of node and initialize it with partner
          // and the atom in DCLinkedHedron or DCLinkedCycle or DCCloseCycle
          // Atoms will be build from prev(atom) to focus(partner)
          Edge e = Edge(partner,
                        new DCLinkedHedron(&data, setupKind, partner, atom));
          node.edges.push_back(e);
        }
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

void DCCyclic::InitCrankShaft(const mol_setup::MolKind &kind) {
  using namespace mol_setup;
  std::vector<Angle> angles = AngsAll(kind);
  std::vector<uint> bondCount(totAtom, 0);
  std::vector<Bond> allBonds = BondsAll(kind);
  // Count the number of bonds for each atom
  for (uint b = 0; b < allBonds.size(); ++b) {
    ++bondCount[allBonds[b].a0];
    ++bondCount[allBonds[b].a1];
  }

  for (uint a = 0; a < angles.size(); a++) {
    // find the atomindex in the angle
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
    // Check to see if atoms that are bonded to a1 belongs to same ring or not
    bool sameRing = false;
    if (isRing[a1]) {
      // FInd the atoms that are bonded to a1
      std::vector<Bond> bonds = AtomBonds(kind, a1);
      for (uint b = 0; b < bonds.size(); b++) {
        uint partner = bonds[b].a1;
        if ((partner == a0) || (partner == a2)) {
          continue;
        }
        if (isRing[partner]) {
          sameRing |= (ringIdx[a1] == ringIdx[partner]);
        }
      }
    }

    // If there was no fix angles and atom a1 and any atom bonded to a1 are not
    // in the same ring, we create DCCrankShaftAngle
    if (!fixAngle && !sameRing) {
      crankshaft.push_back(new DCCrankShaftAng(&data, kind, a0, a1, a2));
    }
  }

  // find the atoms that attached to the edge of the ring
  for (uint atom = 0; atom < totAtom; ++atom) {
    // If this atom is in the ring
    if (isRing[atom]) {
      // Find all the angle that forms x-atom-x
      std::vector<Angle> angle = AtomMidAngles(kind, atom);
      for (uint a = 0; a < angle.size(); a++) {
        // find the atomindex in the angle
        uint a0 = angle[a].a0;
        uint a1 = angle[a].a1;
        uint a2 = angle[a].a2;
        // If number of bonds are less than 3, there is no atom attached
        if (bondCount[a1] < 3) {
          continue;
        }
        // To be on the edge, both a0 and a2 must be in the ring
        if (isRing[a0] && isRing[a2]) {
          bool fixAngle = false;
          bool sameRing = false;
          // Find the atoms that are bonded to a1
          std::vector<Bond> bonds = AtomBonds(kind, a1);
          for (uint b = 0; b < bonds.size(); b++) {
            uint partner = bonds[b].a1;
            if ((partner == a0) || (partner == a2)) {
              continue;
            }
            if (isRing[partner]) {
              sameRing |= (ringIdx[a1] == ringIdx[partner]);
            }
            // Find all the angle that forms partner-a1-x (x is either a0 or a2)
            std::vector<Angle> ang = AtomMidEndAngles(kind, a1, partner);
            // Check to see if any of these angles are fixed or not.
            for (uint i = 0; i < ang.size(); i++) {
              fixAngle |= data.ff.angles->AngleFixed(ang[i].kind);
            }
          }

          // If there was no fix angles and atom a1 and any atom bonded to a1
          // are not
          // in the same ring, we create DCCrankShaftAngle
          if (!fixAngle && !sameRing) {
            crankshaft.push_back(new DCRotateOnAtom(&data, kind, a0, a1, a2));
          }
        }
      }
    }
  }

  hasCrankShaft = (crankshaft.size() != 0);
}

void DCCyclic::CrankShaft(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  if (hasCrankShaft) {
    // Set tCoords to coordinate of actual molecule, it will be modified
    oldMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
    newMol.SetCoords(coords, 0);
    // Pick a random node pair
    uint pick = data.prng.randIntExc(crankshaft.size());
    // Call DCCrankShaftAng and rotate a1 node around a0-a2 shaft
    crankshaft[pick]->PrepareNew(newMol, molIndex);
    crankshaft[pick]->BuildNew(newMol, molIndex);
    crankshaft[pick]->PrepareOld(oldMol, molIndex);
    crankshaft[pick]->BuildOld(oldMol, molIndex);
  } else {
    // No crank shaft move means all angles are fixed.
    // Instead we perform IntraSwap move
    Build(oldMol, newMol, molIndex);
  }
}

void DCCyclic::Build(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  // Set bCoords to unwrap coordinate of actual molecule
  oldMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
  oldMol.GetAxes().UnwrapPBC(coords, oldMol.GetBox(), oldMol.AtomPosition(0));
  oldMol.SetBCoords(coords, 0);
  newMol.SetBCoords(coords, 0);

  // Randomely pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  destVisited.assign(totAtom, false);
  // Visiting the node
  visited[current] = true;
  destVisited[nodes[current].atomIndex] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  BuildEdges(oldMol, newMol, molIndex, current);
}

void DCCyclic::BuildEdges(TrialMol &oldMol, TrialMol &newMol, uint molIndex,
                          const uint cur) {
  uint current = cur;
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  for (uint i = 0; i < fringe.size(); i++) {
    destVisited[fringe[i].atomIndex] = true;
  }
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].connect;
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
      if (!visited[e.destination] && !destVisited[e.atomIndex]) {
        fringe.push_back(e);
        destVisited[e.atomIndex] = true;
      }
    }
  }
}

void DCCyclic::Regrowth(TrialMol &oldMol, TrialMol &newMol, uint molIndex) {
  // check if we want to grow all atoms from node's focus or not
  bool growAll = data.prng() < 1.0 / nodes.size();

  // Randomely pick a node to keep it fix and not grow it
  int current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  destVisited.assign(totAtom, false);
  // Visiting the node
  visited[current] = true;
  // Copy the current node's focus coordinate
  uint seedInx = nodes[current].atomIndex;
  destVisited[seedInx] = true;

  if (isRing[seedInx] && !growAll) {
    // if selected node was part of ring and we did not grow all, perform
    // crankshaft Avoid the case where we dont have any crankshaft move. Its
    // better to regrowth the whole molecule rather than call IntraSwap move
    if (hasCrankShaft) {
      return CrankShaft(oldMol, newMol, molIndex);
    } else {
      growAll = true;
    }
  }

  // Set bCoords to unwrap coordinate of actual molecule
  oldMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
  oldMol.GetAxes().UnwrapPBC(coords, oldMol.GetBox(), oldMol.AtomPosition(0));
  oldMol.SetBCoords(coords, 0);
  newMol.SetBCoords(coords, 0);

  newMol.AddAtom(seedInx, oldMol.AtomPosition(seedInx));
  oldMol.ConfirmOldAtom(seedInx);
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
    // First we pick a edge that will be fix and continue copy the coordinate
    // We continue the same until only one edge left from this node
    // If current is the terminal node, we dont enter to while loop
    // Then continue to build the rest of the molecule from current

    // Copy the edges of the node to currFringe
    currFringe = nodes[current].edges;
    while (currFringe.size() > 1) {
      // randomely pick one one of the edges connected to fixNode
      uint pickFixEdg = data.prng.randIntExc(currFringe.size());
      // Travel to picked edges and make it as new fixNode
      uint fixNode = currFringe[pickFixEdg].destination;
      visited[fixNode] = true;
      destVisited[nodes[fixNode].atomIndex] = true;
      // Copy the all atoms bonded to fixNode's focus
      for (uint b = 0; b < nodes[fixNode].partnerIndex.size(); b++) {
        uint partner = nodes[fixNode].partnerIndex[b];
        newMol.AddAtom(partner, oldMol.AtomPosition(partner));
        oldMol.ConfirmOldAtom(partner);
      }
      // Copy the edges of the new node to fringe
      fringe = nodes[fixNode].edges;
      // remove the edge that we travelled from
      for (uint f = 0; f < fringe.size(); f++) {
        if (fringe[f].destination == current)
          fringe.erase(fringe.begin() + f);
      }

      for (uint i = 0; i < fringe.size(); i++) {
        destVisited[fringe[i].atomIndex] = true;
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
        // Travel to new fixNode, remove traverled edge
        fringe[0] = fringe.back();
        fringe.pop_back();
        visited[fixNode] = true;
        // Add edges to unvisited nodes
        for (uint i = 0; i < nodes[fixNode].edges.size(); ++i) {
          Edge &e = nodes[fixNode].edges[i];
          if (!visited[e.destination] && !destVisited[e.atomIndex]) {
            fringe.push_back(e);
            destVisited[e.atomIndex] = true;
          }
        }
      }
      // Remove the fixed edge from currFring
      currFringe.erase(currFringe.begin() + pickFixEdg);
    }

    for (uint i = 0; i < currFringe.size(); i++) {
      destVisited[currFringe[i].atomIndex] = true;
    }
    // Now Start building the rest of the molecule from current
    // Start with only one left edge
    // Advance along edges, building as we go
    while (!currFringe.empty()) {
      // Randomely pick one of the edges connected to node
      uint pick = data.prng.randIntExc(currFringe.size());
      DCComponent *comp = currFringe[pick].connect;
      // Call DCLinkedCycle and build all Atoms connected to selected edge
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
        if (!visited[e.destination] && !destVisited[e.atomIndex]) {
          currFringe.push_back(e);
          destVisited[e.atomIndex] = true;
        }
      }
    }
  }
}

void DCCyclic::BuildIDNew(TrialMol &newMol, uint molIndex) {
  idExchange->PrepareNew(newMol, molIndex);
  idExchange->BuildNew(newMol, molIndex);
}

void DCCyclic::BuildIDOld(TrialMol &oldMol, uint molIndex) {
  idExchange->PrepareOld(oldMol, molIndex);
  idExchange->BuildOld(oldMol, molIndex);
}

void DCCyclic::BuildOld(TrialMol &oldMol, uint molIndex) {
  // Set bCoords to unwrap coordinate of actual molecule
  oldMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
  // No need to unwrap since the copied coordinates are unwraped
  oldMol.SetBCoords(coords, 0);

  // Randomely pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  destVisited.assign(totAtom, false);
  // Visiting the node
  visited[current] = true;
  destVisited[nodes[current].atomIndex] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  for (uint i = 0; i < fringe.size(); i++) {
    destVisited[fringe[i].atomIndex] = true;
  }
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].connect;
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
      if (!visited[e.destination] && !destVisited[e.atomIndex]) {
        fringe.push_back(e);
        destVisited[e.atomIndex] = true;
      }
    }
  }
}

void DCCyclic::BuildNew(TrialMol &newMol, uint molIndex) {
  // Set bCoords to unwrap coordinate of actual molecule
  newMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
  // No need to unwrap since copied coordinates are unwraped
  newMol.SetBCoords(coords, 0);

  // Randomely pick a node to call DCFreeHedron on it
  uint current = data.prng.randIntExc(nodes.size());
  visited.assign(nodes.size(), false);
  destVisited.assign(totAtom, false);
  // Visiting the node
  visited[current] = true;
  destVisited[nodes[current].atomIndex] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  for (uint i = 0; i < fringe.size(); i++) {
    destVisited[fringe[i].atomIndex] = true;
  }
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].connect;
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
      if (!visited[e.destination] && !destVisited[e.atomIndex]) {
        fringe.push_back(e);
        destVisited[e.atomIndex] = true;
      }
    }
  }
}

void DCCyclic::BuildGrowOld(TrialMol &oldMol, uint molIndex) {
  // Set bCoords to unwrap coordinate of actual molecule
  oldMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
  // No need to unwrap since copied coordinate is unwraped
  oldMol.SetBCoords(coords, 0);

  visited.assign(nodes.size(), false);
  destVisited.assign(totAtom, false);
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
  destVisited[nodes[current].atomIndex] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  for (uint i = 0; i < fringe.size(); i++) {
    destVisited[fringe[i].atomIndex] = true;
  }
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].connect;
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
      if (!visited[e.destination] && !destVisited[e.atomIndex]) {
        fringe.push_back(e);
        destVisited[e.atomIndex] = true;
      }
    }
  }
}

void DCCyclic::BuildGrowNew(TrialMol &newMol, uint molIndex) {
  // Set bCoords to unwrap coordinate of actual molecule
  newMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
  // No need to unwrap since copied coordinate is unwraped
  newMol.SetBCoords(coords, 0);

  visited.assign(nodes.size(), false);
  destVisited.assign(totAtom, false);
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
  destVisited[nodes[current].atomIndex] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  // Advance along edges, building as we go
  // Copy the edges of the node to fringe
  fringe = nodes[current].edges;
  for (uint i = 0; i < fringe.size(); i++) {
    destVisited[fringe[i].atomIndex] = true;
  }
  // Advance along edges, building as we go
  while (!fringe.empty()) {
    // Randomely pick one of the edges connected to node
    uint pick = data.prng.randIntExc(fringe.size());
    DCComponent *comp = fringe[pick].connect;
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
      if (!visited[e.destination] && !destVisited[e.atomIndex]) {
        fringe.push_back(e);
        destVisited[e.atomIndex] = true;
      }
    }
  }
}

void DCCyclic::BuildGrowInCav(TrialMol &oldMol, TrialMol &newMol,
                              uint molIndex) {
  // Set bCoords to unwrap coordinate of actual molecule
  newMol.GetCoords().CopyRange(coords, 0, 0, coords.Count());
  // No need to unwrap since copied coordinate is unwraped and wrapped properly
  newMol.SetBCoords(coords, 0);

  visited.assign(nodes.size(), false);
  destVisited.assign(totAtom, false);
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

  // Use seedIndex atom to start the node
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
  destVisited[nodes[current].atomIndex] = true;
  DCComponent *comp = nodes[current].starting;
  // Call DCFreeHedron to build all Atoms connected to the node
  comp->PrepareNew(newMol, molIndex);
  comp->BuildNew(newMol, molIndex);
  comp->PrepareOld(oldMol, molIndex);
  comp->BuildOld(oldMol, molIndex);
  // Advance along edges, building as we go
  BuildEdges(oldMol, newMol, molIndex, current);
}

DCCyclic::~DCCyclic() {
  delete idExchange;
  for (uint v = 0; v < nodes.size(); ++v) {
    Node &node = nodes[v];
    delete node.starting;
    delete node.restarting;
    for (uint e = 0; e < node.edges.size(); ++e) {
      delete node.edges[e].connect;
    }
  }

  for (uint i = 0; i < crankshaft.size(); i++) {
    delete crankshaft[i];
  }
}

} // namespace cbmc