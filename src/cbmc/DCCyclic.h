/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCCYCLIC_H
#define DCCYCLIC_H
#include <utility>
#include <vector>

#include "CBMC.h"
#include "CircuitFinder.h"
#include "DCComponent.h"
#include "DCData.h"
#include "FloydWarshallCycle.h"
#include "MolSetup.h"
#include "MoleculeKind.h"
#include "Setup.h"
/*CBMC forcyclic molecule
 * The Decoupled/Coupled CBMC algorithm is represented by
 * traversing a spanning tree of the graph.
 * body of the ring is keep rigid during the move.
 */

namespace cbmc {

class DCComponent;

class DCCyclic : public CBMC {
public:
  DCCyclic(System &sys, const Forcefield &ff, const MoleculeKind &kind,
           const Setup &set);

  void Build(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void Regrowth(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void CrankShaft(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void BuildEdges(TrialMol &oldMol, TrialMol &newMol, uint molIndex,
                  const uint current);
  // Used in MEMC moves
  void BuildIDNew(TrialMol &newMol, uint molIndex);
  void BuildIDOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildGrowNew(TrialMol &newMol, uint molIndex);
  void BuildGrowOld(TrialMol &oldMol, uint molIndex);
  // used in TargetedSwap
  void BuildGrowInCav(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  ~DCCyclic();

private:
  // Find the two nodes that are forming dihedral or angle and initialize it.
  void InitCrankShaft(const mol_setup::MolKind &kind);
  // Store edge's atom that are connected to node and has more than 1 bond
  // Each edge is a node as well
  struct Edge {
    int destination; // destination is partner node index.
    uint atomIndex;  // atom index of the edge
    // To build the next segment from prev-focus
    DCComponent *connect;
    Edge(uint d, DCComponent *c) : destination(d), atomIndex(d), connect(c) {}
  };

  // Store the branching atom and all Atoms that are connected to this
  // branching atom
  struct Node {
    uint atomIndex;
    // starting will be initialized with DCFreeCycle or DCFreeHedron using
    // random loc
    DCComponent *starting;
    // starting will be initialized with DCFreeCycleSeed or DCFreeHedronSeed,
    // using specify seed
    DCComponent *restarting;
    // all the atom that are connected to this node and has more than 1 bond
    // will be in edges and initialized with DCLinkedCycle or DCLinkedHedron
    std::vector<Edge> edges;
    std::vector<uint> partnerIndex;
  };

  DCComponent *idExchange;
  DCData data;
  bool hasCrankShaft;
  std::vector<bool> isRing;  // To check if atom is belong to a ring
  std::vector<uint> ringIdx; // index to the row of cyclicAtoms
  std::vector<Node> nodes;
  std::vector<Edge> fringe, currFringe;
  std::vector<bool> visited, destVisited;
  std::vector<DCComponent *> crankshaft;
  std::vector<std::vector<int>> cyclicAtoms;
  XYZArray coords;
  uint totAtom;
};
} // namespace cbmc

#endif