/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCGRAPH_H
#define DCGRAPH_H
#include <utility>
#include <vector>

#include "CBMC.h"
#include "DCComponent.h"
#include "DCData.h"
#include "MolSetup.h"
#include "MoleculeKind.h"
#include "Setup.h"

/*CBMC graph of a branched/cyclic molecule
 * The Decoupled/Coupled CBMC algorithm is represented by
 * traversing a spanning tree of the graph.
 */

namespace cbmc {

class DCComponent;

class DCGraph : public CBMC {
public:
  DCGraph(System &sys, const Forcefield &ff, const MoleculeKind &kind,
          const Setup &set);

  void Build(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void Regrowth(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void CrankShaft(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  void BuildEdges(TrialMol &oldMol, TrialMol &newMol, uint molIndex,
                  const uint current);
  // used in MEMC moves
  void BuildIDNew(TrialMol &newMol, uint molIndex);
  void BuildIDOld(TrialMol &oldMol, uint molIndex);
  void BuildNew(TrialMol &newMol, uint molIndex);
  void BuildOld(TrialMol &oldMol, uint molIndex);
  void BuildGrowNew(TrialMol &newMol, uint molIndex);
  void BuildGrowOld(TrialMol &oldMol, uint molIndex);
  // used in TargetedSwap
  void BuildGrowInCav(TrialMol &oldMol, TrialMol &newMol, uint molIndex);
  ~DCGraph();

private:
  // Find the two nodes that are forming dihedral or angle and initialize it.
  void InitCrankShaft(const mol_setup::MolKind &kind);
  // Store edge's atom that are connected to node and has more than 1 bond
  // Each edge is a node as well
  struct Edge {
    int destination; // destination is partner node index.
    DCComponent *component;
    Edge(uint d, DCComponent *c) : destination(d), component(c) {}
  };

  // Store the branching atom and all Atoms that are connected to this
  // branching atom
  struct Node {
    uint atomIndex;
    // starting will be initialized with DCFreeHedron using random loc
    DCComponent *starting;
    // starting will be initialized with DCFreeHedronSeed, using specify seed
    DCComponent *restarting;
    // all the atom that are connected to this node and has more than 1 bond
    // will be in edges and initialized with DCLinkedHedron
    std::vector<Edge> edges;
    std::vector<uint> partnerIndex;
  };

  DCComponent *idExchange;
  DCData data;
  bool hasCrankShaft;
  std::vector<Node> nodes;
  std::vector<Edge> fringe, currFringe;
  std::vector<bool> visited;
  std::vector<DCComponent *> shaftNodes;
  XYZArray coords;
};
} // namespace cbmc

#endif
