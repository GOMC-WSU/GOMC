/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
/* Courtesy of https://www.softwaretestinghelp.com/graph-implementation-cpp/ */

#ifndef BOND_ADJACENCY_LIST_H
#define BOND_ADJACENCY_LIST_H

#include <limits.h>

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <vector>

#include "BasicTypes.h"

// stores adjacency list items
struct adjNode {
  int val;
  int cost;
  adjNode *next;
};
// structure to store edges
struct graphEdge {
  int start_ver, end_ver, weight;
};

class BondAdjacencyList {
public:
  adjNode **head; // adjacency list as array of pointers
  uint nAtoms;    // number of nodes in the graph
  uint nBonds;    // number of edges in the graph

  BondAdjacencyList(FILE *psf, uint nAtoms, uint nBonds,
                    std::vector<std::vector<uint>> &moleculeXAtomIDY);
  ~BondAdjacencyList();
  adjNode *getAdjListNode(int value, int weight, adjNode *head);
  void display_AdjList(adjNode *ptr, int i);
  void connectedComponents(std::vector<std::vector<uint>> &moleculeXAtomIDY);
  void DFSUtil(int v, adjNode *node, adjNode **head, bool *visited,
               std::vector<uint> &moleculeX);

  graphEdge *edges;
};
#endif /*BOND_ADJACENCY_LIST_H*/