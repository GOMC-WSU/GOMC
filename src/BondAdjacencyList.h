/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

/* Courtesy of https://www.softwaretestinghelp.com/graph-implementation-cpp/ */

#ifndef BONDADJACENCYLIST_H
#define BONDADJACENCYLIST_H

#include <vector>
#include <cstdio>
#include "BasicTypes.h"
#include <iostream>
#include <limits.h>
#include <algorithm>


// stores adjacency list items
struct adjNode {
    int val;
    int cost;
    adjNode* next;
};
// structure to store edges
struct graphEdge {
    int start_ver, end_ver, weight;
};

class BondAdjacencyList{
public:

adjNode **head;                //adjacency list as array of pointers
uint nAtoms;  // number of nodes in the graph
uint nBonds; // number of edges in the graph


BondAdjacencyList(FILE* psf, uint nAtoms, uint nBonds, std::vector< std::vector<int> > & moleculeXAtomIDY);
~BondAdjacencyList();
adjNode* getAdjListNode(int value, int weight, adjNode* head);
void display_AdjList(adjNode* ptr, int i);
void connectedComponents(std::vector< std::vector<int> > & moleculeXAtomIDY);
void DFSUtil(int v, adjNode * node, adjNode ** head, bool * visited, std::vector<int> & moleculeX);

graphEdge * edges;

};
#endif