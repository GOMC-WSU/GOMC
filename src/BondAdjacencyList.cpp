/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

/* Courtesy of https://www.softwaretestinghelp.com/graph-implementation-cpp/ */

#include "BondAdjacencyList.h"


//return if we fail to read anything
const int READERROR = -1;

// insert new nodes into adjacency list from given graph
adjNode* BondAdjacencyList::getAdjListNode(int value, int weight, adjNode* head)   {
    adjNode* newNode = new adjNode;
    newNode->val = value;
    newNode->cost = weight;
    newNode->next = head;   // point new node to current head
    return newNode;
}

BondAdjacencyList::BondAdjacencyList(FILE* psf, uint nAtoms, uint nBonds, std::vector< std::vector<uint> > & moleculeXAtomIDY){

    this->nAtoms = nAtoms;
    this->nBonds = nBonds;

    edges = new graphEdge[nBonds];

    unsigned int atom0, atom1;
    int dummy;
    // Loads all the bonds into edges array
    for (uint i = 0; i < this->nBonds; i++) {
    dummy = fscanf(psf, "%u %u", &atom0, &atom1);
    if(dummy != 2) {
        fprintf(stderr, "ERROR: Incorrect Number of bonds in PSF file ");
    } else if (feof(psf) || ferror(psf)) {
        fprintf(stderr, "ERROR: Could not find all bonds in PSF file ");
    }
    edges[i].weight = 1;
    edges[i].start_ver = atom0-1;
    edges[i].end_ver = atom1-1;
    }

    // allocate new node
    head = new adjNode*[this->nAtoms]();
    // initialize head pointer for all vertices
    for (uint i = 0; i < this->nAtoms; i++)
        head[i] = nullptr;
    // construct directed graph by adding edges to it
    for (uint i = 0; i < this->nBonds; i++)  {
        int start_ver = edges[i].start_ver;
        int end_ver = edges[i].end_ver;
        int weight = edges[i].weight;
        // insert in the beginning
        adjNode* newNode = getAdjListNode(end_ver, weight, head[start_ver]);
           
        // point head pointer to new node
        head[start_ver] = newNode;
 
        start_ver = edges[i].end_ver;
        end_ver = edges[i].start_ver;
        weight = edges[i].weight;
        // insert in the beginning
        newNode = getAdjListNode(end_ver, weight, head[start_ver]);
            
        // point head pointer to new node
        head[start_ver] = newNode;
    
    }

    connectedComponents(moleculeXAtomIDY);
    /* For debugging 
    for (uint i = 0; i < this->nAtoms; i++)
        display_AdjList(head[i], i);
    
    std::cout << "Before sorting" << std::endl;
    for (std::vector< std::vector<uint> >::iterator it = moleculeXAtomIDY.begin();
        it != moleculeXAtomIDY.end(); it++){
        for (std::vector<uint>::iterator it2 = it->begin();
            it2 != it->end(); it2++){
            std::cout << *it2 << ", ";
        }
        std::cout << std::endl;
    }
    */
    /* Sort Atom Indices in N connected components, then sort N connected components by first atom index*/
    for (std::vector< std::vector<uint> >::iterator it = moleculeXAtomIDY.begin();
        it != moleculeXAtomIDY.end(); it++){
        std::sort(it->begin(), it->end());
    }
    std::sort(moleculeXAtomIDY.begin(), moleculeXAtomIDY.end());
    /* For debugging 
    std::cout << "After sorting" << std::endl;
    for (std::vector< std::vector<uint> >::iterator it = moleculeXAtomIDY.begin();
        it != moleculeXAtomIDY.end(); it++){
        for (std::vector<uint>::iterator it2 = it->begin();
            it2 != it->end(); it2++){
            std::cout << *it2 << ", ";
        }
        std::cout << std::endl;
    }
    */
}

// Destructor
BondAdjacencyList::~BondAdjacencyList() {
    for (int i = 0; i < this->nAtoms; i++){
        delete[] head[i];
    }
    delete[] edges;        
    delete[] head;
}

// print all adjacent vertices of given vertex
void BondAdjacencyList::display_AdjList(adjNode* ptr, int i)
{
    while (ptr != nullptr) {
        std::cout << "(" << i << ", " << ptr->val
            << ", " << ptr->cost << ") ";
        ptr = ptr->next;
    }
    std::cout << std::endl;
}

// Method to print connected components in an 
// undirected graph 
void BondAdjacencyList::connectedComponents(std::vector< std::vector<uint> > & moleculeXAtomIDY) 
{ 
    // Mark all the vertices as not visited 
    bool *visited = new bool[this->nAtoms]; 
    for(int v = 0; v < this->nAtoms; v++) 
        visited[v] = false; 
  
    for (int v=0; v<this->nAtoms; v++) 
    { 
        if (visited[v] == false) 
        { 
            // print all reachable vertices 
            // from v 
            /* For debugging
            std::cout << "Calling DFSUtil" << std::endl; */
            std::vector<uint> moleculeX;
            DFSUtil(v, this->head[v], this->head, visited, moleculeX); 
            moleculeXAtomIDY.push_back(moleculeX);
            /* For debugging std::cout << "\n"; */
        } 
    } 
    delete[] visited; 
} 
  
void BondAdjacencyList::DFSUtil(int v, adjNode * node, adjNode ** head, bool * visited, std::vector<uint> & moleculeX) 
{ 
    // Mark the current node as visited and print it 
    visited[v] = true; 
    /* For debugging std::cout << v << " "; */
    moleculeX.push_back(v);
    // Recur for all the vertices 
    // adjacent to this vertex
    while (node != nullptr){
        // outgoing edge : node->val
        v = node->val;
        if (visited[v]==false){
            visited[v] = true; 
            // Evaluate adjacency list of outgoing edge for prev visited
            DFSUtil(v, head[v], head, visited, moleculeX);
        }
        // Go to next node original node's adjacency list
        node = node->next;
    }
}
    