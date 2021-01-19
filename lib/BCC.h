#include <iostream> 
#include <list> 
#include <stack> 

class Edge { 
public: 
	int u; 
	int v; 
	Edge(int u, int v); 
}; 

// A class that represents an undirected graph 
class Graph { 
	int V; // No. of vertices 
	int E; // No. of edges 
	std::list<int>* adj; // A dynamic array of adjacency lists 

	// A Recursive DFS based function used by BCC() 
	void BCCUtil(int u, int disc[], int low[], 
				std::list<Edge>* st, int parent[]); 

public: 
	Graph(int V); // Constructor 
	void addEdge(int v, int w); // function to add an edge to graph 
	void BCC(); // prints strongly connected components 
}; 