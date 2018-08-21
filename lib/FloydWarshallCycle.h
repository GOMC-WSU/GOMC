#pragma once
#include <vector>
#include <algorithm>
#include <assert.h>

class FloydWarshallCycle
{
public:
	// In constructor you have to mention the number of nodes
	FloydWarshallCycle(int numberOfNodes);
	~FloydWarshallCycle();

	// You can add more edges to the graph by calling this function
	// src and dest doesn't matter and you can add the nodes by any order i.e. (i, j) == (j, i)
	void AddEdge(int src, int dest);

	// Get the shortest path for specific node
	std::vector<int> GetShortestCycle(int src);

	std::vector< std::vector<int> > GetAllUniqueCycles();

	std::vector< std::vector<int> > GetAllUniqueCyclesAndCommonCycles();

private:
	// Number of nodes
	int numberOfNodes;

	// This vector will hold the edges as a 2D array. Each edge will be a 2 element vector -> (src, dest)
	std::vector< std::vector<int> > connections;

	// This vector will hold the weights of the edges.
	// We assume every connection has a weight of 1
	// So by running Floyd Warshall on this graph the shortest path
	//     will be equal to number of edges
	std::vector< std::vector<int> > graph;

	// This vector will hold the next jump for a shortest path
	// This vector is required to get the full path
	std::vector< std::vector<int> > next;

	// Run the Floyd-Warshall algorithm
	void floydWarshall();

	// Get the full path after running Floyd-Warshall algorithm
	// You need a source and destination
	// Or simply by giving the connections index
	std::vector<int> getPath(int src, int dest);
	std::vector<int> getPath(int connectionIndex);

	// Set everything to default values
	void setDefaults();

	// Set the values for every edges except index
	void setValues(int exceptThisOne);

	// Return all connections including an index
	std::vector<int> getConnectionsFor(int index);

	// Find the minimum of 2D array and return the shortest path
	std::vector<int> findMinimumPath(std::vector< std::vector<int> > paths);

	// Find the unique cycles 
	std::vector< std::vector<int> > getUniqueVectors(std::vector< std::vector<int> > allCycles);

	// Compare two vectors
	bool isVectorsEqual(std::vector<int> first, std::vector<int> second);

	// See if two vectors have common elements
	bool haveCommonElements(std::vector<int> first, std::vector<int> second);

	// Return the combined set
	std::vector<int> returnCombinedSet(std::vector<int> first, std::vector<int> second);
};

