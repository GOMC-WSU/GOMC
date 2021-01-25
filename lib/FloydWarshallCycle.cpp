/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "FloydWarshallCycle.h"

FloydWarshallCycle::FloydWarshallCycle(int numNodes)
{
  // Confirm that the number of nodes is greater than zero
  assert(numNodes > 0);

  // Store the size of the graph in numberOfNodes
  numberOfNodes = numNodes;

  // Resize both graph and next vectors to the size of N*N
  graph.resize(numberOfNodes);
  for (int i = 0; i < numberOfNodes; i++)
    graph[i].resize(numberOfNodes);
  next.resize(numberOfNodes);
  for (int i = 0; i < numberOfNodes; i++)
    next[i].resize(numberOfNodes);

  // double check
  assert(graph.size() == numberOfNodes);
}


FloydWarshallCycle::~FloydWarshallCycle()
{
}

void FloydWarshallCycle::AddEdge(int src, int dest)
{
  // make sure src and dest are valid node IDs
  assert(src >= 0);
  assert(src < numberOfNodes);
  assert(dest >= 0);
  assert(dest < numberOfNodes);

  std::vector<int> temp;
  temp.push_back(src);
  temp.push_back(dest);
  connections.push_back(temp);
}

std::vector<int> FloydWarshallCycle::GetShortestCycle(int src)
{
  // check src value
  assert(src >= 0);
  assert(src < numberOfNodes);

  std::vector< std::vector<int> > paths;
  std::vector<int> allEdgesIncludingSrc = getConnectionsFor(src);
  for (int i = 0; i < allEdgesIncludingSrc.size(); i++) {
    setDefaults();
    setValues(allEdgesIncludingSrc[i]);
    floydWarshall();
    paths.push_back(getPath(allEdgesIncludingSrc[i]));
  }
  return findMinimumPath(paths);
}

std::vector< std::vector<int> > FloydWarshallCycle::GetAllUniqueCycles()
{
  // check number of nodes
  assert(numberOfNodes >= 0);

  std::vector< std::vector<int> > allCycles;
  for (int i = 0; i < numberOfNodes; i++) {
    allCycles.push_back(GetShortestCycle(i));
  }
  return getUniqueVectors(allCycles);
}

std::vector< std::vector<int> > FloydWarshallCycle::GetAllCommonCycles()
{
  std::vector< std::vector<int> > uniqueCycles = GetAllUniqueCycles();
  std::vector< std::vector<int> > commons;
  int len = uniqueCycles.size();
  std::vector<bool> visited(len, false);

  for (int i = 0; i < len; i++) {
    if (!visited[i]) {
      std::vector<int> combined(uniqueCycles[i]);
      visited[i] = true;
      for (int j = i + 1; j < len; j++) {
        if (haveCommonElements(combined, uniqueCycles[j])) {
          combined = returnCombinedSet(combined, uniqueCycles[j]);
          visited[j] = true;
        }
      }
      commons.push_back(combined);
    }
  }

  return commons;
}

void FloydWarshallCycle::floydWarshall()
{
  for (int k = 0; k < numberOfNodes; k++) {
    for (int i = 0; i < numberOfNodes; i++) {
      for (int j = 0; j < numberOfNodes; j++) {
        if (graph[i][j] > graph[i][k] + graph[k][j]) {
          graph[i][j] = graph[i][k] + graph[k][j];
          next[i][j] = next[i][k];
        }
      }
    }
  }
}

std::vector<int> FloydWarshallCycle::getPath(int src, int dest)
{
  // make sure src and dest are valid node IDs
  assert(src >= 0);
  assert(src < numberOfNodes);
  assert(dest >= 0);
  assert(dest < numberOfNodes);

  std::vector<int> path;
  if (next[src][dest] == -1)
    return path;
  path.push_back(src);
  while (src != dest) {
    src = next[src][dest];
    path.push_back(src);
  }
  return path;
}

std::vector<int> FloydWarshallCycle::getPath(int connectionIndex)
{
  int src = connections[connectionIndex][0];
  int dest = connections[connectionIndex][1];
  std::vector<int> path;
  if (next[src][dest] == -1)
    return path;
  path.push_back(src);
  while (src != dest) {
    src = next[src][dest];
    path.push_back(src);
  }
  return path;
}

void FloydWarshallCycle::setDefaults()
{
  for (int i = 0; i < numberOfNodes; i++)
    for (int j = 0; j < numberOfNodes; j++) {
      graph[i][j] = 10000;
      next[i][j] = -1;
    }
}

void FloydWarshallCycle::setValues(int exceptThisOne)
{
  for (int j = 0; j < connections.size(); j++) {
    if (exceptThisOne != j) {
      int zero = connections[j][0];
      int one = connections[j][1];
      graph[zero][one] = 1;
      graph[one][zero] = 1;
      next[zero][one] = one;
      next[one][zero] = zero;
    }
  }
}

std::vector<int> FloydWarshallCycle::getConnectionsFor(int index)
{
  std::vector<int> conn;
  for (int i = 0; i < connections.size(); i++) {
    if (connections[i][0] == index || connections[i][1] == index)
      conn.push_back(i);
  }
  return conn;
}

std::vector<int> FloydWarshallCycle::findMinimumPath(const std::vector< std::vector<int> > &paths)
{
  //The path shouldn't be more than numberOfNodes without negative cycles
  int min = 2 * numberOfNodes;
  std::vector<int> shortest;
  for (int i = 0; i < paths.size(); i++) {
    if (paths[i].size() != 0 && paths[i].size() < min) {
      min = paths[i].size();
      shortest = paths[i];
    }
  }
  return shortest;
}

std::vector< std::vector<int> > FloydWarshallCycle::getUniqueVectors(std::vector< std::vector<int> > allCycles)
{
  int j;
  std::vector< std::vector<int> > uniqueCycles;
  for (int i = 0; i < allCycles.size(); i++) {
    std::sort(allCycles[i].begin(), allCycles[i].end());
  }
  for (int i = 0; i < allCycles.size(); i++) {
    if (uniqueCycles.size() == 0) {
      uniqueCycles.push_back(allCycles[i]);
    } else {
      for (j = 0; j < uniqueCycles.size(); j++)
        if (isVectorsEqual(uniqueCycles[j], allCycles[i]))
          break;
      if (j == uniqueCycles.size()) {
        uniqueCycles.push_back(allCycles[i]);
      }
    }
  }
  return uniqueCycles;
}

bool FloydWarshallCycle::isVectorsEqual(const std::vector<int> &first, const std::vector<int> &second)
{
  if (first.size() != second.size())
    return false;
  for (int i = 0; i < first.size(); i++) {
    if (first[i] != second[i]) {
      return false;
    }
  }
  return true;
}

bool FloydWarshallCycle::haveCommonElements(const std::vector<int> &first, const std::vector<int> &second)
{
  for (int i = 0; i < first.size(); i++) {
    if (std::find(second.begin(), second.end(), first[i]) != second.end()) {
      return true;
    }
  }
  return false;
}

std::vector<int> FloydWarshallCycle::returnCombinedSet(const std::vector<int> &first, const std::vector<int> &second)
{
  std::vector<int> ret(first);
  for (int i = 0; i < second.size(); i++) {
    if (std::find(ret.begin(), ret.end(), second[i]) == ret.end()) {
      ret.push_back(second[i]);
    }
  }
  return ret;
}

int FloydWarshallCycle::GetCentricNode()
{
  // Got it from:
  // https://codeforces.com/blog/entry/17974

  // values of eccentricity
  std::vector<int> ecc(numberOfNodes, 0);
  int i, j, raduse;
  setDefaults();
  
  for (j = 0; j < connections.size(); j++) {
    int zero = connections[j][0];
    int one = connections[j][1];
    graph[zero][one] = 1;
    graph[one][zero] = 1;
  }

  floydWarshall();

  // Counting values of eccentricity which is
  // maximal distance from that node to other
  for (i = 0; i < numberOfNodes; i++) {
    for (j = 0; j < numberOfNodes; j++) {
      if (ecc[i] < graph[i][j]) {
        ecc[i] = graph[i][j];
      }
    }
  }
  
  raduse = 1000000;  
  // Having values of eccentricity of all nodes, 
  // we can define radius of graph as minimal one among them
  for (i = 0; i < numberOfNodes; i++) {
    if (raduse > ecc[i]) {
      raduse = ecc[i];
    }
  }

  // Center of graph is set of nodes with eccentricity equal to the radius of graph:
  for (i = 0; i < numberOfNodes; i++) {
    if (raduse == ecc[i]) {
      return i;
    }
  }

  return -1;

}