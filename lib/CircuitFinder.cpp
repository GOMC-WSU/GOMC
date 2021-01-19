#include "CircuitFinder.h"

void CircuitFinder::addEdge(int src, int dest)
{
	AK[src].push_back(dest); 
	E++; 
}


void CircuitFinder::unblock(int U)
{
  Blocked[U] = false;

  while (!B[U].empty()) {
    int W = B[U].front();
    B[U].pop_front();

    if (Blocked[W]) {
      unblock(W);
    }
  }
}


bool CircuitFinder::circuit(int V)
{
  bool F = false;
  Stack.push_back(V);
  Blocked[V] = true;

  for (int W : AK[V]) {
    if (W == S) {
      output();
      F = true;
    } else if (W > S && !Blocked[W]) {
      F = circuit(W);
    }
  }

  if (F) {
    unblock(V);
  } else {
    for (int W : AK[V]) {
      auto IT = std::find(B[W].begin(), B[W].end(), V);
      if (IT == B[W].end()) {
        B[W].push_back(V);
      }
    }
  }

  Stack.pop_back();
  return F;
}


void CircuitFinder::output()
{
 if (Stack.size() > 2){
    std::cout << "circuit: ";
    for (auto I = Stack.begin(), E = Stack.end(); I != E; ++I) {
        std::cout << *I << " -> ";
    }
    std::cout << *Stack.begin() << std::endl;
 }
}


void CircuitFinder::run()
{
  Stack.clear();
  S = 0;

  while (S < V) {
    for (int I = S; I < V; ++I) {
      Blocked[I] = false;
      B[I].clear();
    }
    circuit(S);
    ++S;
  }
}