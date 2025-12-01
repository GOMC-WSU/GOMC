//===----------------------------------------------------------------------===//
//
// Written by Xing Mingjie (mingjie.xing@gmail.com).
//
// An implementation of the Johnson's circuit finding algorithm [1].
//
// [1] Donald B. Johnson, Finding all the elementary circuits of a directed
//     graph, SIAM Journal on Computing, 1975.
//
//===----------------------------------------------------------------------===//

#include "CircuitFinder.h"

void CircuitFinder::addEdge(int src, int dest) {
  AK[src].push_back(dest);
  E++;
}

void CircuitFinder::unblock(int U) {
  Blocked[U] = false;

  while (!B[U].empty()) {
    int W = B[U].front();
    B[U].pop_front();

    if (Blocked[W]) {
      unblock(W);
    }
  }
}

bool CircuitFinder::circuit(int V) {
  bool F = false;
  Stack.push_back(V);
  Blocked[V] = true;

  for (int W : AK[V]) {
    if (W == S) {
      SaveCycle();
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

void CircuitFinder::SaveCycle() {
  if (Stack.size() > 2) {

    std::vector<int> cycle;
    // std::cout << "circuit: ";
    for (auto I = Stack.begin(), E = Stack.end(); I != E; ++I) {
      // std::cout << *I << " -> ";
      cycle.push_back(*I);
    }
    uniqueCycles.push_back(cycle);
    // std::cout << *Stack.begin() << std::endl;
  }
}

void CircuitFinder::run() {
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

bool CircuitFinder::haveCommonElements(std::vector<int> first,
                                       std::vector<int> second) {
  for (int i = 0; i < (int)first.size(); i++) {
    if (std::find(second.begin(), second.end(), first[i]) != second.end()) {
      return true;
    }
  }
  return false;
}

std::vector<int> CircuitFinder::returnCombinedSet(std::vector<int> first,
                                                  std::vector<int> second) {
  std::vector<int> ret(first);
  for (int i = 0; i < (int)second.size(); i++) {
    if (std::find(ret.begin(), ret.end(), second[i]) == ret.end()) {
      ret.push_back(second[i]);
    }
  }
  return ret;
}

std::vector<std::vector<int>> CircuitFinder::GetAllCommonCycles() {
  run();
  std::vector<std::vector<int>> commons;
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