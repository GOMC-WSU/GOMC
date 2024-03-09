/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef SUBDIV_ARRAY
#define SUBDIV_ARRAY

#include <cstddef>

#include "BasicTypes.h"

// Common class used for dihedral, sorted kind array, topology arrays, etc.

class SubdividedArray {
public:
  SubdividedArray() : start(NULL), subdivCount(0) {}
  SubdividedArray(SubdividedArray const &other) {
    subdivCount = other.subdivCount;
    start = new uint[other.subdivCount + 1];
    for (uint i = 0; i <= subdivCount; ++i)
      start[i] = other.start[i];
  }
  void Init(const uint subdiv) {
    if (start != NULL)
      Cleanup();
    subdivCount = subdiv;
    start = new uint[subdiv + 1];
  }
  void Set(const uint div, const uint first, const uint len) {
    start[div] = first;
    start[div + 1] = first + len;
  }

  ~SubdividedArray(void) { Cleanup(); }
  void Cleanup(void) {
    delete[] start;
    start = NULL;
  }

  // returns index of the offset-th element of kind
  uint Index(const uint kind, const uint offset) const {
    return start[kind] + offset;
  }

  // Return first el, last el, or length of subdiv
  uint Begin(const uint kind) const { return start[kind]; }
  uint End(const uint kind) const { return start[kind + 1]; }
  uint Length(const uint kind) const { return End(kind) - Begin(kind); }

  SubdividedArray &operator=(SubdividedArray other) {
    subdivCount = other.subdivCount;
    unsigned int *tmp = other.start;
    other.start = start;
    start = tmp;
    return *this;
  }

private:
  // note start is one longer than subdivCount to account for length of last
  // actual element
  uint *start, subdivCount;
};

#endif /*SUBDIV_ARRAY*/
