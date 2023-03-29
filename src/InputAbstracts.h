/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef INPUT_ABSTRACTS_H
#define INPUT_ABSTRACTS_H

#include <string>
#include <vector>

#include "BasicTypes.h" //For ulong
#include "BitLib.h"     //For bit check, mask generating functions
#include "FFConst.h"
#include "FixedWidthReader.h"
#include "Reader.h"
#include "StrLib.h" //For string comparison wrapper.

////////////////////////////
///  READ
//

struct ReadableBase {
  virtual void Read(Reader &file) = 0;
  virtual ~ReadableBase() {}
};

struct FWReadableBase {
  virtual void Read(FixedWidthReader &file) = 0;
  virtual ~FWReadableBase() {}
};
struct ReadableBaseWithFirst {
  virtual void Read(Reader &file, std::string const &firstVar) = 0;
  virtual ~ReadableBaseWithFirst() {}
};

// N is # of elements involved in search
struct ReadableStepDependentBase {
  virtual void Read(Reader &file, const ulong totalSteps,
                    std::string const &kindName) = 0;
  virtual ~ReadableStepDependentBase() {}
};

// template <uint N>
class SearchableBase {
public:
  SearchableBase(uint lN) : N(lN), mask(bits::GetMasks(lN)) {}
  // ubiquitous find for element identified by N maskable strings.
  int Find(std::string const *const kindNames,
           std::vector<std::string> const &readGeomMerged,
           std::string const &wild = ff::part::WILD) const {
    int found = NOT_FOUND;
    // NOTE: should be < N... was wrong before, wouldn't work for
    // particle case...
    for (uint w = 0; w < N && found == NOT_FOUND; w++)
      for (uint i = 0; i < readGeomMerged.size() && found == NOT_FOUND; i++)
        if (Match(kindNames, readGeomMerged[i], w, wild))
          found = i;
    return found;
  }

private:
  // Check a particular previously read item for a match, if wild
  // cards are active, loop through all of them.
  bool Match(std::string const *const kindNames, std::string const &merged,
             const unsigned numberOfWCBits, std::string const &wild) const {
    bool result = false;
    for (uint m = 0; m < mask[numberOfWCBits].size(); m++)
      result |= Match(kindNames, mask[numberOfWCBits][m], merged, wild);
    return result;
  }

  // Match one particular wild card or wildcardless set and its flip
  bool Match(std::string const *const kindNames, const uint wildCardMask,
             std::string const &readGeomMerged, std::string const &wild) const {
    return str::compare(SearchStr(kindNames, wildCardMask, false, wild),
                        readGeomMerged) ||
           str::compare(SearchStr(kindNames, wildCardMask, true, wild),
                        readGeomMerged);
  }

  // Gets string to use in search
  std::string SearchStr(std::string const *const v, const uint wildCardMask,
                        const bool flip, std::string const &wild) const {
    std::string s = "";
    for (uint i = 0; i < N; i++) {
      uint idx = (flip ? N - 1 - i : i);
      s += (bits::Check(wildCardMask, idx) ? wild : v[idx]);
    }
    return s;
  }
  uint N;
  const std::vector<std::vector<uint>> mask;
  static const int NOT_FOUND = -1;
};

#endif /*INPUT_ABSTRACTS_H*/
