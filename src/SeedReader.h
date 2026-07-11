/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef SEED_READER_H
#define SEED_READER_H

#include "Reader.h"

struct SeedReader : Reader {
  inline SeedReader(std::string const &nm, std::string const &als,
                    const bool crit, const bool note)
      : Reader(nm, als, crit, note) {}

  // Keep going until desired var hit.
  inline Reader &GotoStep(const ulong step) {
    while (!file.eof() && !RdLn().StrtComp()("STEP")(step).FinComp()) {
    }
    if (file.eof()) {
      printf("ERROR: Step #%u seed value not found in file %s,\n"
             "while loading system to restart.\n",
             step, fileNm.c_str());
      exit(1);
    }
    return *this;
  }
};

#endif /*SEED_READER_H*/
