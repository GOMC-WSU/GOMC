/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SEED_READER_H
#define SEED_READER_H

#include "Reader.h"

struct SeedReader : Reader {

  inline SeedReader(std::string const& nm, std::string const& als,
                    const bool crit, const bool note) :
    Reader(nm, als, crit, note)
  {
  }

  //Keep going until desired var hit.
  inline Reader & GotoStep(const unsigned int step)
  {
    while (!file.eof() && !RdLn().StrtComp()("STEP")(step).FinComp()) {}
    if (file.eof()) {
      printf("ERROR: Step #%u seed value not found in file %s,\n"
             "while loading system to restart.\n", step, fileNm.c_str());
      exit(1);
    }
    return *this;
  }
};

#endif /*SEED_READER_H*/
