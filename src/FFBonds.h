/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FF_BONDS_H
#define FF_BONDS_H

#include "BasicTypes.h" //For "uint"
#include "FFSetup.h"    //For initialization data
#include "NumLib.h"     //For "Sq" function
#include "VectorLib.h"  //For transfer vect --> array function

class FFBonds {
public:
  FFBonds(void) : Kb(NULL), b0(NULL), fixed(NULL) {}
  ~FFBonds(void) {
    delete[] Kb;
    delete[] b0;
    delete[] fixed;
  }

  double Calc(const uint kind, const double dist) const {
    return (fixed[kind] ? 0.0 : Kb[kind] * num::Sq(dist - b0[kind]));
  }

  double Length(const uint kind) const { return b0[kind]; }

  void Init(ff_setup::Bond const &bond) {
    count = bond.getKbcnt();
    Kb = bond.CopyKb();
    b0 = bond.Copyb0();
    fixed = bond.Copyfixed();
  }

private:
  double *Kb, *b0;
  bool *fixed;
  uint count;
};

#endif /*FF_BONDS_H*/
