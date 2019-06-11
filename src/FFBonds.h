/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_BONDS_H
#define FF_BONDS_H

#include "FFSetup.h" //For initialization data
#include "BasicTypes.h" //For "uint"
#include "VectorLib.h" //For transfer vect --> array function
#include "NumLib.h" //For "Sq" function

class FFBonds
{
public:
  FFBonds(void) : Kb(NULL), b0(NULL), fixed(NULL) {}
  ~FFBonds(void)
  {
    delete[] Kb;
    delete[] b0;
    delete[] fixed;
  }

  real Calc(const uint kind, const real dist) const
  {
    return (fixed[kind] ? 0.0 : Kb[kind] * num::Sq(dist - b0[kind]));
  }

  real Length(const uint kind) const
  {
    return b0[kind];
  }

  void Init(ff_setup::Bond const& bond)
  {
    count = bond.getKbcnt();
    Kb = bond.CopyKb();
    b0 = bond.Copyb0();
    fixed = bond.Copyfixed();
  }
private:
  real * Kb, * b0;
  bool * fixed;
  uint count;
};

#endif /*FF_BONDS_H*/
