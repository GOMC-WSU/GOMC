/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_BONDS_H
#define FF_BONDS_H

#include "FFSetup.h" //For initialization data
#include "../lib/BasicTypes.h" //For "uint"
#include "../lib/VectorLib.h" //For transfer vect --> array function
#include "../lib/NumLib.h" //For "Sq" function

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

   double Calc(const uint kind, const double dist) const
   {
      return (fixed[kind] ? 0.0 : Kb[kind] * num::Sq(dist-b0[kind]));
   }

   double Length(const uint kind) const
   {
      return b0[kind];
   }

   void Init(ff_setup::Bond const& bond)
   {
      count = bond.Kb.size();
      Kb = vect::transfer(bond.Kb);
      b0 = vect::transfer(bond.b0);
      fixed = vect::transfer(bond.fixed);
   }
 private:
   double * Kb, * b0;
   bool * fixed;
   uint count;
};

#endif /*FF_BONDS_H*/

