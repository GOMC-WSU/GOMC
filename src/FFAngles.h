/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_ANGLES_H
#define FF_ANGLES_H

#include "FFSetup.h" //For initialization data
#include "BasicTypes.h" //For "uint"
#include "VectorLib.h" //For transfer vect --> array function
#include "NumLib.h" //For "Sq" function
#include <math.h>

class FFAngles
{
public:
  FFAngles(void) : Ktheta(NULL), theta0(NULL), fixed(NULL) {}
  ~FFAngles(void)
  {
    delete[] Ktheta;
    delete[] theta0;
    delete[] fixed;
  }

  real Angle(const uint kind) const
  {
    return theta0[kind];
  }

  real AngleEnergy(const uint kind) const
  {
    return Ktheta[kind];
  }

  bool AngleFixed(const uint kind) const
  {
    return fixed[kind];
  }

  virtual real Calc(const uint kind, const real ang) const
  {
    return (fixed[kind] ? 0.0 : Ktheta[kind] * num::Sq(ang - theta0[kind]));
  }


  void Init(ff_setup::Angle const& angle)
  {
    count = angle.getKthetacnt();
    Ktheta = angle.CopyKtheta();
    theta0 = angle.Copytheta0();
    fixed = angle.Copyfixed();
  }

protected:
  real * Ktheta, * theta0;
  bool * fixed;
  uint count;
};

class FFAngleMartini : public FFAngles
{
  virtual real Calc(const uint kind, const real ang) const
  {
    return (fixed[kind] ? 0.0 : Ktheta[kind] *
            num::Sq(cos(ang) - cos(theta0[kind])));
  }

};
#endif /*FF_ANGLES_H*/
