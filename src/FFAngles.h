/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef FF_ANGLES_H
#define FF_ANGLES_H

#include "BasicTypes.h" //For "uint"
#include "FFSetup.h"    //For initialization data
#include "NumLib.h"     //For "Sq" function
#include "VectorLib.h"  //For transfer vect --> array function

class FFAngles {
public:
  FFAngles(void) : Ktheta(NULL), theta0(NULL), fixed(NULL) {}
  virtual ~FFAngles(void) {
    delete[] Ktheta;
    delete[] theta0;
    delete[] fixed;
  }

  double Angle(const uint kind) const { return theta0[kind]; }

  double AngleEnergy(const uint kind) const { return Ktheta[kind]; }

  bool AngleFixed(const uint kind) const { return fixed[kind]; }

  virtual double Calc(const uint kind, const double ang) const {
    return (fixed[kind] ? 0.0 : Ktheta[kind] * num::Sq(ang - theta0[kind]));
  }

  void Init(ff_setup::Angle const &angle) {
    count = angle.getKthetacnt();
    Ktheta = angle.CopyKtheta();
    theta0 = angle.Copytheta0();
    fixed = angle.Copyfixed();
  }

protected:
  double *Ktheta, *theta0;
  bool *fixed;
  uint count;
};

class FFAngleMartini : public FFAngles {
  virtual double Calc(const uint kind, const double ang) const {
    return (fixed[kind] ? 0.0
                        : Ktheta[kind] * num::Sq(cos(ang) - cos(theta0[kind])));
  }
};
#endif /*FF_ANGLES_H*/
