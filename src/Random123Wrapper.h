/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef RANDOM123_WRAPPER_H
#define RANDOM123_WRAPPER_H

#include "BasicTypes.h"
#include "Random123/philox.h"
typedef r123::Philox4x64 RNG;

class Random123Wrapper {
public:
  Random123Wrapper();

  void SetStep(ulong step);
  void SetRandomSeed(ulong seedValue);
  void SetKey(unsigned int key);

  unsigned int GetStep() const;
  unsigned int GetKeyValue() const;
  unsigned int GetSeedValue() const;
  double operator()(unsigned int counter);

  double GetRandomNumber(unsigned int counter);
  XYZ GetRandomCoords(unsigned int counter);

  double GetSymRandom(unsigned int counter, double bound);
  XYZ GetSymRandomCoords(unsigned int counter, double bound);
  XYZ GetRandomCoordsOnSphere(unsigned int counter);

  double GetGaussian(unsigned int counter);
  double GetGaussianNumber(unsigned int counter, double mean, double stdDev);
  XYZ GetGaussianCoords(unsigned int counter, double mean, double stdDev);

private:
  inline RNG::ctr_type getRNG(unsigned int counter);

  RNG::ctr_type c;
  RNG::key_type uk;
  RNG rng;
};

#endif /*RANDOM123_WRAPPER_H*/
