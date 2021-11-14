#pragma once

#include <climits>
#include "BasicTypes.h"
#include "Random123/philox.h"
typedef r123::Philox4x64 RNG;

class Random123Wrapper
{
public:
  Random123Wrapper();

  void SetStep(ulong step);
  void SetRandomSeed(ulong seedValue);
  void SetKey(unsigned int key);

  double GetRandomNumber(unsigned int counter);
  XYZ GetRandomCoords(unsigned int counter);

  unsigned int GetStep();
  unsigned int GetKeyValue();
  unsigned int GetSeedValue();
  double operator() (unsigned int counter);

  double GetGaussian(unsigned int counter);
  double GetGaussianNumber(unsigned int counter, double mean, double stdDev);
  XYZ GetGaussianCoords(unsigned int counter, double mean, double stdDev);

  double GetSymRandom(unsigned int counter, double bound);
  XYZ GetSymRandomCoords(unsigned int counter, double bound);
  XYZ GetRandomCoordsOnSphere(unsigned int counter);

private:
  RNG::ctr_type c;
  RNG::key_type uk;
  RNG rng;
};