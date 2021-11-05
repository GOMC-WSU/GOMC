#pragma once

#include <climits>
#include "BasicTypes.h"
#include "Random123/philox.h"
typedef r123::Philox4x64 RNG;
static const double RAND_INTERVAL = 1.0/static_cast<double>(ULONG_MAX);

class Random123Wrapper
{
public:
  Random123Wrapper();

  void SetStep(ulong step);
  void SetRandomSeed(ulong seedValue);

  double GetRandomNumber(unsigned int counter);
  XYZ GetRandomCoords(unsigned int counter);

  unsigned int GetStep();
  unsigned int GetSeedValue();
  double operator() (unsigned int counter);

  double GetGaussian(unsigned int counter);
  double GetGaussianNumber(unsigned int counter, double mean, double stdDev);
  XYZ GetGaussianCoords(unsigned int counter, double mean, double stdDev);

private:
  RNG::ctr_type c;
  RNG::key_type uk;
  RNG rng;
};