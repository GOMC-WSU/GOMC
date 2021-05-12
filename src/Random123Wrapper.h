#pragma once

#include "Random123/philox.h"
typedef r123::Philox4x64 RNG;

class Random123Wrapper
{
public:
  Random123Wrapper();

  void SetStep(ulong step);
  void SetRandomSeed(ulong seedValue);

  double GetRandomNumber(unsigned int counter);

  unsigned int GetStep();
  unsigned int GetSeedValue();
  double operator() (unsigned int counter);

  double GetGaussian(unsigned int counter);
  double GetGaussianNumber(unsigned int counter, double mean, double stdDev);

private:
  RNG::ctr_type c;
  RNG::key_type uk;
  RNG rng;
};