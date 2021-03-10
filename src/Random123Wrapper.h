#pragma once

#include "Random123/philox.h"
typedef r123::Philox4x32 RNG;

class Random123Wrapper
{
public:
  Random123Wrapper();

  void SetStep(unsigned int step);
  void SetRandomSeed(unsigned int seedValue);

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