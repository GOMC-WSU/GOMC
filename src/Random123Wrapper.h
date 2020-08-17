#pragma once

#include "Random123/philox.h"
typedef r123::Philox4x32 RNG;

class Random123Wrapper {
public:
  Random123Wrapper() {
    c = {{}};
    uk = {{}};
  }

  void SetStep(unsigned int step) { uk[0] = step; }
  void SetRandomSeed(unsigned int seedValue) { uk[1] = seedValue; }
  double GetRandomNumber(unsigned int counter) {
    c[0] = counter;
    RNG::key_type k = uk;
    RNG::ctr_type r = rng(c, k);
    double r01 = r[0];
    r01 /= UINT_MAX;
    return r01;
  }
  unsigned int GetStep() { return uk[0]; }
  unsigned int GetSeedValue() { return uk[1]; }
  double operator() (unsigned int counter) { return GetRandomNumber(counter); }

private:
  RNG::ctr_type c;
  RNG::key_type uk;
  RNG rng;
};