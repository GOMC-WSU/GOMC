#include "Random123Wrapper.h"
#include "Random123/boxmuller.hpp"
#include <climits>

Random123Wrapper::Random123Wrapper()
{
  c = {{}};
  uk = {{}};
}

void Random123Wrapper::SetStep(unsigned int step)
{
  uk[0] = step;
}
void Random123Wrapper::SetRandomSeed(unsigned int seedValue)
{
  uk[1] = seedValue;
}

double Random123Wrapper::GetRandomNumber(unsigned int counter)
{
  c[0] = counter;
  RNG::key_type k = uk;
  RNG::ctr_type r = rng(c, k);
  double r01 = r[0];
  r01 /= UINT_MAX;
  return r01;
}

unsigned int Random123Wrapper::GetStep()
{
  return uk[0];
}

unsigned int Random123Wrapper::GetSeedValue()
{
  return uk[1];
}

double Random123Wrapper::operator() (unsigned int counter)
{
  return GetRandomNumber(counter);
}

double Random123Wrapper::GetGaussian(unsigned int counter)
{
  c[0] = counter;
  RNG::key_type k = uk;
  RNG::ctr_type r = rng(c, k);
  r123::float2 normalf2 = r123::boxmuller(r[0], r[1]);
  return double(normalf2.x);
}

double Random123Wrapper::GetGaussianNumber(unsigned int counter, double mean, double stdDev)
{
  double gNum = this->GetGaussian(counter);
  return (mean + gNum * stdDev);
}
