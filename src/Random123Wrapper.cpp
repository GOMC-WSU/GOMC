#include "Random123Wrapper.h"
#include "Random123/boxmuller.hpp"

Random123Wrapper::Random123Wrapper()
{
  c = {{}};
  uk = {{}};
}

void Random123Wrapper::SetStep(ulong step)
{
  uk[0] = step;
}
void Random123Wrapper::SetRandomSeed(ulong seedValue)
{
  uk[1] = seedValue;
}

double Random123Wrapper::GetRandomNumber(unsigned int counter)
{
  c[0] = counter;
  RNG::key_type k = uk;
  RNG::ctr_type r = rng(c, k);
  double r01 = r[0] * RAND_INTERVAL;
  return r01;
}

XYZ Random123Wrapper::GetRandomCoords(unsigned int counter)
{
  c[0] = counter;
  RNG::key_type k = uk;
  RNG::ctr_type r = rng(c, k);
  XYZ r01;
  r01.x = static_cast<double>(r[0]) * RAND_INTERVAL;
  r01.y = static_cast<double>(r[1]) * RAND_INTERVAL;
  r01.z = static_cast<double>(r[2]) * RAND_INTERVAL;
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
  r123::double2 normal2 = r123::boxmuller(r[0], r[1]);
  return normal2.x;
}

double Random123Wrapper::GetGaussianNumber(unsigned int counter, double mean, double stdDev)
{
  double gNum = this->GetGaussian(counter);
  return (mean + gNum * stdDev);
}

XYZ Random123Wrapper::GetGaussianCoords(unsigned int counter, double mean, double stdDev)
{
  c[0] = counter;
  RNG::key_type k = uk;
  RNG::ctr_type r = rng(c, k);
  r123::double2 normal1 = r123::boxmuller(r[0], r[1]);
  r123::double2 normal2 = r123::boxmuller(r[2], r[3]);

  XYZ normals(mean + normal1.x * stdDev, mean + normal1.y * stdDev, mean + normal2.x * stdDev);
  return normals;
}