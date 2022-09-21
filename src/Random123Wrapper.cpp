#include "Random123Wrapper.h"

#include "Random123/boxmuller.hpp"
#include "Random123/uniform.hpp"

Random123Wrapper::Random123Wrapper() {
  c = {{}};
  uk = {{}};
}

inline RNG::ctr_type Random123Wrapper::getRNG(unsigned int counter) {
  // Need to use the localc variable to avoid OpenMP race conditions
  RNG::ctr_type localc = {{}};
  localc[0] = counter;
  localc[1] = GetKeyValue();
  return rng(localc, uk);
}

void Random123Wrapper::SetStep(ulong step) { uk[0] = step; }

void Random123Wrapper::SetRandomSeed(ulong seedValue) { uk[1] = seedValue; }

void Random123Wrapper::SetKey(unsigned int key) { c[1] = key; }

unsigned int Random123Wrapper::GetStep() const { return uk[0]; }

unsigned int Random123Wrapper::GetKeyValue() const { return c[1]; }

unsigned int Random123Wrapper::GetSeedValue() const { return uk[1]; }

double Random123Wrapper::operator()(unsigned int counter) {
  return GetRandomNumber(counter);
}

double Random123Wrapper::GetRandomNumber(unsigned int counter) {
  RNG::ctr_type r = getRNG(counter);
  double r01 = r123::u01<double>(r[0]);
  return r01;
}

XYZ Random123Wrapper::GetRandomCoords(unsigned int counter) {
  RNG::ctr_type r = getRNG(counter);
  XYZ r01;
  r01.x = r123::u01<double>(r[0]);
  r01.y = r123::u01<double>(r[1]);
  r01.z = r123::u01<double>(r[2]);
  return r01;
}

double Random123Wrapper::GetSymRandom(unsigned int counter, double bound) {
  RNG::ctr_type r = getRNG(counter);
  double r01;
  r01 = bound * r123::uneg11<double>(r[0]);
  return r01;
}

XYZ Random123Wrapper::GetSymRandomCoords(unsigned int counter, double bound) {
  RNG::ctr_type r = getRNG(counter);
  XYZ r01;
  r01.x = bound * r123::uneg11<double>(r[0]);
  r01.y = bound * r123::uneg11<double>(r[1]);
  r01.z = bound * r123::uneg11<double>(r[2]);
  return r01;
}

// Returns a uniformly random point on the unit sphere
XYZ Random123Wrapper::GetRandomCoordsOnSphere(unsigned int counter) {
  RNG::ctr_type r = getRNG(counter);
  double r01;
  // picking phi uniformly will cluster points at poles
  // pick u = cos(phi) uniformly instead
  // start from r[1] because I used r[0] in GetSymRandom when called in
  // multiparticle
  double u = r123::uneg11<double>(r[1]);
  // theta must be [0, 2pi) !
  double theta = 2.0 * M_PI * r123::u01<double>(r[2]);
  double rootTerm = sqrt(1.0 - u * u);
  return XYZ(rootTerm * cos(theta), rootTerm * sin(theta), u);
}

double Random123Wrapper::GetGaussian(unsigned int counter) {
  RNG::ctr_type r = getRNG(counter);
  r123::double2 normal2 = r123::boxmuller(r[0], r[1]);
  return normal2.x;
}

double Random123Wrapper::GetGaussianNumber(unsigned int counter, double mean,
                                           double stdDev) {
  double gNum = this->GetGaussian(counter);
  return (mean + gNum * stdDev);
}

XYZ Random123Wrapper::GetGaussianCoords(unsigned int counter, double mean,
                                        double stdDev) {
  RNG::ctr_type r = getRNG(counter);
  r123::double2 normal1 = r123::boxmuller(r[0], r[1]);
  r123::double2 normal2 = r123::boxmuller(r[2], r[3]);

  XYZ normals(mean + normal1.x * stdDev, mean + normal1.y * stdDev,
              mean + normal2.x * stdDev);
  return normals;
}