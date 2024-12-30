/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "FFParticle.h"

#include "NumLib.h" //For Sq, Cb, and MeanA/G functions.
#ifdef GOMC_CUDA
#include "ConstantDefinitionsCUDAKernel.cuh"
#endif

FFParticle::FFParticle(Forcefield &ff)
    : forcefield(ff), mass(NULL), nameFirst(NULL), nameSec(NULL), n(NULL),
      n_1_4(NULL), sigmaSq(NULL), sigmaSq_1_4(NULL), epsilon(NULL),
      epsilon_1_4(NULL), epsilon_cn(NULL), epsilon_cn_1_4(NULL),
      epsilon_cn_6(NULL), epsilon_cn_6_1_4(NULL), nOver6(NULL), nOver6_1_4(NULL)
#ifdef GOMC_CUDA
      ,
      varCUDA(NULL)
#endif
{
  exp6 = ff.exp6;
}

FFParticle::~FFParticle(void) {
  delete[] mass;
  delete[] nameFirst;
  delete[] nameSec;

  delete[] sigmaSq;
  delete[] n;
  delete[] epsilon;
  delete[] epsilon_cn;
  delete[] epsilon_cn_6;
  delete[] nOver6;
  // parameter for 1-4 interaction,
  // same one will be used for 1-3 interaction
  delete[] sigmaSq_1_4;
  delete[] n_1_4;
  delete[] epsilon_1_4;
  delete[] epsilon_cn_1_4;
  delete[] epsilon_cn_6_1_4;
  delete[] nOver6_1_4;

#ifdef GOMC_CUDA
  DestroyCUDAVars(varCUDA);
  delete varCUDA;
#endif
}

void FFParticle::Init(ff_setup::Particle const &mie,
                      ff_setup::NBfix const &nbfix) {
#ifdef GOMC_CUDA
  // Variables for GPU stored in here
  varCUDA = new VariablesCUDA();
#endif
  count = mie.epsilon.size(); // Get # particles read
  // Size LJ particle kind arrays
  mass = new double[count];
  // Size LJ-LJ pair arrays
  uint size = num::Sq(count);
  nameFirst = new std::string[size];
  nameSec = new std::string[size];

  n = new double[size];
  n_1_4 = new double[size];
  epsilon = new double[size];
  epsilon_cn = new double[size];
  epsilon_cn_6 = new double[size];
  nOver6 = new double[size];
  sigmaSq = new double[size];

  epsilon_1_4 = new double[size];
  epsilon_cn_1_4 = new double[size];
  epsilon_cn_6_1_4 = new double[size];
  nOver6_1_4 = new double[size];
  sigmaSq_1_4 = new double[size];

  // Combining VDW parameter
  Blend(mie);
  // Adjusting VDW parameter using NBFIX
  AdjNBfix(nbfix);

#ifdef GOMC_CUDA
  double diElectric_1 = 1.0 / forcefield.dielectric;
  InitGPUForceField(*varCUDA, sigmaSq, epsilon_cn, n, forcefield.vdwKind,
                    forcefield.isMartini, count, forcefield.rCut,
                    forcefield.rCutSq, forcefield.rCutCoulomb,
                    forcefield.rCutCoulombSq, forcefield.rCutLow,
                    forcefield.rswitch, forcefield.alpha, forcefield.alphaSq,
                    forcefield.ewald, diElectric_1);
#endif
}

double FFParticle::EnergyLRC(const uint kind1, const uint kind2) const {
  uint idx = FlatIndex(kind1, kind2);
  double tailCorrection = 1.0;
  double sigma = sqrt(sigmaSq[idx]);
  double rRat = sigma / forcefield.rCut;
  double N = (double)(n[idx]);
  tailCorrection *= sigma * sigmaSq[idx];
  tailCorrection *= 2.0 * M_PI * epsilon_cn[idx] / (N - 3.0);
  tailCorrection *=
      (pow(rRat, n[idx] - 3.0) - ((N - 3.0) / 3.0) * rRat * rRat * rRat);
  /*
    if(forcefield.freeEnergy) {
      //For free energy calc, we only consider attraction part
      tailCorrection *= (-((N - 3.0) / 3.0) * pow(rRat, 3));
    } else {
      tailCorrection *= (pow(rRat, n[idx] - 3) - ((N - 3.0) / 3.0) * pow(rRat,
    3));
    }
  */
  return tailCorrection;
}

double FFParticle::VirialLRC(const uint kind1, const uint kind2) const {
  uint idx = FlatIndex(kind1, kind2);
  double tailCorrection = 1.0;
  double sigma = sqrt(sigmaSq[idx]);
  double rRat = sigma / forcefield.rCut;
  double N = (double)(n[idx]);
  tailCorrection *= sigma * sigmaSq[idx];
  tailCorrection *= 2.0 * M_PI * epsilon_cn[idx] / (N - 3.0);
  tailCorrection *=
      (N * pow(rRat, n[idx] - 3.0) - (N - 3.0) * 2.0 * rRat * rRat * rRat);
  /*
    if(forcefield.freeEnergy) {
      //For free energy calc, we only consider attraction part
      tailCorrection *= (- (N - 3.0) * 2.0 * pow(rRat, 3));
    } else {
      tailCorrection *= (N * pow(rRat, n[idx] - 3) - (N - 3.0) * 2.0 * pow(rRat,
    3));
    }
  */
  return tailCorrection;
}

double FFParticle::ImpulsePressureCorrection(const uint kind1,
                                             const uint kind2) const {
  uint idx = FlatIndex(kind1, kind2);
  double tailCorrection = 1.0;
  double sigma = sqrt(sigmaSq[idx]);
  double rRat = sigma / forcefield.rCut;
  double N = (double)(n[idx]);
  tailCorrection *=
      2.0 * M_PI * epsilon_cn[idx] * forcefield.rCutSq * forcefield.rCut;
  tailCorrection *= (pow(rRat, N) - pow(rRat, 6.0));
  return tailCorrection;
}

void FFParticle::Blend(ff_setup::Particle const &mie) {
  for (uint i = 0; i < count; ++i) {
    for (uint j = 0; j < count; ++j) {
      uint idx = FlatIndex(i, j);
      // get all name combination for using in nbfix
      nameFirst[idx] = mie.getname(i);
      nameFirst[idx] += mie.getname(j);
      nameSec[idx] = mie.getname(j);
      nameSec[idx] += mie.getname(i);

      if (exp6) {
        n[idx] = num::MeanG(mie.n, mie.n, i, j);
        n_1_4[idx] = num::MeanG(mie.n_1_4, mie.n_1_4, i, j);
      } else {
        n[idx] = num::MeanA(mie.n, mie.n, i, j);
        n_1_4[idx] = num::MeanA(mie.n_1_4, mie.n_1_4, i, j);
      }
      // Warning: Causes division by zero if n (or n_1_4) = 0 or 6.
      // Handled through error checking in FFSetup::Read() to require
      // n (and n_1_4) > 6.
      double cn =
          n[idx] / (n[idx] - 6.0) * pow(n[idx] / 6.0, (6.0 / (n[idx] - 6.0)));
      double cn_1_4 = n_1_4[idx] / (n_1_4[idx] - 6.0) *
                      pow(n_1_4[idx] / 6.0, (6.0 / (n_1_4[idx] - 6.0)));

      double sigma, sigma_1_4;
      sigma = sigma_1_4 = 0.0;

      if (forcefield.vdwGeometricSigma) {
        sigma = num::MeanG(mie.sigma, mie.sigma, i, j);
        sigma_1_4 = num::MeanG(mie.sigma_1_4, mie.sigma_1_4, i, j);
      } else {
        sigma = num::MeanA(mie.sigma, mie.sigma, i, j);
        sigma_1_4 = num::MeanA(mie.sigma_1_4, mie.sigma_1_4, i, j);
      }

      sigmaSq_1_4[idx] = sigma_1_4 * sigma_1_4;
      sigmaSq[idx] = sigma * sigma;
      epsilon[idx] = num::MeanG(mie.epsilon, mie.epsilon, i, j);
      epsilon_cn[idx] = cn * epsilon[idx];
      epsilon_1_4[idx] = num::MeanG(mie.epsilon_1_4, mie.epsilon_1_4, i, j);
      epsilon_cn_1_4[idx] = cn_1_4 * epsilon_1_4[idx];
      epsilon_cn_6[idx] = epsilon_cn[idx] * 6;
      epsilon_cn_6_1_4[idx] = epsilon_cn_1_4[idx] * 6;
      nOver6[idx] = n[idx] / 6;
      nOver6_1_4[idx] = n_1_4[idx] / 6;
    }
  }
}

void FFParticle::AdjNBfix(ff_setup::NBfix const &nbfix) {
  uint size = num::Sq(count);
  for (uint i = 0; i < nbfix.epsilon.size(); i++) {
    for (uint j = 0; j < size; j++) {
      if (nbfix.getname(i) == nameFirst[j] || nbfix.getname(i) == nameSec[j]) {
        n[j] = nbfix.n[i];
        n_1_4[j] = nbfix.n_1_4[i];
        sigmaSq[j] = nbfix.sigma[i] * nbfix.sigma[i];
        sigmaSq_1_4[j] = nbfix.sigma_1_4[i] * nbfix.sigma_1_4[i];
        double cn = n[j] / (n[j] - 6) * pow(n[j] / 6, (6 / (n[j] - 6)));
        double cn_1_4 =
            n_1_4[j] / (n_1_4[j] - 6) * pow(n_1_4[j] / 6, (6 / (n_1_4[j] - 6)));

        epsilon[j] = nbfix.epsilon[i];
        epsilon_cn[j] = cn * nbfix.epsilon[i];
        epsilon_1_4[j] = nbfix.epsilon_1_4[i];
        epsilon_cn_1_4[j] = cn_1_4 * nbfix.epsilon_1_4[i];
        epsilon_cn_6[j] = epsilon_cn[j] * 6;
        epsilon_cn_6_1_4[j] = epsilon_cn_1_4[j] * 6;
        nOver6[j] = n[j] / 6;
        nOver6_1_4[j] = n_1_4[j] / 6;
      }
    }
  }
}

double FFParticle::GetEpsilon(const uint i, const uint j) const {
  uint idx = FlatIndex(i, j);
  return epsilon[idx];
}
double FFParticle::GetEpsilon_1_4(const uint i, const uint j) const {
  uint idx = FlatIndex(i, j);
  return epsilon_1_4[idx];
}
double FFParticle::GetSigma(const uint i, const uint j) const {
  uint idx = FlatIndex(i, j);
  return sqrt(sigmaSq[idx]);
}
double FFParticle::GetSigma_1_4(const uint i, const uint j) const {
  uint idx = FlatIndex(i, j);
  return sqrt(sigmaSq_1_4[idx]);
}
double FFParticle::GetN(const uint i, const uint j) const {
  uint idx = FlatIndex(i, j);
  return n[idx];
}
double FFParticle::GetN_1_4(const uint i, const uint j) const {
  uint idx = FlatIndex(i, j);
  return n_1_4[idx];
}
inline double FFParticle::GetRmin(const uint i, const uint j) const {
  return 0.0;
}
inline double FFParticle::GetRmax(const uint i, const uint j) const {
  return 0.0;
}
inline double FFParticle::GetRmin_1_4(const uint i, const uint j) const {
  return 0.0;
}
inline double FFParticle::GetRmax_1_4(const uint i, const uint j) const {
  return 0.0;
}

// Defining the functions

inline void FFParticle::CalcAdd_1_4(double &en, const double distSq,
                                    const uint kind1, const uint kind2) const {
  if (forcefield.rCutSq < distSq)
    return;

  uint index = FlatIndex(kind1, kind2);
  double rRat2 = sigmaSq_1_4[index] / distSq;
  double attract = rRat2 * rRat2 * rRat2;
  double repulse = pow(sqrt(rRat2), n_1_4[index]);

  en += epsilon_cn_1_4[index] * (repulse - attract);
}

inline void FFParticle::CalcCoulombAdd_1_4(double &en, const double distSq,
                                           const double qi_qj_Fact,
                                           const bool NB) const {
  if (forcefield.rCutSq < distSq)
    return;

  double dist = sqrt(distSq);
  if (NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * forcefield.scaling_14 / dist;
}

// mie potential
inline double FFParticle::CalcEn(const double distSq, const uint kind1,
                                 const uint kind2, const double lambda) const {
  if (forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  if (lambda >= 0.999999) {
    // save computation time
    return CalcEn(distSq, index);
  }
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef =
      forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);

  double en = lambda * CalcEn(softRsq, index);
  return en;
}

inline double FFParticle::CalcEn(const double distSq, const uint index) const {
  double rRat2 = sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double n_ij = n[index];
  double repulse = pow(rRat2, (n_ij * 0.5));

  return (epsilon_cn[index] * (repulse - attract));
}

inline double FFParticle::CalcVir(const double distSq, const uint kind1,
                                  const uint kind2, const double lambda) const {
  if (forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  if (lambda >= 0.999999) {
    // save computation time
    return CalcVir(distSq, index);
  }
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef =
      forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double correction = distSq / softRsq;
  // We need to fix the return value from calcVir
  double vir = lambda * correction * correction * CalcVir(softRsq, index);
  return vir;
}

inline double FFParticle::CalcVir(const double distSq, const uint index) const {
  double rNeg2 = 1.0 / distSq;
  double rRat2 = rNeg2 * sigmaSq[index];
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double n_ij = n[index];
  double repulse = pow(rRat2, (n_ij * 0.5));
  // Virial is F.r = -dE/dr * 1/r
  return epsilon_cn_6[index] * (nOver6[index] * repulse - attract) * rNeg2;
}

inline double FFParticle::CalcCoulomb(const double distSq, const uint kind1,
                                      const uint kind2, const double qi_qj_Fact,
                                      const double lambda, const uint b) const {
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if (lambda >= 0.999999) {
    // save computation time
    return CalcCoulomb(distSq, qi_qj_Fact, b);
  }
  double en = 0.0;
  if (forcefield.sc_coul) {
    uint index = FlatIndex(kind1, kind2);
    double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
    sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    en = lambda * CalcCoulomb(softRsq, qi_qj_Fact, b);
  } else {
    en = lambda * CalcCoulomb(distSq, qi_qj_Fact, b);
  }
  return en;
}

inline double FFParticle::CalcCoulomb(const double distSq,
                                      const double qi_qj_Fact,
                                      const uint b) const {
  double dist = sqrt(distSq);
  double val = qi_qj_Fact;
  if (forcefield.ewald) {
    val *= erfc(forcefield.alpha[b] * dist);
  }
  return val / dist;
}

inline double FFParticle::CalcCoulombVir(const double distSq, const uint kind1,
                                         const uint kind2, const double qi_qj,
                                         const double lambda,
                                         const uint b) const {
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if (lambda >= 0.999999) {
    // save computation time
    return CalcCoulombVir(distSq, qi_qj, b);
  }
  double vir = 0.0;
  if (forcefield.sc_coul) {
    uint index = FlatIndex(kind1, kind2);
    double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
    sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double correction = distSq / softRsq;
    // We need to fix the return value from calcVir
    vir = lambda * correction * correction * CalcCoulombVir(softRsq, qi_qj, b);
  } else {
    vir = lambda * CalcCoulombVir(distSq, qi_qj, b);
  }
  return vir;
}

inline double FFParticle::CalcCoulombVir(const double distSq,
                                         const double qi_qj,
                                         const uint b) const {
  const double dist = std::sqrt(distSq);
  if (forcefield.ewald) {
    // M_2_SQRTPI is 2/sqrt(PI)
    double alphaValue = forcefield.alpha[b] * M_2_SQRTPI;
    double expValue = std::exp(-forcefield.alphaSq[b] * distSq);
    double erfcValue = std::erfc(forcefield.alpha[b] * dist) / dist;
    return qi_qj * (alphaValue * expValue + erfcValue) / distSq;
  } else {
    return qi_qj / (distSq * dist);
  }
}

// Calculate the dE/dlambda for vdw energy
inline double FFParticle::CalcdEndL(const double distSq, const uint kind1,
                                    const uint kind2,
                                    const double lambda) const {
  if (forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef =
      forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = cbrt(softDist6);
  double fCoef = lambda * forcefield.sc_alpha * forcefield.sc_power / 6.0;
  fCoef *= pow(1.0 - lambda, forcefield.sc_power - 1.0) * sigma6 /
           (softRsq * softRsq);
  double dhdl = CalcEn(softRsq, index) + fCoef * CalcVir(softRsq, index);
  return dhdl;
}

// Calculate the dE/dlambda for Coulomb energy
inline double FFParticle::CalcCoulombdEndL(const double distSq,
                                           const uint kind1, const uint kind2,
                                           const double qi_qj_Fact,
                                           const double lambda, uint b) const {
  if (forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  double dhdl = 0.0;
  if (forcefield.sc_coul) {
    uint index = FlatIndex(kind1, kind2);
    double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
    sigma6 = std::max(sigma6, forcefield.sc_sigma_6);
    double dist6 = distSq * distSq * distSq;
    double lambdaCoef =
        forcefield.sc_alpha * pow((1.0 - lambda), forcefield.sc_power);
    double softDist6 = lambdaCoef * sigma6 + dist6;
    double softRsq = cbrt(softDist6);
    double fCoef = lambda * forcefield.sc_alpha * forcefield.sc_power / 6.0;
    fCoef *= pow(1.0 - lambda, forcefield.sc_power - 1.0) * sigma6 /
             (softRsq * softRsq);
    dhdl = CalcCoulomb(softRsq, qi_qj_Fact, b) +
           fCoef * CalcCoulombVir(softRsq, qi_qj_Fact, b);
  } else {
    dhdl = CalcCoulomb(distSq, qi_qj_Fact, b);
  }
  return dhdl;
}