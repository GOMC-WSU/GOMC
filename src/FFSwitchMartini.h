/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef FF_SWITCH_MARTINI_H
#define FF_SWITCH_MARTINI_H

#include "BasicTypes.h" //for uint
#include "FFConst.h"    //constants related to particles.
#include "FFParticle.h"
#include "NumLib.h" //For Cb, Sq

///////////////////////////////////////////////////////////////////////
////////////////////////// LJ Switch Martini Style ////////////////////////////
///////////////////////////////////////////////////////////////////////
// LJ potential calculation:
// Eij = cn * eps_ij * ( sig_ij^n * (1/rij^n + phi(n)) -
//       sig_ij^6 * (1/rij^6 + phi(6)))
// cn = n/(n-6) * ((n/6)^(6/(n-6)))
//
// Eelec = qi*qj*(1/rij + phi(1))
// Welec = qi*qj*(1/rij^3 + phiW(1)/r)
//
// phi(x) = -Cx , if r < rswitch
// phi(x) = -Ax *(r - rswitch)^3/3 - Bx * (r - rswitch)^4 / 4 - Cx ,if r>rswitch
//
// Ax = x * ((x + 1) * rswitch - (x + 4) * rcut) /
//      (rcut^(x + 2)) * (rcut - rswitch)^2)
// Bx = x * ((x + 1) * rswitch - (x + 3) * rcut) /
//      (rcut^(x + 2)) * (rcut - rswitch)^3)
// Cx = 1/rcut^x - Ax * (rcut - rswitch)^3 / 3 - Bx * (rcut - rswitch)^4 / 4
//
// Virial Calculation
//
// Wij = cn * eps_ij * ( sig_ij^n * (n/rij^(n+2) + phiW(n)/rij) -
//       sig_ij^6 * (6/rij^(6+2) + phiW(6)/r))
//
// phiW(x) = 0 , if r < rswitch
// phiW(x) = Ax *(r - rswitch)^2 + Bx * (r - rswitch)^3 ,if r > rswitch
//
//

struct FF_SWITCH_MARTINI : public FFParticle {
public:
  FF_SWITCH_MARTINI(Forcefield &ff)
      : FFParticle(ff), An(NULL), Bn(NULL), Cn(NULL), An_1_4(NULL),
        Bn_1_4(NULL), Cn_1_4(NULL), sig6(NULL), sign(NULL), sig6_1_4(NULL),
        sign_1_4(NULL) {
    A1 = B1 = C1 = A6 = B6 = C6 = 0.0;
  }
  virtual ~FF_SWITCH_MARTINI() {
    delete[] An;
    delete[] Bn;
    delete[] Cn;
    delete[] An_1_4;
    delete[] Bn_1_4;
    delete[] Cn_1_4;
    delete[] sig6;
    delete[] sign;
    delete[] sig6_1_4;
    delete[] sign_1_4;
  }

  virtual void Init(ff_setup::Particle const &mie,
                    ff_setup::NBfix const &nbfix);

  virtual double CalcEn(const double distSq, const uint kind1, const uint kind2,
                        const double lambda) const;
  virtual double CalcVir(const double distSq, const uint kind1,
                         const uint kind2, const double lambda) const;
  virtual void CalcAdd_1_4(double &en, const double distSq, const uint kind1,
                           const uint kind2) const;

  // coulomb interaction functions
  virtual double CalcCoulomb(const double distSq, const uint kind1,
                             const uint kind2, const double qi_qj_Fact,
                             const double lambda, const uint b) const;
  virtual double CalcCoulombVir(const double distSq, const uint kind1,
                                const uint kind2, const double qi_qj,
                                const double lambda, const uint b) const;
  virtual void CalcCoulombAdd_1_4(double &en, const double distSq,
                                  const double qi_qj_Fact, const bool NB) const;

  //! Returns Ezero, no energy correction
  virtual double EnergyLRC(const uint kind1, const uint kind2) const {
    return 0.0;
  }
  //!!Returns Ezero, no virial correction
  virtual double VirialLRC(const uint kind1, const uint kind2) const {
    return 0.0;
  }
  //! Returns zero for impulse pressure correction term for a kind pair
  virtual double ImpulsePressureCorrection(const uint kind1,
                                           const uint kind2) const {
    return 0.0;
  }

  // Calculate the dE/dlambda for vdw energy
  virtual double CalcdEndL(const double distSq, const uint kind1,
                           const uint kind2, const double lambda) const;
  // Calculate the dE/dlambda for Coulomb energy
  virtual double CalcCoulombdEndL(const double distSq, const uint kind1,
                                  const uint kind2, const double qi_qj_Fact,
                                  const double lambda, uint b) const;

protected:
  virtual double CalcEn(const double distSq, const uint index) const;
  virtual double CalcVir(const double distSq, const uint index) const;
  virtual double CalcCoulomb(const double distSq, const double qi_qj_Fact,
                             const uint b) const;
  virtual double CalcCoulombVir(const double distSq, const double qi_qj,
                                uint b) const;

  double *An, *Bn, *Cn, *An_1_4, *Bn_1_4, *Cn_1_4;
  double *sig6, *sign, *sig6_1_4, *sign_1_4;

  double diElectric_1, rOn, rOnSq, rOnCoul, A1, B1, C1, A6, B6, C6;
};

inline void FF_SWITCH_MARTINI::Init(ff_setup::Particle const &mie,
                                    ff_setup::NBfix const &nbfix) {
  // Initialize sigma and epsilon
  FFParticle::Init(mie, nbfix);
  uint size = num::Sq(count);
  // allocate memory
  An = new double[size];
  Bn = new double[size];
  Cn = new double[size];
  sign = new double[size];
  sig6 = new double[size];
  An_1_4 = new double[size];
  Bn_1_4 = new double[size];
  Cn_1_4 = new double[size];
  sign_1_4 = new double[size];
  sig6_1_4 = new double[size];
  // Set martini constant
  diElectric_1 = 1.0 / forcefield.dielectric;
  rOn = forcefield.rswitch;
  rOnSq = rOn * rOn;
  // in Martini, Coulomb switching distance is zero
  rOnCoul = 0.0;
  double rCut = forcefield.rCut;
  // Original Unoptimized computation
  // Set LJ constants
  // A6 = 6.0 * ((6.0 + 1) * rOn - (6.0 + 4) * rCut) / (pow(rCut, 6.0 + 2) *
  // pow(rCut - rOn, 2));
  // B6 = -6.0 * ((6.0 + 1) * rOn - (6.0 + 3) * rCut) / (pow(rCut, 6.0 + 2) *
  // pow(rCut - rOn, 3));
  // C6 = 1.0 / pow(rCut, 6.0) - A6 / 3.0 * pow(rCut - rOn, 3) - B6 / 4.0 *
  // pow(rCut - rOn, 4);
  // // Set Coulomb constants
  // A1 = 1.0 * ((1.0 + 1) * rOnCoul - (1.0 + 4) * rCut) / (pow(rCut, 1.0 + 2) *
  // pow(rCut - rOnCoul, 2));
  // B1 = -1.0 * ((1.0 + 1) * rOnCoul - (1.0 + 3) * rCut) / (pow(rCut, 1.0 + 2)
  // * pow(rCut - rOnCoul, 3));
  // C1 = 1.0 / pow(rCut, 1.0) - A1 / 3.0 * pow(rCut - rOnCoul, 3) - B1 / 4.0 *
  // pow(rCut - rOnCoul, 4);

  // Optimized computation
  // Set LJ constants
  A6 = 6.0 * (7.0 * rOn - 10.0 * rCut) /
       (pow(rCut, 8.0) * (rCut - rOn) * (rCut - rOn));
  B6 = -6.0 * (7.0 * rOn - 9.0 * rCut) /
       (pow(rCut, 8.0) * (rCut - rOn) * (rCut - rOn) * (rCut - rOn));
  C6 = pow(rCut, -6.0) - A6 / 3.0 * (rCut - rOn) * (rCut - rOn) * (rCut - rOn) -
       B6 / 4.0 * (rCut - rOn) * (rCut - rOn) * (rCut - rOn) * (rCut - rOn);
  // Set Coulomb constants
  A1 = (2.0 * rOnCoul - 5.0 * rCut) /
       (rCut * rCut * rCut * (rCut - rOnCoul) * (rCut - rOnCoul));
  B1 = -1.0 * (2.0 * rOnCoul - 4.0 * rCut) /
       (rCut * rCut * rCut * (rCut - rOnCoul) * (rCut - rOnCoul) *
        (rCut - rOnCoul));
  C1 = 1.0 / rCut -
       A1 / 3.0 * (rCut - rOnCoul) * (rCut - rOnCoul) * (rCut - rOnCoul) -
       B1 / 4.0 * (rCut - rOnCoul) * (rCut - rOnCoul) * (rCut - rOnCoul) *
           (rCut - rOnCoul);

  for (uint i = 0; i < count; ++i) {
    for (uint j = 0; j < count; ++j) {
      uint idx = FlatIndex(i, j);
      double pn = n[idx];
      An[idx] = pn * ((pn + 1.0) * rOn - (pn + 4.0) * rCut) /
                (pow(rCut, pn + 2.0) * (rCut - rOn) * (rCut - rOn));
      Bn[idx] =
          -pn * ((pn + 1.0) * rOn - (pn + 3.0) * rCut) /
          (pow(rCut, pn + 2.0) * (rCut - rOn) * (rCut - rOn) * (rCut - rOn));
      Cn[idx] = 1.0 / pow(rCut, pn) -
                An[idx] / 3.0 * (rCut - rOn) * (rCut - rOn) * (rCut - rOn) -
                Bn[idx] / 4.0 * (rCut - rOn) * (rCut - rOn) * (rCut - rOn) *
                    (rCut - rOn);
      double sigma = sqrt(sigmaSq[idx]);
      sig6[idx] = pow(sigma, 6.0);
      sign[idx] = pow(sigma, pn);

      // for 1-4 interaction
      double pn_1_4 = n_1_4[idx];
      An_1_4[idx] = pn_1_4 * ((pn_1_4 + 1.0) * rOn - (pn_1_4 + 4.0) * rCut) /
                    (pow(rCut, pn_1_4 + 2.0) * (rCut - rOn) * (rCut - rOn));
      Bn_1_4[idx] = -pn_1_4 * ((pn_1_4 + 1.0) * rOn - (pn_1_4 + 3.0) * rCut) /
                    (pow(rCut, pn_1_4 + 2.0) * (rCut - rOn) * (rCut - rOn) *
                     (rCut - rOn));
      Cn_1_4[idx] =
          1.0 / pow(rCut, pn_1_4) -
          An_1_4[idx] / 3.0 * (rCut - rOn) * (rCut - rOn) * (rCut - rOn) -
          Bn_1_4[idx] / 4.0 * (rCut - rOn) * (rCut - rOn) * (rCut - rOn) *
              (rCut - rOn);
      double sigma_1_4 = sqrt(sigmaSq_1_4[idx]);
      sig6_1_4[idx] = pow(sigma_1_4, 6.0);
      sign_1_4[idx] = pow(sigma_1_4, pn_1_4);
    }
  }
}

inline void FF_SWITCH_MARTINI::CalcAdd_1_4(double &en, const double distSq,
                                           const uint kind1,
                                           const uint kind2) const {
  if (forcefield.rCutSq < distSq)
    return;

  uint index = FlatIndex(kind1, kind2);
  double r_2 = 1.0 / distSq;
  double r_6 = r_2 * r_2 * r_2;
  double r_n = pow(sqrt(r_2), n_1_4[index]);

  double rij_ron = sqrt(distSq) - rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double shifttempRep = -(An_1_4[index] / 3.0) * rij_ron_3 -
                        (Bn_1_4[index] / 4.0) * rij_ron_4 - Cn_1_4[index];
  double shifttempAtt = -(A6 / 3.0) * rij_ron_3 - (B6 / 4.0) * rij_ron_4 - C6;

  const double shiftRep = (distSq > rOnSq ? shifttempRep : -Cn_1_4[index]);
  const double shiftAtt = (distSq > rOnSq ? shifttempAtt : -C6);

  en += epsilon_cn_1_4[index] * (sign_1_4[index] * (r_n + shiftRep) -
                                 sig6_1_4[index] * (r_6 + shiftAtt));
}

inline void FF_SWITCH_MARTINI::CalcCoulombAdd_1_4(double &en,
                                                  const double distSq,
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

inline double FF_SWITCH_MARTINI::CalcEn(const double distSq, const uint kind1,
                                        const uint kind2,
                                        const double lambda) const {
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

inline double FF_SWITCH_MARTINI::CalcEn(const double distSq,
                                        const uint index) const {
  double r_2 = 1.0 / distSq;
  double r_4 = r_2 * r_2;
  double r_6 = r_4 * r_2;
  double n_ij = n[index];
  double r_n = pow(r_2, (n_ij * 0.5));

  double rij_ron = sqrt(distSq) - rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double shifttempRep = -(An[index] / 3.0) * rij_ron_3 -
                        (Bn[index] / 4.0) * rij_ron_4 - Cn[index];
  double shifttempAtt = -(A6 / 3.0) * rij_ron_3 - (B6 / 4.0) * rij_ron_4 - C6;

  const double shiftRep = (distSq > rOnSq ? shifttempRep : -Cn[index]);
  const double shiftAtt = (distSq > rOnSq ? shifttempAtt : -C6);

  double Eij = epsilon_cn[index] * (sign[index] * (r_n + shiftRep) -
                                    sig6[index] * (r_6 + shiftAtt));
  return Eij;
}

inline double FF_SWITCH_MARTINI::CalcVir(const double distSq, const uint kind1,
                                         const uint kind2,
                                         const double lambda) const {
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

inline double FF_SWITCH_MARTINI::CalcVir(const double distSq,
                                         const uint index) const {
  double n_ij = n[index];
  double r_1 = 1.0 / sqrt(distSq);
  double r_8 = distSq * distSq * distSq * distSq;
  double r_n2 = pow(r_1, n_ij + 2.0);

  double rij_ron = sqrt(distSq) - rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;

  double dshifttempRep = An[index] * rij_ron_2 + Bn[index] * rij_ron_3;
  double dshifttempAtt = A6 * rij_ron_2 + B6 * rij_ron_3;

  const double dshiftRep = (distSq > rOnSq ? dshifttempRep * r_1 : 0);
  const double dshiftAtt = (distSq > rOnSq ? dshifttempAtt * r_1 : 0);

  double Wij = epsilon_cn[index] * (sign[index] * (n_ij * r_n2 + dshiftRep) -
                                    sig6[index] * (6.0 * r_8 + dshiftAtt));
  return Wij;
}

inline double FF_SWITCH_MARTINI::CalcCoulomb(const double distSq,
                                             const uint kind1, const uint kind2,
                                             const double qi_qj_Fact,
                                             const double lambda,
                                             const uint b) const {
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

inline double FF_SWITCH_MARTINI::CalcCoulomb(const double distSq,
                                             const double qi_qj_Fact,
                                             const uint b) const {
  if (forcefield.ewald) {
    double dist = sqrt(distSq);
    double val = forcefield.alpha[b] * dist;
    return qi_qj_Fact * erfc(val) / dist;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    double dist = sqrt(distSq);
    double rij_ronCoul_3 = dist * distSq;
    double rij_ronCoul_4 = distSq * distSq;

    double coul = -(A1 / 3.0) * rij_ronCoul_3 - (B1 / 4.0) * rij_ronCoul_4 - C1;
    return qi_qj_Fact * diElectric_1 * (1.0 / dist + coul);
  }
}

inline double
FF_SWITCH_MARTINI::CalcCoulombVir(const double distSq, const uint kind1,
                                  const uint kind2, const double qi_qj,
                                  const double lambda, const uint b) const {
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

inline double FF_SWITCH_MARTINI::CalcCoulombVir(const double distSq,
                                                const double qi_qj,
                                                const uint b) const {
  if (forcefield.ewald) {
    double dist = sqrt(distSq);
    // M_2_SQRTPI is 2/sqrt(PI)
    double constValue = forcefield.alpha[b] * M_2_SQRTPI;
    double expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    double temp = erfc(forcefield.alpha[b] * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    // in Martini, the Coulomb switching distance is zero, so we will have
    // sqrt(distSq) - rOnCoul =  sqrt(distSq)
    double dist = sqrt(distSq);
    double rij_ronCoul_2 = distSq;
    double rij_ronCoul_3 = dist * distSq;

    double virCoul = A1 / rij_ronCoul_2 + B1 / rij_ronCoul_3;
    return qi_qj * diElectric_1 * (1.0 / (dist * distSq) + virCoul / dist);
  }
}

inline double FF_SWITCH_MARTINI::CalcdEndL(const double distSq,
                                           const uint kind1, const uint kind2,
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
inline double
FF_SWITCH_MARTINI::CalcCoulombdEndL(const double distSq, const uint kind1,
                                    const uint kind2, const double qi_qj_Fact,
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

#endif /*FF_SWITCH_MARTINI_H*/
