/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef FF_PARTICLE_H
#define FF_PARTICLE_H

#include "BasicTypes.h" //for uint
#include "FFConst.h"    //constants related to particles.
#include "Forcefield.h"
#include "NumLib.h" //For Cb, Sq
#include "Setup.h"
#ifdef GOMC_CUDA
#include "VariablesCUDA.cuh"
#endif

// Virial and LJ potential calculation:
// U(rij) = cn * eps_ij * ( (sig_ij/rij)^n - (sig_ij/rij)^6)
//
// cn = n/(n-6) * ((n/6)^(6/(n-6)))
//
// eps_E_cn = cn * eps_ij
//                  __________const__________________
// U_lrc = density * 0.5 * 4.0 / (n-3) * cn * pi * eps_ij * sig_ij^3 *
//          ( (sig_ij/rij)^(n-3) - (n-3)/3*(sig_ij/rij)^3)
//
// Vir(r) = cn * eps_ij * n * (sig_ij/rij)^n - cn * eps_ij * 6 * (sig_ij/rij)^6
// Vir(r) = cn * eps_ij * n * repulse - cn * eps_ij * 6 * attract
// Vir(r) = cn * eps_ij * (n * repulse - 6 * attract)
// Vir(r) = cn * eps_ij * 6 * ((n/6) * repulse - attract)
//
// Vir_lrc = density * 0.5 * 4.0 * 2/3 * cn * pi * eps_ij * sig_ij^3 *
//          ( n/(n-3) * 3/2 * (sig_ij/rij)^(n-3) - 3*(sig_ij/rij)^3)

namespace ff_setup {
class Particle;
class NBfix;
} // namespace ff_setup

class Forcefield;

struct FFParticle {
public:
  FFParticle(Forcefield &ff);
  virtual ~FFParticle(void);

  virtual void Init(ff_setup::Particle const &mie,
                    ff_setup::NBfix const &nbfix);

  double GetEpsilon(const uint i, const uint j) const;
  double GetEpsilon_1_4(const uint i, const uint j) const;
  double GetSigma(const uint i, const uint j) const;
  double GetSigma_1_4(const uint i, const uint j) const;
  double GetN(const uint i, const uint j) const;
  double GetN_1_4(const uint i, const uint j) const;
  virtual double GetRmin(const uint i, const uint j) const;
  virtual double GetRmax(const uint i, const uint j) const;
  virtual double GetRmin_1_4(const uint i, const uint j) const;
  virtual double GetRmax_1_4(const uint i, const uint j) const;

  // LJ interaction functions
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
                                const double lambda, uint b) const;
  virtual void CalcCoulombAdd_1_4(double &en, const double distSq,
                                  const double qi_qj_Fact, const bool NB) const;

  //! Returns Energy long-range correction term for a kind pair
  virtual double EnergyLRC(const uint kind1, const uint kind2) const;
  //! Returns Energy long-range correction term for a kind pair
  virtual double VirialLRC(const uint kind1, const uint kind2) const;
  //! Returns impulse pressure correction term for a kind pair
  virtual double ImpulsePressureCorrection(const uint kind1,
                                           const uint kind2) const;

  // Calculate the dE/dlambda for vdw energy
  virtual double CalcdEndL(const double distSq, const uint kind1,
                           const uint kind2, const double lambda) const;
  // Calculate the dE/dlambda for Coulomb energy
  virtual double CalcCoulombdEndL(const double distSq, const uint kind1,
                                  const uint kind2, const double qi_qj_Fact,
                                  const double lambda, uint b) const;

  uint NumKinds() const { return count; }
  double GetMass(const uint kind) const { return mass[kind]; }

#ifdef GOMC_CUDA
  VariablesCUDA *getCUDAVars() { return varCUDA; }
#endif

protected:
  virtual double CalcEn(const double distSq, const uint index) const;
  virtual double CalcVir(const double distSq, const uint index) const;
  virtual double CalcCoulomb(const double distSq, const double qi_qj_Fact,
                             const uint b) const;
  virtual double CalcCoulombVir(const double distSq, const double qi_qj,
                                uint b) const;
  // Find the index of the pair kind
  uint FlatIndex(const uint i, const uint j) const { return i + j * count; }
  // Combining sigma, epsilon, and n value for different kind
  void Blend(ff_setup::Particle const &mie);
  // Use NBFIX to adjust sigma, epsilon, and n value for different kind
  void AdjNBfix(ff_setup::NBfix const &nbfix);
  // To access rcut and other forcefield data
  const Forcefield &forcefield;

  double *mass;
  std::string *nameFirst;
  std::string *nameSec;

  // vars for LJ-LJ pairs
  double *n, *n_1_4;
  // For LJ eps_cn(en) --> 4eps, eps_cn_6 --> 24eps, eps_cn_n --> 48eps
  double *sigmaSq, *sigmaSq_1_4, *epsilon, *epsilon_1_4, *epsilon_cn,
      *epsilon_cn_1_4, *epsilon_cn_6, *epsilon_cn_6_1_4, *nOver6, *nOver6_1_4;
#ifdef GOMC_CUDA
  VariablesCUDA *varCUDA;
#endif
  uint count;
  bool exp6;
};

#endif /*FF_PARTICLE_H*/
