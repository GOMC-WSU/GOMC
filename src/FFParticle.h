/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.11
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_PARTICLE_H
#define FF_PARTICLE_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "BasicTypes.h" //for uint
#include "NumLib.h" //For Cb, Sq
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

namespace ff_setup
{
class Particle;
class NBfix;
}
namespace config_setup
{
  struct SystemVals;
  struct FFValues;
  struct FFKind;
}

struct FFParticle
{
public:

  FFParticle();
  ~FFParticle(void);

  double GetMass(const uint kind) const
  {
    return mass[kind];
  }
  // LJ interaction functions
  virtual double CalcEn(const double distSq,
                        const uint kind1, const uint kind2) const;
  virtual double CalcVir(const double distSq,
                         const uint kind1, const uint kind2) const;
  virtual void CalcAdd_1_4(double& en, const double distSq,
                           const uint kind1, const uint kind2) const;

  // coulomb interaction functions
  virtual double CalcCoulomb(const double distSq,
			     const double qi_qj_Fact) const;
  virtual double CalcCoulombEn(const double distSq,
                               const double qi_qj_Fact) const;
  virtual double CalcCoulombVir(const double distSq,
                                const double qi_qj) const;
  virtual void CalcCoulombAdd_1_4(double& en, const double distSq,
                                  const double qi_qj_Fact,
				  const bool NB) const;

  void Init(ff_setup::Particle const& mie,
            ff_setup::NBfix const& nbfix,
            config_setup::SystemVals const& sys,
            config_setup::FFKind const& ffKind);
  //!Returns Energy long-range correction term for a kind pair
  virtual double EnergyLRC(const uint kind1, const uint kind2) const;
  //!Returns Energy long-range correction term for a kind pair
  virtual double VirialLRC(const uint kind1, const uint kind2) const;

  uint NumKinds() const
  {
    return count;
  }

#ifdef GOMC_CUDA
  VariablesCUDA *getCUDAVars() {
    return varCUDA;
  }
#endif

protected:

  uint FlatIndex(const uint i, const uint j) const
  {
    return i + j * count;
  }
  void Blend(ff_setup::Particle const& mie, const double rCut);
  void AdjNBfix(ff_setup::Particle const& mie, ff_setup::NBfix const& nbfix,
                const double rCut);

  //vars for lj particles.
  double* mass;
  std::string *nameFirst;
  std::string *nameSec;

  //vars for LJ-LJ pairs
#ifdef MIE_INT_ONLY
  uint* n, *n_1_4;
#else
  double *n, *n_1_4;
#endif
  //For LJ eps_cn(en) --> 4eps, eps_cn_6 --> 24eps, eps_cn_n --> 48eps
  double * sigmaSq, * epsilon_cn, * epsilon_cn_6, * nOver6,
         * sigmaSq_1_4, * epsilon_cn_1_4, * epsilon_cn_6_1_4, * nOver6_1_4,
         * enCorrection, * virCorrection, *shiftConst, *An, *Bn, *Cn, *sig6, *sign,
         *shiftConst_1_4, *An_1_4, *Bn_1_4, *Cn_1_4, *sig6_1_4, *sign_1_4;

  double rCut, rCutSq, rOn, rOnSq, rOnCoul, A1, B1, C1, A6, B6, C6,
         factor1, factor2, scaling_14, alpha, diElectric_1;
  double rCutLow, rCutLowSq;

  uint count, vdwKind;
  bool isMartini, ewald;
#ifdef GOMC_CUDA
  VariablesCUDA *varCUDA;
#endif
};



inline void FFParticle::CalcAdd_1_4(double& en, const double distSq,
                                    const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  double rRat2 = sigmaSq_1_4[index]/distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n_1_4[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n_1_4[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  en += epsilon_cn_1_4[index] * (repulse-attract);
}

inline void FFParticle::CalcCoulombAdd_1_4(double& en, const double distSq,
					   const double qi_qj_Fact,
					   const bool NB) const
{
  double dist = sqrt(distSq);
  if(NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * scaling_14 / dist;
}



//mie potential
inline double FFParticle::CalcEn(const double distSq,
                                 const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  double rRat2 = sigmaSq[index]/distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  return epsilon_cn[index] * (repulse-attract);
}

inline double FFParticle::CalcCoulomb(const double distSq,
				      const double qi_qj_Fact) const
{
  if(ewald)
  {
     double dist = sqrt(distSq);
     double val = alpha * dist;
     return  qi_qj_Fact * erfc(val)/ dist;
  }
  else
  {
     double dist = sqrt(distSq);
     return  qi_qj_Fact / dist;
  }
}

//will be used in energy calculation after each move
inline double FFParticle::CalcCoulombEn(const double distSq,
                                        const double qi_qj_Fact) const
{
  if(distSq <= rCutLowSq)
    return num::BIGNUM;

  if(ewald)
  {
     double dist = sqrt(distSq);
     double val = alpha * dist;
     return  qi_qj_Fact * erfc(val)/ dist;
  }
  else
  {
     double dist = sqrt(distSq);
     return  qi_qj_Fact / dist;
  }
}


inline double FFParticle::CalcVir(const double distSq,
                                  const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  double rNeg2 = 1.0/distSq;
  double rRat2 = rNeg2 * sigmaSq[index];
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  //Virial is the derivative of the pressure... mu
  return epsilon_cn_6[index] * (nOver6[index]*repulse-attract)*rNeg2;
}

inline double FFParticle::CalcCoulombVir(const double distSq,
					 const double qi_qj) const
{
  if(ewald)
  {
    double dist = sqrt(distSq);
    double constValue = 2.0 * alpha / sqrt(M_PI);
    double expConstValue = exp(-1.0 * alpha * alpha * distSq);
    double temp = 1.0 - erf(alpha * dist);
    return  qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  }
  else
  {
     double dist = sqrt(distSq);
     return qi_qj/(distSq * dist);
  }
}

#endif /*FF_PARTICLE_H*/
