#ifndef FF_SHIFT_H
#define FF_SHIFT_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "BasicTypes.h" //for uint
#include "NumLib.h" //For Cb, Sq
#include "FFParticle.h"

//////////////////////////////////////////////////////////////////////
////////////////////////// LJ Shift Style ////////////////////////////
//////////////////////////////////////////////////////////////////////
// Virial and LJ potential calculation:
// U(rij) = cn * eps_ij * ( (sig_ij/rij)^n - (sig_ij/rij)^6) + shiftConst
// shiftConst = cn * eps_ij * ( (sig_ij/rcut)^n - (sig_ij/rcut)^6)
// cn = n/(n-6) * ((n/6)^(6/(n-6)))
//
// Vir(r) = cn * eps_ij * 6 * ((n/6) * repulse - attract)/rij^2
// U_lrc = 0
// Vir_lrc = 0
//
// Eelect = qi * qj * (1/r - 1/rcut)
// Welect = qi * qj * 1/rij^3

struct FF_SHIFT : public FFParticle
{
public:

  virtual double CalcEn(const double distSq,
                        const uint kind1, const uint kind2) const;
  virtual double CalcVir(const double distSq,
                         const uint kind1, const uint kind2) const;
  virtual void CalcAdd_1_4(double& en, const double distSq,
                           const uint kind1, const uint kind2) const;

  // coulomb interaction functions
  virtual double CalcCoulombEn(const double distSq,
                               const double qi_qj_Fact) const;
  virtual double CalcCoulombVir(const double distSq,
                                const double qi_qj_Fact) const;
  virtual void CalcCoulombAdd_1_4(double& en, const double distSq,
                                  const double qi_qj_Fact) const;

  //!Returns Ezero, no energy correction
  virtual double EnergyLRC(const uint kind1, const uint kind2) const
  {
    return 0.0;
  }
  //!!Returns Ezero, no virial correction
  virtual double VirialLRC(const uint kind1, const uint kind2) const
  {
    return 0.0;
  }
};


inline void FF_SHIFT::CalcAdd_1_4(double& en, const double distSq,
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

  en += (epsilon_cn_1_4[index] * (repulse-attract) - shiftConst_1_4[index]);
}

inline void FF_SHIFT::CalcCoulombAdd_1_4(double& en, const double distSq,
    const double qi_qj_Fact) const
{
  double dist = sqrt(distSq);
  if(!ewald)
  {
     en += scaling_14 * qi_qj_Fact * (1.0/dist - 1.0/rCut);
  }
  else
  {
    double erfc = alpha * dist;
    en += scaling_14 * qi_qj_Fact * (1 - erf(erfc))/ dist;
  }
}


//mie potential
inline double FF_SHIFT::CalcEn(const double distSq,
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

  return (epsilon_cn[index] * (repulse-attract) - shiftConst[index]);
}

inline double FF_SHIFT::CalcCoulombEn(const double distSq,
                                      const double qi_qj_Fact) const
{
  double dist = sqrt(distSq);
  if(!ewald)
  {
     return  qi_qj_Fact * (1.0/dist - 1.0/rCut);
  }
  else
  {
     double erfc = alpha * dist;
     return  qi_qj_Fact * (1 - erf(erfc))/ dist;
  }
}

//mie potential
inline double FF_SHIFT::CalcVir(const double distSq,
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

inline double FF_SHIFT::CalcCoulombVir(const double distSq,
                                       const double qi_qj) const
{
  double dist = sqrt(distSq);
  if(!ewald)
  { 
     return qi_qj/(distSq * dist);
  }
  else
  {
     double constValue = 2.0 * alpha / sqrt(M_PI);
     double expConstValue = exp(-1.0 * alpha * alpha * distSq);
     double temp = 1.0 - erf(alpha * dist);
     return  qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  }
}


#endif /*FF_SHIFT_H*/
