/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_SWITCH_H
#define FF_SWITCH_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "BasicTypes.h" //for uint
#include "NumLib.h" //For Cb, Sq
#include "FFParticle.h"


///////////////////////////////////////////////////////////////////////
////////////////////////// LJ Switch Style ////////////////////////////
///////////////////////////////////////////////////////////////////////
// LJ potential calculation:
// Eij = (cn * eps_ij * ( (sig_ij/rij)^n - (sig_ij/rij)^6)) * FE
// cn = n/(n-6) * ((n/6)^(6/(n-6)))
// FE = 1 , if rij < rswitch
// FE = (rcut^2 - rij^2)^2 *(factor1 + 2 * rij^2)* factor2 , if rcut>rij>rswitch
// factor1 = (rcut^2 - 3 * rswitch^2)
// factor2 = (rcut^2 - rswitch^2)^-3
//
// Virial calculation
// Vir(r) = FE * cn * eps_ij * 6 * ((n/6) * repulse - attract)/rij^2 -
//          cn * eps_ij * (repulse - attract) * FW
//
// FW = 0 , if rij < rswitch
// FW = 12 * (rcut^2 - rij^2)(rswitch^2-rij^2) * factor2 , if rcut>rij >rswitch
//
// Eelec = qi * qj * (rij^2/rcut^2 - 1.0)^2 / rij
// Welect = -1 * qi * qj * (dSwitchVal/rij^2 - (rij^2/rcut^2 - 1.0)^2/(rij^3))
// dSwitchVa = 2.0 * (rij^2/rcut^2 - 1.0) * 2.0 * rij/rcut^2

struct FF_SWITCH : public FFParticle {
public:

  FF_SWITCH(Forcefield &ff) : FFParticle(ff) 
  {
    rOnSq = rOn = factor1 = factor2 = 0.0;
  }

  virtual void Init(ff_setup::Particle const& mie,
                    ff_setup::NBfix const& nbfix);

  virtual double CalcEn(const double distSq,
                        const uint kind1, const uint kind2,
                        const double lambda) const;
  virtual double CalcVir(const double distSq, const uint kind1, 
                         const uint kind2, const double lambda) const;
  virtual void CalcAdd_1_4(double& en, const double distSq,
                           const uint kind1, const uint kind2) const;

  // coulomb interaction functions
  virtual double CalcCoulomb(const double distSq,
                             const double qi_qj_Fact, 
                             const double lambda,
                             const uint b) const;
  virtual double CalcCoulombVir(const double distSq, const double qi_qj, 
                                const double lambda, const uint b) const;
  virtual void CalcCoulombAdd_1_4(double& en, const double distSq,
                                  const double qi_qj_Fact,
                                  const bool NB) const;

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

  //Calculate the dE/dlambda for vdw energy
  virtual double CalcdEndL(const double distSq, const uint kind1,
                           const uint kind2, const 
                           double lambda) const;

  //Calculate the dE/dlambda for Coulomb energy
  virtual double CalcCoulombdEndL(const double distSq, const double qi_qj_Fact,
                                  const double lambda, uint b) const;

  protected:
  virtual double CalcEn(const double distSq, const uint index) const;
  virtual double CalcVir(const double distSq, const uint index) const;
  virtual double CalcCoulomb(const double distSq, const double qi_qj_Fact, 
                             const uint b) const;
  virtual double CalcCoulombVir(const double distSq, const double qi_qj,
                                uint b) const;

  double rOn, rOnSq, factor1, factor2;

};

inline void FF_SWITCH::Init(ff_setup::Particle const& mie,
                            ff_setup::NBfix const& nbfix)
{
  //Initializ sigma and epsilon
  FFParticle::Init(mie, nbfix);
  rOn = forcefield.rswitch;
  rOnSq = rOn * rOn;
  //calculate switch constant
  factor1 = forcefield.rCutSq - 3 * rOnSq;
  factor2 = pow((forcefield.rCutSq - rOnSq), -3);
}

inline void FF_SWITCH::CalcAdd_1_4(double& en, const double distSq,
                                   const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  double rCutSq_rijSq = forcefield.rCutSq - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  double rRat2 = sigmaSq_1_4[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n_1_4[index];
  double repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  double n_ij = n_1_4[index];
  double repulse = pow(sqrt(rRat2), n_ij);
#endif

  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

  const double factE = ( distSq > rOnSq ? fE : 1.0);

  en += (epsilon_cn_1_4[index] * (repulse - attract)) * factE;
}

inline void FF_SWITCH::CalcCoulombAdd_1_4(double& en, const double distSq,
                                          const double qi_qj_Fact,
                                          const bool NB) const
{
  double dist = sqrt(distSq);
  if(NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * forcefield.scaling_14 / dist;
}

inline double FF_SWITCH::CalcEn(const double distSq, const uint kind1,
                                const uint kind2,
                                const double lambda) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = 0.5 * (1.0 - lambda) * (1.0 - lambda);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, 1.0/3.0);

  double en = lambda * CalcEn(softRsq, index);
  return en;
}

inline double FF_SWITCH::CalcEn(const double distSq, const uint index) const
{
  double rCutSq_rijSq = forcefield.rCutSq - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;
  double rRat2 = sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double n_ij = n[index];
  double repulse = pow(rRat2, (n_ij * 0.5));

  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);
  const double factE = ( distSq > rOnSq ? fE : 1.0);

  return (epsilon_cn[index] * (repulse - attract)) * factE;
}

inline double FF_SWITCH::CalcVir(const double distSq, const uint kind1, 
                                 const uint kind2, const double lambda) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;
    
  uint index = FlatIndex(kind1, kind2);
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = 0.5 * (1.0 - lambda) * (1.0 - lambda);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, 1.0/3.0);
  double correction = distSq / softRsq;
  //We need to fix the return value from calcVir
  double vir = lambda * correction * correction * CalcVir(softRsq, index);
  return vir;
}

inline double FF_SWITCH::CalcVir(const double distSq, const uint index) const
{
  double rCutSq_rijSq = forcefield.rCutSq - distSq;
  double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  double rNeg2 = 1.0 / distSq;
  double rRat2 = rNeg2 * sigmaSq[index];
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double n_ij = n[index];
  double repulse = pow(rRat2, (n_ij * 0.5));

  double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);
  double fW = 12.0 * factor2 * rCutSq_rijSq * (rOnSq - distSq);

  const double factE = ( distSq > rOnSq ? fE : 1.0);
  const double factW = ( distSq > rOnSq ? fW : 0.0);
  double Wij = epsilon_cn_6[index] * (nOver6[index] * repulse - attract) * rNeg2;
  double Eij = epsilon_cn[index] * (repulse - attract);

  return (Wij * factE - Eij * factW);
}

inline double FF_SWITCH::CalcCoulomb(const double distSq,
                                    const double qi_qj_Fact, 
                                    const double lambda,
                                    const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  //Use amber scheme for soft core potential in Coulomb interaction
  // 12 is default value for Beta according to amber tutorial (Eq. 21.7)
  // http://ambermd.org/doc12/Amber18.pdf

  double lambdaCoef = 12.00 * (1.0 - lambda);
  double softRsq = lambdaCoef + distSq;

  double en = lambda * CalcCoulomb(softRsq, qi_qj_Fact, b);
  return en;
}

inline double FF_SWITCH::CalcCoulomb(const double distSq,
                                     const double qi_qj_Fact, 
                                     const uint b) const
{
  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double val = forcefield.alpha[b] * dist;
    return qi_qj_Fact * erfc(val) / dist;
  } else {
    double dist = sqrt(distSq);
    double switchVal = distSq / forcefield.rCutSq - 1.0;
    switchVal *= switchVal;
    return qi_qj_Fact * switchVal / dist;
  }
}

inline double FF_SWITCH::CalcCoulombVir(const double distSq,
                                        const double qi_qj, 
                                        const double lambda,
                                        const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  //Use amber scheme for soft core potential in Coulomb interaction
  // 12 is default value for Beta according to amber tutorial (Eq. 21.7)
  // http://ambermd.org/doc12/Amber18.pdf

  double lambdaCoef = 12.00 * (1.0 - lambda);
  double softRsq = lambdaCoef + distSq;
  //The only correction is to multiply by lambda
  double vir = lambda * CalcCoulombVir(softRsq, qi_qj, b);
  return vir;
}

inline double FF_SWITCH::CalcCoulombVir(const double distSq, const double qi_qj,
                                        const uint b) const
{
  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double constValue = 2.0 * forcefield.alpha[b] / sqrt(M_PI);
    double expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    double temp = erfc(forcefield.alpha[b] * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    double dist = sqrt(distSq);
    double switchVal = distSq / forcefield.rCutSq - 1.0;
    switchVal *= switchVal;
    double dSwitchVal = 2.0 * (distSq / forcefield.rCutSq - 1.0) * 2.0 *
                        dist / forcefield.rCutSq;
    return -qi_qj * (dSwitchVal / distSq - switchVal / (distSq * dist));
  }
}

inline double FF_SWITCH::CalcdEndL(const double distSq, const uint kind1,
                                  const uint kind2, const double lambda) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  double sigma6 = sigmaSq[index] * sigmaSq[index] * sigmaSq[index];
  double dist6 = distSq * distSq * distSq;
  double lambdaCoef = 0.5 * (1.0 - lambda) * (1.0 - lambda);
  double softDist6 = lambdaCoef * sigma6 + dist6;
  double softRsq = pow(softDist6, 1.0/3.0);
  double fCoef = lambda * (1.0 - lambda) * sigma6 / (6.0 * softRsq * softRsq);
  double dhdl = CalcEn(softRsq, index) + fCoef * CalcVir(softRsq, index);
  return dhdl;
}

//Calculate the dE/dlambda for Coulomb energy
inline double FF_SWITCH::CalcCoulombdEndL(const double distSq,
                                          const double qi_qj_Fact,
                                          const double lambda, uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  //Use amber scheme for soft core potential in Coulomb interaction
  // 12 is default value for Beta according to amber tutorial (Eq. 21.7)
  // http://ambermd.org/doc12/Amber18.pdf

  double lambdaCoef = 12.00 * (1.0 - lambda);
  double softRsq = lambdaCoef + distSq;
  //dE/dlambda = E + 6.0 * lambda + vir
  double dhdl = CalcCoulomb(softRsq, qi_qj_Fact, b) + 
                6.0 * lambda * CalcCoulombVir(softRsq, qi_qj_Fact, b);
  return dhdl; 
}

#endif /*FF_SWITCH_H*/
