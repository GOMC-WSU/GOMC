/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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

  virtual real CalcEn(const real distSq,
                        const uint kind1, const uint kind2) const;
  virtual real CalcVir(const real distSq,
                         const uint kind1, const uint kind2) const;
  virtual void CalcAdd_1_4(real& en, const real distSq,
                           const uint kind1, const uint kind2) const;

  // coulomb interaction functions
  virtual real CalcCoulomb(const real distSq,
                             const real qi_qj_Fact, const uint b) const;
  virtual real CalcCoulombVir(const real distSq,
                                const real qi_qj, const uint b) const;
  virtual void CalcCoulombAdd_1_4(real& en, const real distSq,
                                  const real qi_qj_Fact,
                                  const bool NB) const;

  //!Returns Ezero, no energy correction
  virtual real EnergyLRC(const uint kind1, const uint kind2) const
  {
    return 0.0;
  }
  //!!Returns Ezero, no virial correction
  virtual real VirialLRC(const uint kind1, const uint kind2) const
  {
    return 0.0;
  }

protected:

  real rOn, rOnSq, factor1, factor2;

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

inline void FF_SWITCH::CalcAdd_1_4(real& en, const real distSq,
                                   const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  real rCutSq_rijSq = forcefield.rCutSq - distSq;
  real rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  real rRat2 = sigmaSq_1_4[index] / distSq;
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n_1_4[index];
  real repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  real n_ij = n_1_4[index];
  real repulse = pow(sqrt(rRat2), n_ij);
#endif

  real fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

  const real factE = ( distSq > rOnSq ? fE : 1.0);

  en += (epsilon_cn_1_4[index] * (repulse - attract)) * factE;
}

inline void FF_SWITCH::CalcCoulombAdd_1_4(real& en, const real distSq,
    const real qi_qj_Fact,
    const bool NB) const
{
  real dist = sqrt(distSq);
  if(NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * forcefield.scaling_14 / dist;
}


//mie potential
inline real FF_SWITCH::CalcEn(const real distSq,
                                const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  real rCutSq_rijSq = forcefield.rCutSq - distSq;
  real rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  real rRat2 = sigmaSq[index] / distSq;
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n[index];
  real repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  real n_ij = n[index];
  real repulse = pow(sqrt(rRat2), n_ij);
#endif

  real fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

  const real factE = ( distSq > rOnSq ? fE : 1.0);

  return (epsilon_cn[index] * (repulse - attract)) * factE;
}

inline real FF_SWITCH::CalcCoulomb(const real distSq,
                                     const real qi_qj_Fact, const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if(forcefield.ewald) {
    real dist = sqrt(distSq);
    real val = forcefield.alpha[b] * dist;
    return  qi_qj_Fact * erfc(val) / dist;
  } else {
    real dist = sqrt(distSq);
    real switchVal = distSq / forcefield.rCutSq - 1.0;
    switchVal *= switchVal;
    return  qi_qj_Fact * switchVal / dist;
  }
}

//mie potential
inline real FF_SWITCH::CalcVir(const real distSq,
                                 const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
  real rCutSq_rijSq = forcefield.rCutSq - distSq;
  real rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

  real rNeg2 = 1.0 / distSq;
  real rRat2 = rNeg2 * sigmaSq[index];
  real rRat4 = rRat2 * rRat2;
  real attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
  uint n_ij = n[index];
  real repulse = num::POW(rRat2, rRat4, attract, n_ij);
#else
  real n_ij = n[index];
  real repulse = pow(sqrt(rRat2), n_ij);
#endif

  real fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);
  real fW = 12.0 * factor2 * rCutSq_rijSq * (rOnSq - distSq);

  const real factE = ( distSq > rOnSq ? fE : 1.0);
  const real factW = ( distSq > rOnSq ? fW : 0.0);

  real Wij = epsilon_cn_6[index] * (nOver6[index] * repulse - attract) * rNeg2;
  real Eij = epsilon_cn[index] * (repulse - attract);

  return (Wij * factE - Eij * factW);
}

inline real FF_SWITCH::CalcCoulombVir(const real distSq,
                                        const real qi_qj, const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if(forcefield.ewald) {
    real dist = sqrt(distSq);
    real constValue = 2.0 * forcefield.alpha[b] / sqrt(M_PI);
    real expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    real temp = erfc(forcefield.alpha[b] * dist);
    return  qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    real dist = sqrt(distSq);
    real switchVal = distSq / forcefield.rCutSq - 1.0;
    switchVal *= switchVal;
    real dSwitchVal = 2.0 * (distSq / forcefield.rCutSq - 1.0) * 2.0 *
                        dist / forcefield.rCutSq;
    return -1.0 * qi_qj * (dSwitchVal / distSq - switchVal / (distSq * dist));
  }
}


#endif /*FF_SWITCH_H*/
