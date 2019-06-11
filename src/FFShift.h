/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
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


struct FF_SHIFT : public FFParticle {
public:
  FF_SHIFT(Forcefield &ff) : FFParticle(ff), shiftConst(NULL), shiftConst_1_4(NULL) {}
  virtual ~FF_SHIFT()
  {
    delete[] shiftConst;
    delete[] shiftConst_1_4;
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

  real *shiftConst, *shiftConst_1_4;

};

inline void FF_SHIFT::Init(ff_setup::Particle const& mie,
                           ff_setup::NBfix const& nbfix)
{
  //Initializ sigma and epsilon
  FFParticle::Init(mie, nbfix);
  uint size = num::Sq(count);
  //allocate memory
  shiftConst = new real [size];
  shiftConst_1_4 = new real [size];
  //calculate shift constant
  for(uint i = 0; i < count; ++i) {
    for(uint j = 0; j < count; ++j) {
      uint idx = FlatIndex(i, j);
      real rRat2 = sigmaSq[idx] / forcefield.rCutSq;
      real rRat4 = rRat2 * rRat2;
      real attract = rRat4 * rRat2;
      //for 1-4 interaction
      real rRat2_1_4 = sigmaSq_1_4[idx] / forcefield.rCutSq;
      real rRat4_1_4 = rRat2_1_4 * rRat2_1_4;
      real attract_1_4 = rRat4_1_4 * rRat2_1_4;
      real repulse = pow(sqrt(rRat2), n[idx]);
      real repulse_1_4 = pow(sqrt(rRat2_1_4), n_1_4[idx]);

      shiftConst[idx] =  epsilon_cn[idx] * (repulse - attract);
      shiftConst_1_4[idx] =  epsilon_cn_1_4[idx] *
                             (repulse_1_4 - attract_1_4);

    }
  }
}


inline void FF_SHIFT::CalcAdd_1_4(real& en, const real distSq,
                                  const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
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

  en += (epsilon_cn_1_4[index] * (repulse - attract) - shiftConst_1_4[index]);
}

inline void FF_SHIFT::CalcCoulombAdd_1_4(real& en, const real distSq,
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
inline real FF_SHIFT::CalcEn(const real distSq,
                               const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
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

  return (epsilon_cn[index] * (repulse - attract) - shiftConst[index]);
}

inline real FF_SHIFT::CalcCoulomb(const real distSq,
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
    return  qi_qj_Fact * (1.0 / dist - 1.0 / forcefield.rCut);
  }
}

//mie potential
inline real FF_SHIFT::CalcVir(const real distSq,
                                const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint index = FlatIndex(kind1, kind2);
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

  //Virial is the derivative of the pressure... mu
  return epsilon_cn_6[index] * (nOver6[index] * repulse - attract) * rNeg2;
}

inline real FF_SHIFT::CalcCoulombVir(const real distSq,
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
    return qi_qj / (distSq * dist);
  }
}


#endif /*FF_SHIFT_H*/
