/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
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

  double *shiftConst, *shiftConst_1_4;

};

inline void FF_SHIFT::Init(ff_setup::Particle const& mie,
                          ff_setup::NBfix const& nbfix)
{
  //Initializ sigma and epsilon
  FFParticle::Init(mie, nbfix);
  uint size = num::Sq(count);
  //allocate memory 
  shiftConst = new double [size];
  shiftConst_1_4 = new double [size];
  //calculate shift constant
  for(uint i = 0; i < count; ++i) {
    for(uint j = 0; j < count; ++j) {
      uint idx = FlatIndex(i, j);
      double rRat2 = sigmaSq[idx] / forcefield.rCutSq;
      double rRat4 = rRat2 * rRat2;
      double attract = rRat4 * rRat2;
      //for 1-4 interaction
      double rRat2_1_4 = sigmaSq_1_4[idx] / forcefield.rCutSq;
      double rRat4_1_4 = rRat2_1_4 * rRat2_1_4;
      double attract_1_4 = rRat4_1_4 * rRat2_1_4;
      double repulse = pow(sqrt(rRat2), n[idx]);
      double repulse_1_4 = pow(sqrt(rRat2_1_4), n_1_4[idx]);

      shiftConst[idx] =  epsilon_cn[idx] * (repulse - attract);
      shiftConst_1_4[idx] =  epsilon_cn_1_4[idx] *
                            (repulse_1_4 - attract_1_4);
    
    }
  }
}
                  

inline void FF_SHIFT::CalcAdd_1_4(double& en, const double distSq,
                                  const uint kind1, const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
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

  en += (epsilon_cn_1_4[index] * (repulse - attract) - shiftConst_1_4[index]);
}

inline void FF_SHIFT::CalcCoulombAdd_1_4(double& en, const double distSq,
    const double qi_qj_Fact,
    const bool NB) const
{
  double dist = sqrt(distSq);
  if(NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * forcefield.scaling_14 / dist;
}

inline double FF_SHIFT::CalcEn(const double distSq, const uint kind1,
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

inline double FF_SHIFT::CalcEn(const double distSq, const uint index) const
{
  double rRat2 = sigmaSq[index] / distSq;
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double n_ij = n[index];
  double repulse = pow(rRat2, (n_ij * 0.5));

  return (epsilon_cn[index] * (repulse - attract) - shiftConst[index]);
}

inline double FF_SHIFT::CalcVir(const double distSq, const uint kind1, 
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

inline double FF_SHIFT::CalcVir(const double distSq, const uint index) const
{
  double rNeg2 = 1.0 / distSq;
  double rRat2 = rNeg2 * sigmaSq[index];
  double rRat4 = rRat2 * rRat2;
  double attract = rRat4 * rRat2;
  double n_ij = n[index];
  double repulse = pow(rRat2, (n_ij * 0.5));

  //Virial is the derivative of the pressure... mu
  return epsilon_cn_6[index] * (nOver6[index] * repulse - attract) * rNeg2;
}

inline double FF_SHIFT::CalcCoulomb(const double distSq,
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

inline double FF_SHIFT::CalcCoulomb(const double distSq, const double qi_qj_Fact, 
                                    const uint b) const
{
  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double val = forcefield.alpha[b] * dist;
    return qi_qj_Fact * erfc(val) / dist;
  } else {
    double dist = sqrt(distSq);
    return qi_qj_Fact * (1.0 / dist - 1.0 / forcefield.rCut);
  }
}

inline double FF_SHIFT::CalcCoulombVir(const double distSq, const double qi_qj,
                                      const double lambda, const uint b) const
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

inline double FF_SHIFT::CalcCoulombVir(const double distSq, const double qi_qj,
                                       uint b) const
{ 
  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double constValue = 2.0 * forcefield.alpha[b] / sqrt(M_PI);
    double expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    double temp = erfc(forcefield.alpha[b] * dist);
    return qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    double dist = sqrt(distSq);
    return qi_qj / (distSq * dist);
  }
}

inline double FF_SHIFT::CalcdEndL(const double distSq, const uint kind1,
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
inline double FF_SHIFT::CalcCoulombdEndL(const double distSq,
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

#endif /*FF_SHIFT_H*/
