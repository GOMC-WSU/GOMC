/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_EXP6_H
#define FF_EXP6_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "BasicTypes.h" //for uint
#include "NumLib.h" //For Cb, Sq
#include "FFParticle.h"

//////////////////////////////////////////////////////////////////////
//////////////////////////// Exp-6 Style /////////////////////////////
//////////////////////////////////////////////////////////////////////
// Virial and LJ potential calculation:
// U(rij)= expConst * ( (6/alpha) * exp( alpha * [1-(r/rmin)] )-(rmin/r)^6) )
// expConst=( eps-ij * alpha )/(alpha-6)
//
// Vir(r)= -du/dr * r
// Vir(r)= 6 * expConst * [(r/rmin) * exp(alpha * [1-(r/rmin)])-(rmin/r)^6]/ r^2
//



struct FF_EXP6 : public FFParticle {
public:
  FF_EXP6(Forcefield &ff) : FFParticle(ff), expConst(NULL), rMin(NULL),
    rMaxSq(NULL), expConst_1_4(NULL), rMin_1_4(NULL), rMaxSq_1_4(NULL) {}
  virtual ~FF_EXP6()
  {
    delete[] expConst;
    delete[] expConst_1_4;
    delete[] rMin;
    delete[] rMin_1_4;
    delete[] rMaxSq;
    delete[] rMaxSq_1_4;
  }

  virtual void Init(ff_setup::Particle const& mie,
                    ff_setup::NBfix const& nbfix);

  virtual double GetRmin(const uint i, const uint j) const;
  virtual double GetRmax(const uint i, const uint j) const;
  virtual double GetRmin_1_4(const uint i, const uint j) const;
  virtual double GetRmax_1_4(const uint i, const uint j) const;

  virtual double CalcEn(const double distSq,
                        const uint kind1, const uint kind2) const;
  virtual double CalcVir(const double distSq,
                         const uint kind1, const uint kind2) const;
  virtual void CalcAdd_1_4(double& en, const double distSq,
                           const uint kind1, const uint kind2) const;

  // coulomb interaction functions
  virtual double CalcCoulomb(const double distSq,
                             const double qi_qj_Fact, const uint b) const;
  virtual double CalcCoulombVir(const double distSq,
                                const double qi_qj, const uint b) const;
  virtual void CalcCoulombAdd_1_4(double& en, const double distSq,
                                  const double qi_qj_Fact,
                                  const bool NB) const;

  //!Returns energy correction
  virtual double EnergyLRC(const uint kind1, const uint kind2) const
  {
    return enCorrection[FlatIndex(kind1, kind2)];
  }
  //!!Returns virial correction
  virtual double VirialLRC(const uint kind1, const uint kind2) const
  {
    return virCorrection[FlatIndex(kind1, kind2)];
  }

  protected:

  double *expConst, *expConst_1_4, *rMaxSq, *rMin, *rMaxSq_1_4, *rMin_1_4;

};

inline void FF_EXP6::Init(ff_setup::Particle const& mie,
                          ff_setup::NBfix const& nbfix)
{
  //Initializ sigma and epsilon
  FFParticle::Init(mie, nbfix);
  uint size = num::Sq(count);
  //allocate memory 
  expConst = new double [size];
  rMin = new double [size];
  rMaxSq = new double [size];
  expConst_1_4 = new double [size];
  rMin_1_4 = new double [size];
  rMaxSq_1_4 = new double [size];
  //calculate exp-6 parameter
  for(uint i = 0; i < count; ++i) {
    for(uint j = 0; j < count; ++j) {
      uint idx = FlatIndex(i, j);
      //We use n as alpha for exp-6, with geometric combining
      expConst[idx] = epsilon[idx] * n[idx] / (n[idx] - 6.0);
      expConst_1_4[idx] = epsilon_1_4[idx] * n_1_4[idx] / (n_1_4[idx] - 6.0);

      //Find the Rmin(well depth) 
      double sigma = sqrt(sigmaSq[idx]);
      num::Exp6Fun* func1 = new num::RminFun(n[idx], sigma);
      rMin[idx] = num::Zbrent(func1, sigma, 3.0 * sigma, 1.0e-7);
      //Find the Rmax(du/dr = 0) 
      num::Exp6Fun* func2 = new num::RmaxFun(n[idx], sigma, rMin[idx]);
      double rmax = Zbrent(func2, 0.0, sigma, 1.0e-7);
      rMaxSq[idx] = rmax * rmax;

      //Find the Rmin(well depth) 
      double sigma_1_4 = sqrt(sigmaSq_1_4[idx]);
      num::Exp6Fun* func3 = new num::RminFun(n_1_4[idx], sigma_1_4);
      rMin_1_4[idx] = num::Zbrent(func3, sigma_1_4, 3.0 * sigma_1_4, 1.0e-7);
      //Find the Rmax(du/dr = 0) 
      num::Exp6Fun* func4 = new num::RmaxFun(n_1_4[idx], sigma_1_4, rMin_1_4[idx]);
      double rmax_1_4 = Zbrent(func4, 0.0, sigma_1_4, 1.0e-7);
      rMaxSq_1_4[idx] = rmax_1_4 * rmax_1_4;
      //A,B and C for energy equation
      double A = 6.0 * epsilon[idx] * exp(n[idx]) / ((double)n[idx] - 6.0);
      double B = rMin[idx] / (double)n[idx];
      double C = epsilon[idx] * (double)n[idx] *  pow(rMin[idx], 6) / 
                  ((double)n[idx] - 6.0);

      enCorrection[idx] = 2.0 * M_PI * (A * B * exp(-rCut/B) *
                          (2.0 * pow(B,2) + (2.0 * B * rCut) + pow(rCut,2)) -
                          C / (3.0 * pow(rCut,3)));
      virCorrection[idx] = 2.0 * M_PI * (A * exp(-rCut/B) *
                           (6.0 * pow(B,3) + 6.0 * pow(B,2) * rCut + 
                           3.0 * pow(rCut,2) * B + pow(rCut,3)) -
                           2.0 * C / pow(rCut,3));

      delete func1;
      delete func2;
      delete func3;
      delete func4; 
    }
  }
}
                  

inline void FF_EXP6::CalcAdd_1_4(double& en, const double distSq,
                                  const uint kind1, const uint kind2) const
{
  uint idx = FlatIndex(kind1, kind2);
  if(distSq < rMaxSq_1_4[idx]) {
    en += num::BIGNUM;
    return;
  }

  double dist = sqrt(distSq);
  double rRat = rMin_1_4[idx] / dist;
  double rRat2 = rRat * rRat;
  double attract = rRat2 * rRat2 * rRat2;
  double alpha_ij = n_1_4[idx];
  double repulse = (6.0 / alpha_ij) * exp(alpha_ij * 
                  (1.0 - dist/rMin_1_4[idx]));

  en += expConst_1_4[idx] * (repulse - attract);
}

inline void FF_EXP6::CalcCoulombAdd_1_4(double& en, const double distSq,
    const double qi_qj_Fact,
    const bool NB) const
{
  double dist = sqrt(distSq);
  if(NB)
    en += qi_qj_Fact / dist;
  else
    en += qi_qj_Fact * forcefield.scaling_14 / dist;
}


//mie potential
inline double FF_EXP6::CalcEn(const double distSq,
                               const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint idx = FlatIndex(kind1, kind2);
  if(distSq < rMaxSq[idx]) {
    return num::BIGNUM;
  }

  double dist = sqrt(distSq);
  double rRat = rMin[idx] / dist;
  double rRat2 = rRat * rRat;
  double attract = rRat2 * rRat2 * rRat2;
  
  uint alpha_ij = n[idx];
  double repulse = (6.0 / alpha_ij) * exp(alpha_ij * 
                   (1.0 - dist/rMin[idx]));
  return expConst[idx] * (repulse - attract);
}

inline double FF_EXP6::CalcCoulomb(const double distSq,
                                    const double qi_qj_Fact, const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double val = forcefield.alpha[b] * dist;
    return  qi_qj_Fact * erfc(val) / dist;
  } else {
    double dist = sqrt(distSq);
    return  qi_qj_Fact / dist;
  }
}

//mie potential
inline double FF_EXP6::CalcVir(const double distSq,
                                const uint kind1, const uint kind2) const
{
  if(forcefield.rCutSq < distSq)
    return 0.0;

  uint idx = FlatIndex(kind1, kind2);

  if(distSq < rMaxSq[idx]) {
    return num::BIGNUM;
  }

  double dist = sqrt(distSq);
  double rRat = rMin[idx] / dist;
  double rRat2 = rRat * rRat ;
  double attract = rRat2 * rRat2 * rRat2;

  uint alpha_ij = n[idx];
  double repulse = (dist / rMin[idx]) * exp(alpha_ij * 
                   (1.0 - dist / rMin[idx]));

  //Virial = -r * du/dr
  return 6.0 * expConst[idx] * (repulse - attract) / distSq; 
}

inline double FF_EXP6::CalcCoulombVir(const double distSq,
                                       const double qi_qj, const uint b) const
{
  if(forcefield.rCutCoulombSq[b] < distSq)
    return 0.0;

  if(forcefield.ewald) {
    double dist = sqrt(distSq);
    double constValue = 2.0 * forcefield.alpha[b] / sqrt(M_PI);
    double expConstValue = exp(-1.0 * forcefield.alphaSq[b] * distSq);
    double temp = erfc(forcefield.alpha[b] * dist);
    return  qi_qj * (temp / dist + constValue * expConstValue) / distSq;
  } else {
    double dist = sqrt(distSq);
    return qi_qj / (distSq * dist);
  }
}

inline double FF_EXP6::GetRmin(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return rMin[idx];
}

inline double FF_EXP6::GetRmax(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return sqrt(rMaxSq[idx]);
}

inline double FF_EXP6::GetRmin_1_4(const uint i, const uint j) const 
{
  uint idx = FlatIndex(i, j);
  return rMin_1_4[idx];
}

inline double FF_EXP6::GetRmax_1_4(const uint i, const uint j) const
{
  uint idx = FlatIndex(i, j);
  return sqrt(rMaxSq_1_4[idx]);
}



#endif /*FF_EXP6_H*/
