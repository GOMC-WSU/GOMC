/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.9
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_SWITCH_MARTINI_H
#define FF_SWITCH_MARTINI_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "BasicTypes.h" //for uint
#include "NumLib.h" //For Cb, Sq
#include "FFParticle.h"


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


struct FF_SWITCH_MARTINI : public FFParticle
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
                                const double qi_qj) const;
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


inline void FF_SWITCH_MARTINI::CalcAdd_1_4(double& en, const double distSq,
    const uint kind1,
    const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  double r_2 = 1.0/distSq;
  double r_4 = r_2 * r_2;
  double r_6 = r_4 * r_2;
#ifdef MIE_INT_ONLY
  uint n_ij = n_1_4[index];
  double r_n = num::POW(r_2, r_4, attract, n_ij);
#else
  double n_ij = n_1_4[index];
  double r_n = pow(sqrt(r_2), n_ij);
#endif

  double rij_ron = sqrt(distSq) - rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double shifttempRep = -(An_1_4[index]/3.0)*rij_ron_3 -
                        (Bn_1_4[index]/4.0)*rij_ron_4 - Cn_1_4[index];
  double shifttempAtt = -(A6/3.0)*rij_ron_3 - (B6/4.0)*rij_ron_4 - C6;

  const double shiftRep = ( distSq > rOnSq ? shifttempRep : -Cn_1_4[index]);
  const double shiftAtt = ( distSq > rOnSq ? shifttempAtt : -C6);

  en += epsilon_cn_1_4[index] * (sign_1_4[index] * (r_n + shiftRep) -
                                 sig6_1_4[index] * (r_6 + shiftAtt));
}

inline void FF_SWITCH_MARTINI::CalcCoulombAdd_1_4(double& en,
    const double distSq,
    const double qi_qj_Fact) const
{
  if(ewald)
  {
     double dist = sqrt(distSq);
     double erfc = alpha * dist;
     en += scaling_14 * qi_qj_Fact * (1 - erf(erfc))/ dist;
  }
  else
  {
     // in Martini, the Coulomb switching distance is zero, so we will have
     // sqrt(distSq) - rOnCoul =  sqrt(distSq)
     double dist = sqrt(distSq);
     double rij_ronCoul_3 = dist * distSq;
     double rij_ronCoul_4 = distSq * distSq;
     
     double coul = -(A1/3.0) * rij_ronCoul_3 - (B1/4.0) * rij_ronCoul_4 - C1;
     en += scaling_14 * qi_qj_Fact * diElectric_1 * (coul + 1.0/dist);
  }
}


//mie potential
inline double FF_SWITCH_MARTINI::CalcEn(const double distSq,
                                        const uint kind1,
					const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);

  double r_2 = 1.0/distSq;
  double r_4 = r_2 * r_2;
  double r_6 = r_4 * r_2;
#ifdef MIE_INT_ONLY
  uint n_ij = n[index];
  double r_n = num::POW(r_2, r_4, attract, n_ij);
#else
  double n_ij = n[index];
  double r_n = pow(sqrt(r_2), n_ij);
#endif

  double rij_ron = sqrt(distSq) - rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;
  double rij_ron_4 = rij_ron_2 * rij_ron_2;

  double shifttempRep = -(An[index]/3.0)*rij_ron_3 -
                        (Bn[index]/4.0)*rij_ron_4 - Cn[index];
  double shifttempAtt = -(A6/3.0)*rij_ron_3 - (B6/4.0)*rij_ron_4 - C6;

  const double shiftRep = ( distSq > rOnSq ? shifttempRep : -Cn[index]);
  const double shiftAtt = ( distSq > rOnSq ? shifttempAtt : -C6);

  double Eij = epsilon_cn[index] * (sign[index] * (r_n + shiftRep) -
                                    sig6[index] * (r_6 + shiftAtt));
  return Eij;
}

inline double FF_SWITCH_MARTINI::CalcCoulombEn(const double distSq,
					       const double qi_qj_Fact) const
{
  if(distSq <= rCutLowSq)
    return num::BIGNUM;

  if(ewald)
  {
     double dist = sqrt(distSq);
     double erfc = alpha * dist;
     return  qi_qj_Fact * (1 - erf(erfc))/ dist;
  }
  else
  {     
     // in Martini, the Coulomb switching distance is zero, so we will have
     // sqrt(distSq) - rOnCoul =  sqrt(distSq)
     double dist = sqrt(distSq);
     double rij_ronCoul_3 = dist * distSq;
     double rij_ronCoul_4 = distSq * distSq;
     
     double coul = -(A1/3.0) * rij_ronCoul_3 - (B1/4.0) * rij_ronCoul_4 - C1;
     return qi_qj_Fact  * diElectric_1 * (1.0/dist + coul);
  }
}

//mie potential
inline double FF_SWITCH_MARTINI::CalcVir(const double distSq,
					 const uint kind1,
					 const uint kind2) const
{
  uint index = FlatIndex(kind1, kind2);
  double n_ij = n[index];

  double r_1 = 1.0/sqrt(distSq);
  double r_8 = pow(r_1, 8);
  double r_n2 = pow(r_1, n_ij + 2);

  double rij_ron = sqrt(distSq) - rOn;
  double rij_ron_2 = rij_ron * rij_ron;
  double rij_ron_3 = rij_ron_2 * rij_ron;


  double dshifttempRep = An[index] * rij_ron_2 + Bn[index] * rij_ron_3;
  double dshifttempAtt = A6 * rij_ron_2 + B6 * rij_ron_3;

  const double dshiftRep = ( distSq > rOnSq ? dshifttempRep * r_1 : 0);
  const double dshiftAtt = ( distSq > rOnSq ? dshifttempAtt * r_1 : 0);

  double Wij = epsilon_cn[index] * (sign[index] *
                                    (n_ij * r_n2 + dshiftRep) -
                                    sig6[index] * (6.0 * r_8 + dshiftAtt));
  return Wij;

}

inline double FF_SWITCH_MARTINI::CalcCoulombVir(const double distSq,
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
     // in Martini, the Coulomb switching distance is zero, so we will have
     // sqrt(distSq) - rOnCoul =  sqrt(distSq)
     double dist = sqrt(distSq);
     double rij_ronCoul_2 = distSq;
     double rij_ronCoul_3 = dist * distSq;
     double rij_ronCoul_4 = distSq * distSq;

     double virCoul = A1/rij_ronCoul_2 + B1/rij_ronCoul_3;
     return qi_qj * diElectric_1 * ( 1.0/(dist * distSq) + virCoul/dist);
  }
}


#endif /*FF_SWITCH_MARTINI_H*/
