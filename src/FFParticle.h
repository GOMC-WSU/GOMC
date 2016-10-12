/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.8
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FF_PARTICLE_H
#define FF_PARTICLE_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "../lib/BasicTypes.h" //for uint
#include "../lib/NumLib.h" //For Cb, Sq

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

   double GetMass(const uint kind) const { return mass[kind]; }
  // LJ interaction functions
   virtual void CalcAdd(double& en, double& vir, const double distSq,
                const uint kind1, const uint kind2) const;
   virtual void CalcSub(double& en, double& vir, const double distSq,
                const uint kind1, const uint kind2) const;
   virtual double CalcEn(const double distSq,
                 const uint kind1, const uint kind2) const;
   virtual double CalcVir(const double distSq,
                  const uint kind1, const uint kind2) const;
   virtual void CalcAdd_1_4(double& en, const double distSq,
		const uint kind1, const uint kind2) const;

  // coulomb interaction functions
   virtual void CalcCoulombAdd(double& en, double& vir, const double distSq,
			       const double qi_qj_Fact) const;
   virtual void CalcCoulombSub(double& en, double& vir, const double distSq,
			       const double qi_qj_Fact) const;
   virtual double CalcCoulombEn(const double distSq,
				const double qi_qj_Fact) const;
   virtual double CalcCoulombVir(const double distSq,
				 const double qi_qj_Fact) const;
   virtual void CalcCoulombAdd_1_4(double& en, const double distSq,
				   const double qi_qj_Fact) const;

  void Init(ff_setup::Particle const& mie,
	    ff_setup::NBfix const& nbfix,
	    config_setup::SystemVals const& sys,
	    config_setup::FFKind const& ffKind);
   //!Returns Energy long-range correction term for a kind pair
   virtual double EnergyLRC(const uint kind1, const uint kind2) const;
   //!Returns Energy long-range correction term for a kind pair
   virtual double VirialLRC(const uint kind1, const uint kind2) const;

   uint NumKinds() const { return count; }

 protected:

   virtual void Calc(double& en, double& vir, const double distSq, uint index,
#ifdef MIE_INT_ONLY
	     const uint n
#else
	     const double n
#endif
	     ) const;

  virtual void CalcCoulomb(double& en, double& vir, const double distSq,
			 const double qi_qj_Fact)const;




  uint FlatIndex(const uint i, const uint j) const { return i + j * count; }
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

  uint count, vdwKind;
  bool isMartini;
};

inline void FFParticle::CalcAdd(double& en, double& vir, const double distSq,
				const uint kind1, const uint kind2) const
{
   uint idx = FlatIndex(kind1, kind2);
   Calc(en, vir, distSq, idx, n[idx]);
}

inline void FFParticle::CalcCoulombAdd(double& en, double& vir,
					const double distSq,
					const double qi_qj_Fact) const
{
  CalcCoulomb(en, vir, distSq, qi_qj_Fact);
}


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
					    const double qi_qj_Fact) const
{
   double dist = sqrt(distSq);
   double erfc = alpha * dist;
   en += scaling_14 * qi_qj_Fact * (1 - erf(erfc))/ dist;
}

inline void FFParticle::CalcSub(double& en, double& vir, const double distSq,
				const uint kind1, const uint kind2) const
{
   double tempEn=0, tempVir=0;
   uint idx = FlatIndex(kind1, kind2);
   Calc(tempEn, tempVir, distSq, idx, n[idx]);
   en -= tempEn;
   vir = -1.0 * tempVir;
}

inline void FFParticle::CalcCoulombSub(double& en, double& vir,
					const double distSq,
					const double qi_qj_Fact) const
{
  double tempEn = 0.0, tempVir = 0.0;
  CalcCoulomb(tempEn, tempVir, distSq, qi_qj_Fact);
  en  -= tempEn;
  vir = -1.0 * tempVir;
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

inline double FFParticle::CalcCoulombEn(const double distSq,
					const double qi_qj_Fact) const
{
   double dist = sqrt(distSq);
   double erfc = alpha * dist;
   return  qi_qj_Fact * (1 - erf(erfc))/ dist;
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
					 const double qi_qj_Fact) const
{
  // need to figure out -d(erf)/dr
  //double dist = sqrt(distSq);
  //double erfc = alpha / boxSize * dist;
  //return  qi_qj_Fact * (1 - erf(erfc))/ (distSq * dist);
  return 0.0;
}

//mie potential
inline void FFParticle::Calc(double & en, double & vir,
			     const double distSq, const uint index,
#ifdef MIE_INT_ONLY
			     const uint n,
#else
			     const double n
#endif
			     ) const
{
   double rNeg2 = 1.0/distSq;
   double rRat2 = rNeg2 * sigmaSq[index];
   double rRat4 = rRat2 * rRat2;
   double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
   double repulse = num::POW(rRat2, rRat4, attract, n);
#else
   double repulse = pow(sqrt(rRat2), n);
#endif

   en += epsilon_cn[index] * (repulse-attract);
   //Virial is the derivative of the pressure... mu
   vir = epsilon_cn_6[index] * (nOver6[index]*repulse-attract)*rNeg2;
}

inline void FFParticle::CalcCoulomb(double & en, double & vir,
				    const double distSq,
				    const double qi_qj_Fact)const
{
   double dist = sqrt(distSq);
   double erfc = alpha * dist;
   en += qi_qj_Fact * (1 - erf(erfc))/ dist;
   // need to figure out -d(erf)/dr
   //vir = qi_qj_Fact * (1 - erf(erfc))/(distSq * dist)
}


#endif /*FF_PARTICLE_H*/
