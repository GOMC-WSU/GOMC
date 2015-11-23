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
  class FFValues; 
  class FFKind;
}

struct FFParticle
{
 public:

   FFParticle();
   ~FFParticle(void);

   double GetMass(const uint kind) const { return mass[kind]; }
   virtual void CalcAdd(double& en, double& vir, const double distSq,
                const uint kind1, const uint kind2) const;
   virtual void CalcSub(double& en, double& vir, const double distSq,
                const uint kind1, const uint kind2) const;
   virtual double CalcEn(const double distSq,
                 const uint kind1, const uint kind2) const;
   virtual double CalcVir(const double distSq,
                  const uint kind1, const uint kind2) const;
   // Not Implemented
   void CalcAdd_1_4(double& en, double& vir, const double distSq,
		const uint kind1, const uint kind2) const;
   void CalcSub_1_4(double& en, double& vir, const double distSq,
		const uint kind1, const uint kind2) const;
   double CalcEn_1_4(const double distSq, 
		 const uint kind1, const uint kind2) const;
   double CalcVir_1_4(const double distSq,
		  const uint kind1, const uint kind2) const;
   //Not Implemented

  void Init(ff_setup::Particle const& mie,
	    ff_setup::NBfix const& nbfix,
	    config_setup::FFValues const& val,
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
   // Not Implemented
   void Calc_1_4(double& en, double& vir, const double distSq, uint index,
#ifdef MIE_INT_ONLY
		 const uint n
#else
		 const double n
#endif
		 ) const;
   // Not Implemented



  uint FlatIndex(const uint i, const uint j) const { return i + j * count; }
  void Blend(ff_setup::Particle const& mie, const double rCut);
  void AdjNBfix(ff_setup::Particle const& mie, ff_setup::NBfix const& nbfix,
		const double rCut);

   //vars for lj particles.
   double* mass;
  //char * name;
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
     * enCorrection, * virCorrection, *shiftConst, *An, *Bn, *Cn, *sig6, *sign;

  double rCut, rCutSq, rOn, rOnSq, A6, B6, C6, factor1, factor2;

  uint count, vdwKind;
  bool isMartini;
};

inline void FFParticle::CalcAdd(double& en, double& vir, const double distSq,
				const uint kind1, const uint kind2) const
{
   uint idx = FlatIndex(kind1, kind2);
   Calc(en, vir, distSq, idx, n[idx]);
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

//mie potential
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

#endif /*FF_PARTICLE_H*/
