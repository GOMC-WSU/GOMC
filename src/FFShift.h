#ifndef FF_SHIFT_H
#define FF_SHIFT_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "../lib/BasicTypes.h" //for uint
#include "../lib/NumLib.h" //For Cb, Sq
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

   //!Returns Ezero, no energy correction
   virtual double EnergyLRC(const uint kind1, const uint kind2) const
   {return 0.0;}
   //!!Returns Ezero, no virial correction
   virtual double VirialLRC(const uint kind1, const uint kind2) const
   {return 0.0;}

 private:
   virtual void Calc(double& en, double& vir, const double distSq, uint index,
#ifdef MIE_INT_ONLY
	     const uint n
#else
	     const double n
#endif
	     ) const;

   virtual void CalcCoulomb(double& en, double& vir, const double distSq,
			    const double qi_qj_Fact)const;
};

inline void FF_SHIFT::CalcAdd(double& en, double& vir, const double distSq,
			      const uint kind1, const uint kind2) const
{
   uint idx = FlatIndex(kind1, kind2);
   Calc(en, vir, distSq, idx, n[idx]);
} 

inline void FF_SHIFT::CalcCoulombAdd(double& en, double& vir,
				     const double distSq,
				     const double qi_qj_Fact) const
{
  CalcCoulomb(en, vir, distSq, qi_qj_Fact);
}

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
   en += scaling_14 * qi_qj_Fact * (1.0/dist - 1.0/rCut); 
}

inline void FF_SHIFT::CalcSub(double& en, double& vir, const double distSq,
			      const uint kind1, const uint kind2) const
{
   double tempEn=0, tempVir=0;
   uint idx = FlatIndex(kind1, kind2);
   Calc(tempEn, tempVir, distSq, idx, n[idx]);
   en -= tempEn;
   vir = -1.0 * tempVir;
} 

inline void FF_SHIFT::CalcCoulombSub(double& en, double& vir,
				     const double distSq,
				     const double qi_qj_Fact) const
{
  double tempEn = 0.0, tempVir = 0.0;
  CalcCoulomb(tempEn, tempVir, distSq, qi_qj_Fact);
  en  -= tempEn;
  vir -= tempVir;
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
   return  qi_qj_Fact * (1.0/dist - 1.0/rCut);
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
				       const double qi_qj_Fact) const
{  
  double dist = sqrt(distSq);
  return qi_qj_Fact/(distSq * dist);
}

//mie potential
inline void FF_SHIFT::Calc(double & en, double & vir, 
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

   en += (epsilon_cn[index] * (repulse-attract) - shiftConst[index]);
   //Virial is the derivative of the pressure... mu
   vir = epsilon_cn_6[index] * (nOver6[index]*repulse-attract)*rNeg2;
}

inline void FF_SHIFT::CalcCoulomb(double & en, double & vir,
				  const double distSq, 
				  const double qi_qj_Fact)const
{
   double dist = sqrt(distSq);
   en += qi_qj_Fact *(1.0/dist - 1.0/rCut);
   vir = qi_qj_Fact/(distSq * dist);
  
}


#endif /*FF_SHIFT_H*/
