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


struct FF_SWITCH : public FFParticle
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

inline void FF_SWITCH::CalcAdd(double& en, double& vir, const double distSq,
				const uint kind1, const uint kind2) const
{
   uint idx = FlatIndex(kind1, kind2);
   Calc(en, vir, distSq, idx, n[idx]);
} 

inline void FF_SWITCH::CalcCoulombAdd(double& en, double& vir,
					const double distSq,
					const double qi_qj_Fact) const
{
  CalcCoulomb(en, vir, distSq, qi_qj_Fact);
}

inline void FF_SWITCH::CalcAdd_1_4(double& en, const double distSq,
		const uint kind1, const uint kind2) const
{
   uint index = FlatIndex(kind1, kind2);
   double rCutSq_rijSq = rCutSq - distSq;
   double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

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

   double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

   const double factE = ( distSq > rOnSq ? fE : 1.0);

   en += (epsilon_cn_1_4[index] * (repulse-attract)) * factE;
}

inline void FF_SWITCH::CalcCoulombAdd_1_4(double& en, const double distSq,
					    const double qi_qj_Fact) const
{
   double dist = sqrt(distSq);
   double switchVal = distSq/rCutSq - 1.0;
   switchVal *= switchVal;
   en += scaling_14 * qi_qj_Fact * switchVal/dist; 
}

inline void FF_SWITCH::CalcSub(double& en, double& vir, const double distSq,
				const uint kind1, const uint kind2) const
{
   double tempEn=0, tempVir=0;
   uint idx = FlatIndex(kind1, kind2);
   Calc(tempEn, tempVir, distSq, idx, n[idx]);
   en -= tempEn;
   vir = -1.0 * tempVir;
} 

inline void FF_SWITCH::CalcCoulombSub(double& en, double& vir,
					const double distSq,
					const double qi_qj_Fact) const
{
  double tempEn = 0.0, tempVir = 0.0;
  CalcCoulomb(tempEn, tempVir, distSq, qi_qj_Fact);
  en  -= tempEn;
  vir -= tempVir;
}

//mie potential
inline double FF_SWITCH::CalcEn(const double distSq,
                                 const uint kind1, const uint kind2) const
{
   uint index = FlatIndex(kind1, kind2);
   
   double rCutSq_rijSq = rCutSq - distSq;
   double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

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

   double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);

   const double factE = ( distSq > rOnSq ? fE : 1.0);

   return (epsilon_cn[index] * (repulse-attract)) * factE;
}

inline double FF_SWITCH::CalcCoulombEn(const double distSq,
					const double qi_qj_Fact) const
{
   double dist = sqrt(distSq);
   double switchVal = distSq/rCutSq - 1.0;
   switchVal *= switchVal;
   return  qi_qj_Fact * switchVal/dist;
}

//mie potential
inline double FF_SWITCH::CalcVir(const double distSq,
                                  const uint kind1, const uint kind2) const
{
   uint index = FlatIndex(kind1, kind2);

   double rCutSq_rijSq = rCutSq - distSq;
   double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

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

   double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);
   double fW = 12.0 * factor2 * rCutSq_rijSq * (rOnSq - distSq);

   const double factE = ( distSq > rOnSq ? fE : 1.0);
   const double factW = ( distSq > rOnSq ? fW : 0.0);

   double Wij = epsilon_cn_6[index] * (nOver6[index]*repulse-attract)*rNeg2;
   double Eij = epsilon_cn[index] * (repulse-attract);
   
   return (Wij * factE - Eij * factW);
}

inline double FF_SWITCH::CalcCoulombVir(const double distSq,
					 const double qi_qj_Fact) const
{  
   double dist = sqrt(distSq);
   double switchVal = distSq/rCutSq - 1.0;
   switchVal *= switchVal;
   double dSwitchVal = 2.0 * (distSq/rCutSq - 1.0) * 2.0 * dist/rCutSq;
   return -1.0 * qi_qj_Fact * (dSwitchVal/distSq - switchVal/(distSq * dist));
}

//mie potential
inline void FF_SWITCH::Calc(double & en, double & vir, 
			     const double distSq, const uint index,
#ifdef MIE_INT_ONLY
			     const uint n,
#else
			     const double n
#endif
			     ) const
{
   double rCutSq_rijSq = rCutSq - distSq;
   double rCutSq_rijSq_Sq = rCutSq_rijSq * rCutSq_rijSq;

   double rNeg2 = 1.0/distSq;
   double rRat2 = rNeg2 * sigmaSq[index];
   double rRat4 = rRat2 * rRat2;
   double attract = rRat4 * rRat2;
#ifdef MIE_INT_ONLY
   double repulse = num::POW(rRat2, rRat4, attract, n);
#else
   double repulse = pow(sqrt(rRat2), n);
#endif

   double fE = rCutSq_rijSq_Sq * factor2 * (factor1 + 2 * distSq);
   double fW = 12.0 * factor2 * rCutSq_rijSq * (rOnSq - distSq);

   const double factE = ( distSq > rOnSq ? fE : 1.0);
   const double factW = ( distSq > rOnSq ? fW : 0.0);

   double Wij = epsilon_cn_6[index] * (nOver6[index]*repulse-attract)*rNeg2;
   double Eij = epsilon_cn[index] * (repulse-attract);
   
   en += Eij * factE;
   vir = Wij * factE - Eij * factW;  
}

inline void FF_SWITCH::CalcCoulomb(double & en, double & vir,
				    const double distSq, 
				    const double qi_qj_Fact)const
{
   double dist = sqrt(distSq);
   double switchVal = distSq/rCutSq - 1.0;
   switchVal *= switchVal;
   double dSwitchVal = 2.0 * (distSq/rCutSq - 1.0) * 2.0 * dist/rCutSq;

   en += qi_qj_Fact * switchVal/dist;
   vir = -1.0 * qi_qj_Fact * (dSwitchVal/distSq - switchVal/(distSq * dist));
 
}

#endif /*FF_SWITCH_H*/
