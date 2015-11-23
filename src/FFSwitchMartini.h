#ifndef FF_SWITCH_MARTINI_H
#define FF_SWITCH_MARTINI_H

#include "EnsemblePreprocessor.h" //For "MIE_INT_ONLY" preprocessor.
#include "FFConst.h" //constants related to particles.
#include "../lib/BasicTypes.h" //for uint
#include "../lib/NumLib.h" //For Cb, Sq
#include "FFParticle.h"


///////////////////////////////////////////////////////////////////////
////////////////////////// LJ Switch Style ////////////////////////////
///////////////////////////////////////////////////////////////////////
// Virial and LJ potential calculation:
// Eij = cn * eps_ij * ( (sig_ij/rij)^n - (sig_ij/rij)^6)
// cn = n/(n-6) * ((n/6)^(6/(n-6)))
//

struct FF_SWITCH_MARTINI : public FFParticle
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

};

inline void FF_SWITCH_MARTINI::CalcAdd(double& en, double& vir, const double distSq,
				const uint kind1, const uint kind2) const
{
   uint idx = FlatIndex(kind1, kind2);
   Calc(en, vir, distSq, idx, n[idx]);
} 

inline void FF_SWITCH_MARTINI::CalcSub(double& en, double& vir, const double distSq,
				const uint kind1, const uint kind2) const
{
   double tempEn=0, tempVir=0;
   uint idx = FlatIndex(kind1, kind2);
   Calc(tempEn, tempVir, distSq, idx, n[idx]);
   en -= tempEn;
   vir = -1.0 * tempVir;
} 

//mie potential
inline double FF_SWITCH_MARTINI::CalcEn(const double distSq,
                                 const uint kind1, const uint kind2) const
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

//mie potential
inline double FF_SWITCH_MARTINI::CalcVir(const double distSq,
                                  const uint kind1, const uint kind2) const
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


//mie potential
inline void FF_SWITCH_MARTINI::Calc(double & en, double & vir, 
			     const double distSq, const uint index,
#ifdef MIE_INT_ONLY
			     const uint n,
#else
			     const double n
#endif
			     ) const
{

   double r_1 = 1/sqrt(distSq);
   double r_2 = 1.0/distSq;
   double r_4 = r_2 * r_2;
   double r_6 = r_4 * r_2;  
   double r_8 = r_4 * r_4;
   double r_n = pow(r_1, n);
   double r_n2 = pow(r_1, n + 2);

   double rij_ron = sqrt(distSq) - rOn;
   double rij_ron_2 = rij_ron * rij_ron;
   double rij_ron_3 = rij_ron_2 * rij_ron;
   double rij_ron_4 = rij_ron_2 * rij_ron_2;

   //energy
   double shifttempRep = -(An[index]/3.0)*rij_ron_3 -
     (Bn[index]/4.0)*rij_ron_4 - Cn[index];
   double shifttempAtt = -(A6/3.0)*rij_ron_3 - (B6/4.0)*rij_ron_4 - C6;
   
   const double shiftRep = ( distSq > rOnSq ? shifttempRep : -Cn[index]);
   const double shiftAtt = ( distSq > rOnSq ? shifttempAtt : -C6);
   
   en += epsilon_cn[index] * (sign[index] * (r_n + shiftRep) - 
				     sig6[index] * (r_6 + shiftAtt));

   //virial
   double dshifttempRep = An[index] * rij_ron_2 + Bn[index] * rij_ron_3;
   double dshifttempAtt = A6 * rij_ron_2 + B6 * rij_ron_3;

   const double dshiftRep = ( distSq > rOnSq ? dshifttempRep * r_1 : 0);
   const double dshiftAtt = ( distSq > rOnSq ? dshifttempAtt * r_1 : 0);
   
   vir = epsilon_cn[index] * (sign[index] *
				     (n * r_n2 + dshiftRep) - 
				     sig6[index] * (6.0 * r_8 + dshiftAtt));
   
}

#endif /*FF_SWITCH_MARTINI_H*/
