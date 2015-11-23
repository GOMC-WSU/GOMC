#ifndef ENERGYTYPES_H
#define ENERGYTYPES_H

/*
 *    EnergyTypes.h
 *    Defines structs containing energy values
 *
 */

#include "EnsemblePreprocessor.h" //For # box const.
#include "../lib/BasicTypes.h" //For uint

#ifndef NDEBUG
#include <iostream>
#endif

#ifndef BOXES_WITH_U_NB
#if ENSEMBLE == GCMC || ENSEMBLE == NVT
#define BOXES_WITH_U_NB 1
#elif ENSEMBLE == GEMC
//case for NVT, GCMC
#define BOXES_WITH_U_NB 2
#endif
#endif

#ifndef BOXES_WITH_U_B
#if ENSEMBLE == NVT
#define BOXES_WITH_U_B 1
#elif ENSEMBLE == GEMC || ENSEMBLE == GCMC
//case for NVT, GCMC
#define BOXES_WITH_U_B 2
#endif
#endif

//long-range interactions between particles and all their infinite periodic images
struct Elect
{
	double virial, energy;

	Elect(): virial(0.0), energy(0.0) {}
	Elect(const double vir, const double en): virial(vir), energy(en) {}
	void Zero() {virial = 0;	energy = 0;}

	Elect& operator=(const Elect& rhs)
	{virial = rhs.virial; energy = rhs.energy; return *this;}
	Elect& operator+=(const Elect& rhs)
	{virial += rhs.virial; energy += rhs.energy; return *this;}
	Elect& operator-=(const Elect& rhs)
	{virial -= rhs.virial; energy -= rhs.energy; return *this;}
	Elect operator-(const Elect& rhs)
	{return Elect(virial - rhs.virial, energy - rhs.energy);}
	Elect operator+(const Elect& rhs)
	{return Elect(virial + rhs.virial, energy + rhs.energy);}
};

struct Intermolecular
{ 
   //MEMBERS
   double virial, energy;

   //CONSTRUCTORS
   Intermolecular() : virial(0.0), energy(0.0) {}
   Intermolecular(const double vir, const double en) : 
      virial(vir), energy(en) {} 

   //VALUE SETTER
   void Zero() { virial = energy = 0.0; }

   //OPERATORS
   Intermolecular& operator=(Intermolecular const& rhs) 
   { virial = rhs.virial; energy = rhs.energy; return *this; }
   Intermolecular& operator-=(Intermolecular const& rhs) 
   { virial -= rhs.virial; energy -= rhs.energy; return *this; }
   Intermolecular& operator+=(Intermolecular const& rhs) 
   { virial += rhs.virial; energy += rhs.energy; return *this; }
   Intermolecular& operator+=(Elect const& rhs) 
   { virial += rhs.virial; energy += rhs.energy; return *this; }
   Intermolecular& operator-=(Elect const& rhs) 
   { virial -= rhs.virial; energy -= rhs.energy; return *this; }
   Intermolecular operator-(Intermolecular const& rhs) 
   { return Intermolecular(virial - rhs.virial, energy - rhs.energy); }
   Intermolecular operator+(Intermolecular const& rhs) 
   { return Intermolecular(virial + rhs.virial, energy + rhs.energy); }
};




struct Energy
{
   //MEMBERS
   double intraBond, intraNonbond, inter, tc, total, real, recip, self, correction, elect;

   Energy() : intraBond(0.0), intraNonbond(0.0), inter(0.0), 
      tc(0.0), real(0.0), recip(0.0), self(0.0), correction(0.0), elect(0.0), total(0.0) {}
   Energy(double bond, double nonbond, double inter, double real, double recip, double self, double correction) :
      intraBond(bond), intraNonbond(nonbond), inter(inter), real(real), recip(recip),
	 self(self), correction(correction), tc(0.0), elect(real+recip+self+correction), total(0.0) {}

   //VALUE SETTERS
   double Total() 
   { total = intraBond + intraNonbond + inter + tc + real + recip + self + correction; return total; }
   void Zero() { 
      intraBond = 0.0;
      intraNonbond = 0.0;
      inter = 0.0;
      tc = 0.0;
	  real = 0.0;
	  recip = 0.0;
	  self = 0.0;
	  correction = 0.0;
	  elect = 0.0;
      total = 0.0; 
   }

   //OPERATORS
   Energy& operator-=(Intermolecular const& rhs)
   { inter -= rhs.energy; return *this; }
   Energy& operator+=(Intermolecular const& rhs)
   { inter += rhs.energy; return *this; }
   Energy& operator-=(Elect const& rhs)
   { elect -= rhs.energy; return *this; }
   Energy& operator+=(Elect const& rhs)
   { elect += rhs.energy; return *this; }
   Energy& operator-=(Energy const& rhs);
   Energy& operator+=(Energy const& rhs);
};

inline Energy& Energy::operator-=(Energy const& rhs)
{ 
   inter -= rhs.inter;
   intraBond -= rhs.intraBond;
   intraNonbond -= rhs.intraNonbond;
   tc -= rhs.tc;
   real -= rhs.real;
   recip -= rhs.recip;
   self -= rhs.self;
   correction -= rhs.correction;
   total -= rhs.total;
   elect -= rhs.elect;
   return *this; 
}

inline Energy& Energy::operator+=(Energy const& rhs)
{ 
   inter += rhs.inter;
   intraBond += rhs.intraBond;
   intraNonbond += rhs.intraNonbond;
   tc += rhs.tc;
   real += rhs.real;
   recip += rhs.recip;
   self += rhs.self;
   correction += rhs.correction;
   total += rhs.total;
   elect += rhs.elect;
   return *this; 
}

struct Virial
{
   //MEMBERS
   double inter, elect, tc, total;

   Virial() { Zero(); }

   //VALUE SETTERS
   double Total() { return total = inter + tc; }
   void Zero() { inter = tc = total = 0.0; }

   //OPERATORS
   Virial& operator-=(Virial const& rhs)
   { inter -= rhs.inter; tc -= rhs.tc; total -= rhs.total; return *this; }
   Virial& operator+=(Virial const& rhs)
   { inter += rhs.inter; tc += rhs.tc; total += rhs.total; return *this; }
   Virial& operator-=(Intermolecular const& rhs)
   { inter -= rhs.virial; return *this; }
   Virial& operator+=(Intermolecular const& rhs)
   { inter += rhs.virial; return *this; }
   Virial& operator-=(Elect const& rhs)
   { elect -= rhs.virial; return *this; }
   Virial& operator+=(Elect const& rhs)
   { elect += rhs.virial; return *this; }
   //For accounting for dimensionality
   Virial& operator=(Virial const& rhs)
   { inter = rhs.inter; tc = rhs.tc; total = rhs.total; return *this; } 
   Virial& operator/=(const double rhs)
   { inter /= rhs; tc /= rhs; total /= rhs; return *this; }

};


struct SystemPotential
{
   void Zero();
   double Total();
   void Add(const uint b, Intermolecular const& rhs)
   { boxVirial[b] += rhs; boxEnergy[b] += rhs; } 
   void Add(const uint b, Elect const& rhs)
   { boxVirial[b] += rhs; boxEnergy[b] += rhs; }
   void Add(const uint b, Energy const& en, Virial vir)
   { boxVirial[b] += vir; boxEnergy[b] += en; } 
   void Sub(const uint b, Energy const& en, Virial vir)
   { boxVirial[b] -= vir; boxEnergy[b] -= en; } 
   SystemPotential& operator=(SystemPotential const& rhs)
   {
      for (uint b = 0; b < BOX_TOTAL; b++)
      {
	 boxVirial[b] = rhs.boxVirial[b];
	 boxEnergy[b] = rhs.boxEnergy[b]; 
      }
      totalEnergy = rhs.totalEnergy;
      totalVirial = rhs.totalVirial;
      return *this;
   } 
   SystemPotential& operator+=(SystemPotential const& rhs)
   {
      for (uint b = 0; b < BOX_TOTAL; b++)
      {
	 boxVirial[b] += rhs.boxVirial[b];
	 boxEnergy[b] += rhs.boxEnergy[b]; 
      }
      Total();
      return *this;
   } 
   SystemPotential& operator-=(SystemPotential const& rhs)
   {
      for (uint b = 0; b < BOX_TOTAL; b++)
	 Add(b, rhs.boxEnergy[b], rhs.boxVirial[b]);
      Total();      return *this;
   } 

   Virial boxVirial[BOX_TOTAL], totalVirial; 
   Energy boxEnergy[BOX_TOTAL], totalEnergy; 
};

inline void SystemPotential::Zero()
{
   for (uint b = 0; b < BOX_TOTAL; b++)
   {
      boxEnergy[b].Zero();
      boxVirial[b].Zero();
   }
   totalEnergy.Zero();
   totalVirial.Zero();
}

inline double SystemPotential::Total()
{
   totalEnergy.Zero();
   totalVirial.Zero();
   for (uint b = 0; b < BOX_TOTAL; b++)
   {
      boxEnergy[b].Total();
      totalEnergy += boxEnergy[b];
      boxVirial[b].Total();
      totalVirial += boxVirial[b];
   }
   return totalEnergy.total;
}

#ifndef NDEBUG
inline std::ostream& operator << (std::ostream& out, const Energy& en)
{
   out << "Total: " << en.total << "   Inter: " << en.inter
       << "   IntraB: " << en.intraBond << "   IntraNB: "
       << en.intraNonbond << '\n';
   return out;
}
#endif

#endif
