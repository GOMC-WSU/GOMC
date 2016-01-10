/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
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
   Intermolecular operator-(Intermolecular const& rhs) 
   { return Intermolecular(virial - rhs.virial, energy - rhs.energy); }
   Intermolecular operator+(Intermolecular const& rhs) 
   { return Intermolecular(virial + rhs.virial, energy + rhs.energy); }
};

class Energy
{
public:
   Energy() : intraBond(0.0), intraNonbond(0.0), inter(0.0), 
      tc(0.0), total(0.0) {}
   Energy(double bond, double nonbond, double inter) :
      intraBond(bond), intraNonbond(nonbond), inter(inter),
	 tc(0.0), total(0.0) {}

   //VALUE SETTERS
   double Total() 
   { total = intraBond + intraNonbond + inter + tc; return total; }
   void Zero() { 
      intraBond = 0.0;
      intraNonbond = 0.0;
      inter = 0.0;
      tc = 0.0;
      total = 0.0; 
   }

   //OPERATORS
   Energy& operator-=(Intermolecular const& rhs)
   { inter -= rhs.energy; return *this; }
   Energy& operator+=(Intermolecular const& rhs)
   { inter += rhs.energy; return *this; }
   Energy& operator-=(Energy const& rhs);
   Energy& operator+=(Energy const& rhs);

//private:
   //MEMBERS
   double intraBond, intraNonbond, inter, tc, total;
};

inline Energy& Energy::operator-=(Energy const& rhs)
{ 
   inter -= rhs.inter;
   intraBond -= rhs.intraBond;
   intraNonbond -= rhs.intraNonbond;
   tc -= rhs.tc;
   total -= rhs.total;
   return *this; 
}

inline Energy& Energy::operator+=(Energy const& rhs)
{ 
   inter += rhs.inter;
   intraBond += rhs.intraBond;
   intraNonbond += rhs.intraNonbond;
   tc += rhs.tc;
   total += rhs.total;
   return *this; 
}

class Virial
{
public:
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
   //For accounting for dimensionality
   Virial& operator=(Virial const& rhs)
   { inter = rhs.inter; tc = rhs.tc; total = rhs.total; return *this; } 
   Virial& operator/=(const double rhs)
   { inter /= rhs; tc /= rhs; total /= rhs; return *this; }

//private:
   //MEMBERS
   double inter, tc, total;
};


class SystemPotential
{
public:
   void Zero();
   double Total();
   void Add(const uint b, Intermolecular const& rhs)
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

