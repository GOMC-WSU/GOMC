/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
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
#include "BasicTypes.h" //For uint

#ifndef NDEBUG
#include <iostream>
#endif

#ifndef BOXES_WITH_U_NB
#if ENSEMBLE == GCMC || ENSEMBLE == NVT || ENSEMBLE == NPT
#define BOXES_WITH_U_NB 1
#elif ENSEMBLE == GEMC
//case for NVT, GCMC
#define BOXES_WITH_U_NB 2
#endif
#endif

#ifndef BOXES_WITH_U_B
#if ENSEMBLE == NVT || ENSEMBLE == NPT
#define BOXES_WITH_U_B 1
#elif ENSEMBLE == GEMC || ENSEMBLE == GCMC
//case for NVT, GCMC
#define BOXES_WITH_U_B 2
#endif
#endif

struct Intermolecular {
  //MEMBERS
  real virial, energy;

  //CONSTRUCTORS
  Intermolecular() : virial(0.0), energy(0.0) {}
  Intermolecular(const real vir, const real en) :
    virial(vir), energy(en) {}

  //VALUE SETTER
  void Zero()
  {
    virial = energy = 0.0;
  }

  //OPERATORS
  Intermolecular& operator=(Intermolecular const& rhs)
  {
    virial = rhs.virial;
    energy = rhs.energy;
    return *this;
  }
  Intermolecular& operator-=(Intermolecular const& rhs)
  {
    virial -= rhs.virial;
    energy -= rhs.energy;
    return *this;
  }
  Intermolecular& operator+=(Intermolecular const& rhs)
  {
    virial += rhs.virial;
    energy += rhs.energy;
    return *this;
  }
  Intermolecular operator-(Intermolecular const& rhs)
  {
    return Intermolecular(virial - rhs.virial, energy - rhs.energy);
  }
  Intermolecular operator+(Intermolecular const& rhs)
  {
    return Intermolecular(virial + rhs.virial, energy + rhs.energy);
  }
};

class Energy
{
public:
  Energy() : intraBond(0.0), intraNonbond(0.0), inter(0.0),
    tc(0.0), total(0.0), real_en(0.0), recip(0.0), self(0.0),
    correction(0.0), totalElect(0.0) {}
  Energy(real bond, real nonbond, real inter, real real_en,
         real recip, real self, real correc) :
    intraBond(bond), intraNonbond(nonbond), inter(inter),
    tc(0.0), real_en(real_en), recip(recip), self(self), correction(correc),
    totalElect(0.0), total(0.0) {}

  //VALUE SETTERS
  real Total()
  {
    total = intraBond + intraNonbond + inter + tc + real_en + recip + self +
            correction;
    return total;
  }

  real TotalElect()
  {
    totalElect = real_en + recip + self + correction;
    return totalElect;
  }

  void Zero()
  {
    intraBond = 0.0;
    intraNonbond = 0.0;
    inter = 0.0;
    tc = 0.0;
    real_en = 0.0;
    recip = 0.0;
    self = 0.0;
    correction = 0.0;
    totalElect = 0.0;
    total = 0.0;
  }

  //OPERATORS
  Energy& operator-=(Intermolecular const& rhs)
  {
    inter -= rhs.energy;
    return *this;
  }
  Energy& operator+=(Intermolecular const& rhs)
  {
    inter += rhs.energy;
    return *this;
  }
  Energy& operator-=(Energy const& rhs);
  Energy& operator+=(Energy const& rhs);

//private:
  //MEMBERS
  real intraBond, intraNonbond, inter, tc, total, real_en, recip, self,
         correction, totalElect;
};

inline Energy& Energy::operator-=(Energy const& rhs)
{
  inter -= rhs.inter;
  intraBond -= rhs.intraBond;
  intraNonbond -= rhs.intraNonbond;
  tc -= rhs.tc;
  real_en -= rhs.real_en;
  recip -= rhs.recip;
  self -= rhs.self;
  correction -= rhs.correction;
  totalElect -= rhs.totalElect;
  total -= rhs.total;

  return *this;
}

inline Energy& Energy::operator+=(Energy const& rhs)
{
  inter += rhs.inter;
  intraBond += rhs.intraBond;
  intraNonbond += rhs.intraNonbond;
  tc += rhs.tc;
  real_en += rhs.real_en;
  recip += rhs.recip;
  self += rhs.self;
  correction += rhs.correction;
  totalElect += rhs.totalElect;
  total += rhs.total;

  return *this;
}

class Virial
{
public:
  Virial()
  {
    Zero();
  }

  //VALUE SETTERS
  real Total()
  {
    TotalElect();
    total = inter + tc + real_en + recip + self + correction;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] = interTens[i][j] + realTens[i][j] + recipTens[i][j] +
                          corrTens[i][j];
      }
    }
    return total;
  }

  real TotalElect()
  {
    totalElect = real_en + recip + self + correction;
    return totalElect;
  }

  void Zero()
  {
    inter = 0.0;
    tc = 0.0;
    real_en = 0.0;
    recip = 0.0;
    self = 0.0;
    correction = 0.0;
    totalElect = 0.0;
    total = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] = 0.0;
        interTens[i][j] = 0.0;
        realTens[i][j]  = 0.0;
        recipTens[i][j] = 0.0;
        corrTens[i][j]  = 0.0;
      }
    }
  }

  //OPERATORS
  Virial& operator-=(Virial const& rhs)
  {
    inter -= rhs.inter;
    tc -= rhs.tc;
    real_en -= rhs.real_en;
    recip -= rhs.recip;
    self -= rhs.self;
    correction -= rhs.correction;
    totalElect -= rhs.totalElect;
    total -= rhs.total;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] -= rhs.totalTens[i][j];
        interTens[i][j] -= rhs.interTens[i][j];
        realTens[i][j]  -= rhs.realTens[i][j];
        recipTens[i][j] -= rhs.recipTens[i][j];
        corrTens[i][j]  -= rhs.corrTens[i][j];
      }
    }

    return *this;
  }


  Virial& operator+=(Virial const& rhs)
  {
    inter += rhs.inter;
    tc += rhs.tc;
    real_en += rhs.real_en;
    recip += rhs.recip;
    self += rhs.self;
    correction += rhs.correction;
    totalElect += rhs.totalElect;
    total += rhs.total;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] += rhs.totalTens[i][j];
        interTens[i][j] += rhs.interTens[i][j];
        realTens[i][j]  += rhs.realTens[i][j];
        recipTens[i][j] += rhs.recipTens[i][j];
        corrTens[i][j]  += rhs.corrTens[i][j];
      }
    }

    return *this;
  }

  Virial& operator-=(Intermolecular const& rhs)
  {
    inter -= rhs.virial;
    return *this;
  }
  Virial& operator+=(Intermolecular const& rhs)
  {
    inter += rhs.virial;
    return *this;
  }
  //For accounting for dimensionality
  Virial& operator=(Virial const& rhs)
  {
    inter = rhs.inter;
    tc = rhs.tc;
    real_en = rhs.real_en;
    recip = rhs.recip;
    self = rhs.self;
    correction = rhs.correction;
    totalElect = rhs.totalElect;
    total = rhs.total;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] = rhs.totalTens[i][j];
        interTens[i][j] = rhs.interTens[i][j];
        realTens[i][j]  = rhs.realTens[i][j];
        recipTens[i][j] = rhs.recipTens[i][j];
        corrTens[i][j]  = rhs.corrTens[i][j];
      }
    }

    return *this;
  }

  Virial& operator/=(const real rhs)
  {
    inter /= rhs;
    tc /= rhs;
    real_en /= rhs;
    recip /= rhs;
    self /= rhs;
    correction /= rhs;
    totalElect /= rhs;
    total /= rhs;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] /= rhs;
        interTens[i][j] /= rhs;
        realTens[i][j]  /= rhs;
        recipTens[i][j] /= rhs;
        corrTens[i][j]  /= rhs;
      }
    }

    return *this;
  }

//private:
  //MEMBERS
  real inter, tc, real_en, recip, self, correction, totalElect, total;
  //Store the pressure tensor
  real interTens[3][3], realTens[3][3], recipTens[3][3], totalTens[3][3],
         corrTens[3][3];
};


class SystemPotential
{
public:
  void Zero();
  real Total();
  void Add(const uint b, Intermolecular const& rhs)
  {
    boxVirial[b] += rhs;
    boxEnergy[b] += rhs;
  }
  void Add(const uint b, Energy const& en, Virial vir)
  {
    boxVirial[b] += vir;
    boxEnergy[b] += en;
  }
  void Sub(const uint b, Energy const& en, Virial vir)
  {
    boxVirial[b] -= vir;
    boxEnergy[b] -= en;
  }
  SystemPotential& operator=(SystemPotential const& rhs)
  {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      boxVirial[b] = rhs.boxVirial[b];
      boxEnergy[b] = rhs.boxEnergy[b];
    }
    totalEnergy = rhs.totalEnergy;
    totalVirial = rhs.totalVirial;
    return *this;
  }
  SystemPotential& operator+=(SystemPotential const& rhs)
  {
    for (uint b = 0; b < BOX_TOTAL; b++) {
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
    Total();
    return *this;
  }

  Virial boxVirial[BOX_TOTAL], totalVirial;
  Energy boxEnergy[BOX_TOTAL], totalEnergy;
};

inline void SystemPotential::Zero()
{
  for (uint b = 0; b < BOX_TOTAL; b++) {
    boxEnergy[b].Zero();
    boxVirial[b].Zero();
  }
  totalEnergy.Zero();
  totalVirial.Zero();
}

inline real SystemPotential::Total()
{
  totalEnergy.Zero();
  totalVirial.Zero();
  for (uint b = 0; b < BOX_TOTAL; b++) {
    boxEnergy[b].Total();
    boxEnergy[b].TotalElect();
    totalEnergy += boxEnergy[b];

    boxVirial[b].Total();
    boxVirial[b].TotalElect();
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
