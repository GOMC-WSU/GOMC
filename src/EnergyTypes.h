/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef ENERGYTYPES_H
#define ENERGYTYPES_H

/*
 *    EnergyTypes.h
 *    Defines structs containing energy values
 *
 */

#include "BasicTypes.h"           //For uint
#include "EnsemblePreprocessor.h" //For # box const.

// GJS
//#ifndef NDEBUG
#include <iostream>
//#endif

#ifndef BOXES_WITH_U_NB
#if ENSEMBLE == GCMC || ENSEMBLE == NVT || ENSEMBLE == NPT
#define BOXES_WITH_U_NB 1
#elif ENSEMBLE == GEMC
#define BOXES_WITH_U_NB 2
#endif
#endif

#ifndef BOXES_WITH_U_B
#if ENSEMBLE == NVT || ENSEMBLE == NPT
#define BOXES_WITH_U_B 1
#elif ENSEMBLE == GEMC || ENSEMBLE == GCMC
// case for NVT, GCMC
#define BOXES_WITH_U_B 2
#endif
#endif

struct Intermolecular {
  // MEMBERS
  double virial, energy;

  // CONSTRUCTORS
  Intermolecular() : virial(0.0), energy(0.0) {}
  Intermolecular(const double vir, const double en) : virial(vir), energy(en) {}

  // VALUE SETTER
  void Zero() { virial = energy = 0.0; }

  // OPERATORS
  Intermolecular &operator=(Intermolecular const &rhs) {
    virial = rhs.virial;
    energy = rhs.energy;
    return *this;
  }
  Intermolecular &operator-=(Intermolecular const &rhs) {
    virial -= rhs.virial;
    energy -= rhs.energy;
    return *this;
  }
  Intermolecular &operator+=(Intermolecular const &rhs) {
    virial += rhs.virial;
    energy += rhs.energy;
    return *this;
  }
  Intermolecular operator-(Intermolecular const &rhs) {
    return Intermolecular(virial - rhs.virial, energy - rhs.energy);
  }
  Intermolecular operator+(Intermolecular const &rhs) {
    return Intermolecular(virial + rhs.virial, energy + rhs.energy);
  }
};

class Energy {
public:
  Energy()
      : intraBond(0.0), intraNonbond(0.0), inter(0.0), tailCorrection(0.0),
        total(0.0), real(0.0), recip(0.0), self(0.0), correction(0.0),
        totalElect(0.0) {}
  Energy(double bond, double nonbond, double inter, double real, double recip,
         double self, double correc)
      : intraBond(bond), intraNonbond(nonbond), inter(inter),
        tailCorrection(0.0), total(0.0), real(real), recip(recip), self(self),
        correction(correc), totalElect(0.0) {}

  // VALUE SETTERS
  double Total() {
    total = intraBond + intraNonbond + inter + tailCorrection + real + recip +
            self + correction;
    return total;
  }

  double TotalElect() {
    totalElect = real + recip + self + correction;
    return totalElect;
  }

  void Zero() {
    intraBond = 0.0;
    intraNonbond = 0.0;
    inter = 0.0;
    tailCorrection = 0.0;
    real = 0.0;
    recip = 0.0;
    self = 0.0;
    correction = 0.0;
    totalElect = 0.0;
    total = 0.0;
  }

  // OPERATORS
  Energy &operator-=(Intermolecular const &rhs) {
    inter -= rhs.energy;
    return *this;
  }
  Energy &operator+=(Intermolecular const &rhs) {
    inter += rhs.energy;
    return *this;
  }
  Energy &operator-=(Energy const &rhs);
  Energy &operator+=(Energy const &rhs);
  Energy &operator*=(double const &rhs);

  // private:
  // MEMBERS
  double intraBond, intraNonbond, inter, tailCorrection, total, real, recip,
      self, correction, totalElect;
};

inline Energy &Energy::operator-=(Energy const &rhs) {
  inter -= rhs.inter;
  intraBond -= rhs.intraBond;
  intraNonbond -= rhs.intraNonbond;
  tailCorrection -= rhs.tailCorrection;
  real -= rhs.real;
  recip -= rhs.recip;
  self -= rhs.self;
  correction -= rhs.correction;
  totalElect -= rhs.totalElect;
  total -= rhs.total;

  return *this;
}

inline Energy &Energy::operator+=(Energy const &rhs) {
  inter += rhs.inter;
  intraBond += rhs.intraBond;
  intraNonbond += rhs.intraNonbond;
  tailCorrection += rhs.tailCorrection;
  real += rhs.real;
  recip += rhs.recip;
  self += rhs.self;
  correction += rhs.correction;
  totalElect += rhs.totalElect;
  total += rhs.total;

  return *this;
}

inline Energy &Energy::operator*=(double const &rhs) {
  inter *= rhs;
  intraBond *= rhs;
  intraNonbond *= rhs;
  tailCorrection *= rhs;
  real *= rhs;
  recip *= rhs;
  self *= rhs;
  correction *= rhs;
  totalElect *= rhs;
  total *= rhs;

  return *this;
}

class Virial {
public:
  Virial() { Zero(); }

  // VALUE SETTERS
  double Total() {
    TotalElect();
    total = inter + tailCorrection + real + recip + self + correction;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] =
            interTens[i][j] + realTens[i][j] + recipTens[i][j] + corrTens[i][j];
      }
    }
    return total;
  }

  double TotalElect() {
    totalElect = real + recip + self + correction;
    return totalElect;
  }

  void Zero() {
    inter = 0.0;
    tailCorrection = 0.0;
    real = 0.0;
    recip = 0.0;
    self = 0.0;
    correction = 0.0;
    totalElect = 0.0;
    total = 0.0;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] = 0.0;
        interTens[i][j] = 0.0;
        realTens[i][j] = 0.0;
        recipTens[i][j] = 0.0;
        corrTens[i][j] = 0.0;
      }
    }
  }

  // OPERATORS
  Virial &operator-=(Virial const &rhs) {
    inter -= rhs.inter;
    tailCorrection -= rhs.tailCorrection;
    real -= rhs.real;
    recip -= rhs.recip;
    self -= rhs.self;
    correction -= rhs.correction;
    totalElect -= rhs.totalElect;
    total -= rhs.total;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] -= rhs.totalTens[i][j];
        interTens[i][j] -= rhs.interTens[i][j];
        realTens[i][j] -= rhs.realTens[i][j];
        recipTens[i][j] -= rhs.recipTens[i][j];
        corrTens[i][j] -= rhs.corrTens[i][j];
      }
    }

    return *this;
  }

  Virial &operator+=(Virial const &rhs) {
    inter += rhs.inter;
    tailCorrection += rhs.tailCorrection;
    real += rhs.real;
    recip += rhs.recip;
    self += rhs.self;
    correction += rhs.correction;
    totalElect += rhs.totalElect;
    total += rhs.total;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] += rhs.totalTens[i][j];
        interTens[i][j] += rhs.interTens[i][j];
        realTens[i][j] += rhs.realTens[i][j];
        recipTens[i][j] += rhs.recipTens[i][j];
        corrTens[i][j] += rhs.corrTens[i][j];
      }
    }

    return *this;
  }

  Virial &operator-=(Intermolecular const &rhs) {
    inter -= rhs.virial;
    return *this;
  }
  Virial &operator+=(Intermolecular const &rhs) {
    inter += rhs.virial;
    return *this;
  }
  // For accounting for dimensionality
  Virial &operator=(Virial const &rhs) {
    inter = rhs.inter;
    tailCorrection = rhs.tailCorrection;
    real = rhs.real;
    recip = rhs.recip;
    self = rhs.self;
    correction = rhs.correction;
    totalElect = rhs.totalElect;
    total = rhs.total;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] = rhs.totalTens[i][j];
        interTens[i][j] = rhs.interTens[i][j];
        realTens[i][j] = rhs.realTens[i][j];
        recipTens[i][j] = rhs.recipTens[i][j];
        corrTens[i][j] = rhs.corrTens[i][j];
      }
    }

    return *this;
  }

  Virial &operator/=(const double rhs) {
    inter /= rhs;
    tailCorrection /= rhs;
    real /= rhs;
    recip /= rhs;
    self /= rhs;
    correction /= rhs;
    totalElect /= rhs;
    total /= rhs;

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        totalTens[i][j] /= rhs;
        interTens[i][j] /= rhs;
        realTens[i][j] /= rhs;
        recipTens[i][j] /= rhs;
        corrTens[i][j] /= rhs;
      }
    }

    return *this;
  }

  // private:
  // MEMBERS
  double inter, tailCorrection, real, recip, self, correction, totalElect,
      total;
  // Store the pressure tensor
  double interTens[3][3], realTens[3][3], recipTens[3][3], totalTens[3][3],
      corrTens[3][3];
};

class SystemPotential {
public:
  void Zero();
  double Total();
  bool ComparePotentials(SystemPotential &other);
  void Add(const uint b, Intermolecular const &rhs) {
    boxVirial[b] += rhs;
    boxEnergy[b] += rhs;
  }
  void Add(const uint b, Energy const &en, Virial vir) {
    boxVirial[b] += vir;
    boxEnergy[b] += en;
  }
  void Sub(const uint b, Energy const &en, Virial vir) {
    boxVirial[b] -= vir;
    boxEnergy[b] -= en;
  }
  SystemPotential &operator=(SystemPotential const &rhs) {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      boxVirial[b] = rhs.boxVirial[b];
      boxEnergy[b] = rhs.boxEnergy[b];
    }
    totalEnergy = rhs.totalEnergy;
    totalVirial = rhs.totalVirial;
    return *this;
  }
  SystemPotential &operator+=(SystemPotential const &rhs) {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      boxVirial[b] += rhs.boxVirial[b];
      boxEnergy[b] += rhs.boxEnergy[b];
    }
    Total();
    return *this;
  }
  SystemPotential &operator-=(SystemPotential const &rhs) {
    for (uint b = 0; b < BOX_TOTAL; b++)
      Add(b, rhs.boxEnergy[b], rhs.boxVirial[b]);
    Total();
    return *this;
  }

  Virial boxVirial[BOX_TOTAL], totalVirial;
  Energy boxEnergy[BOX_TOTAL], totalEnergy;
};

inline void SystemPotential::Zero() {
  for (uint b = 0; b < BOX_TOTAL; b++) {
    boxEnergy[b].Zero();
    boxVirial[b].Zero();
  }
  totalEnergy.Zero();
  totalVirial.Zero();
}

inline double SystemPotential::Total() {
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

inline bool SystemPotential::ComparePotentials(SystemPotential &other) {
  std::cout.precision(17);
  bool returnVal = true;
  if (totalEnergy.inter != other.totalEnergy.inter) {
    std::cout << "my inter: " << totalEnergy.inter
              << "  other inter: " << other.totalEnergy.inter << std::endl;
    std::cout << "difference: " << totalEnergy.inter - other.totalEnergy.inter
              << std::endl;
    returnVal = false;
  }

  if (totalEnergy.intraBond != other.totalEnergy.intraBond) {
    std::cout << "my intraBond: " << totalEnergy.intraBond
              << "  other intraBond: " << other.totalEnergy.intraBond
              << std::endl;
    std::cout << "difference: "
              << totalEnergy.intraBond - other.totalEnergy.intraBond
              << std::endl;
    returnVal = false;
  }
  if (totalEnergy.intraNonbond != other.totalEnergy.intraNonbond) {
    std::cout << "my intraNonbond: " << totalEnergy.intraNonbond
              << "  other intraNonbond: " << other.totalEnergy.intraNonbond
              << std::endl;
    std::cout << "difference: "
              << totalEnergy.intraNonbond - other.totalEnergy.intraNonbond
              << std::endl;
    returnVal = false;
  }
  if (totalEnergy.tailCorrection != other.totalEnergy.tailCorrection) {
    std::cout << "my LRC: " << totalEnergy.tailCorrection
              << "  other LRC: " << other.totalEnergy.tailCorrection
              << std::endl;
    std::cout << "difference: "
              << totalEnergy.tailCorrection - other.totalEnergy.tailCorrection
              << std::endl;
    returnVal = false;
  }
  if (totalEnergy.real != other.totalEnergy.real) {
    std::cout << "my real: " << totalEnergy.real
              << "  other real: " << other.totalEnergy.real << std::endl;
    std::cout << "difference: " << totalEnergy.real - other.totalEnergy.real
              << std::endl;
    returnVal = false;
  }
  if (totalEnergy.recip != other.totalEnergy.recip) {
    std::cout << "my recip: " << totalEnergy.recip
              << "  other recip: " << other.totalEnergy.recip << std::endl;
    std::cout << "difference: " << totalEnergy.recip - other.totalEnergy.recip
              << std::endl;
    returnVal = false;
  }
  if (totalEnergy.intraNonbond != other.totalEnergy.intraNonbond) {
    std::cout << "my intraNonbond: " << totalEnergy.intraNonbond
              << "  other intraNonbond: " << other.totalEnergy.intraNonbond
              << std::endl;
    std::cout << "difference: "
              << totalEnergy.intraNonbond - other.totalEnergy.intraNonbond
              << std::endl;
    returnVal = false;
  }
  if (totalEnergy.tailCorrection != other.totalEnergy.tailCorrection) {
    std::cout << "my LRC: " << totalEnergy.tailCorrection
              << "  other LRC: " << other.totalEnergy.tailCorrection
              << std::endl;
    std::cout << "difference: "
              << totalEnergy.tailCorrection - other.totalEnergy.tailCorrection
              << std::endl;
    returnVal = false;
  }
  if (totalEnergy.totalElect != other.totalEnergy.totalElect) {
    std::cout << "my totalElect: " << totalEnergy.totalElect
              << "  other totalElect: " << other.totalEnergy.totalElect
              << std::endl;
    std::cout << "difference: "
              << totalEnergy.totalElect - other.totalEnergy.totalElect
              << std::endl;
    if (totalEnergy.real != other.totalEnergy.real) {
      std::cout << "my real: " << totalEnergy.real
                << "  other real: " << other.totalEnergy.real << std::endl;
      std::cout << "difference: " << totalEnergy.real - other.totalEnergy.real
                << std::endl;
    }
    if (totalEnergy.recip != other.totalEnergy.recip) {
      std::cout << "my recip: " << totalEnergy.recip
                << "  other recip: " << other.totalEnergy.recip << std::endl;
      std::cout << "difference: " << totalEnergy.recip - other.totalEnergy.recip
                << std::endl;
    }
    if (totalEnergy.self != other.totalEnergy.self) {
      std::cout << "my self: " << totalEnergy.self
                << "  other self: " << other.totalEnergy.self << std::endl;
      std::cout << "difference: " << totalEnergy.self - other.totalEnergy.self
                << std::endl;
    }
    if (totalEnergy.correction != other.totalEnergy.correction) {
      std::cout << "my correction: " << totalEnergy.correction
                << "  other correction: " << other.totalEnergy.correction
                << std::endl;
      std::cout << "difference: "
                << totalEnergy.correction - other.totalEnergy.correction
                << std::endl;
    }
    returnVal = false;
  }
  if (totalEnergy.total != other.totalEnergy.total) {
    std::cout << "my total: " << totalEnergy.total
              << "  other total: " << other.totalEnergy.total << std::endl;
    std::cout << "difference: " << totalEnergy.total - other.totalEnergy.total
              << std::endl;
    returnVal = false;
  }
  return returnVal;
}

#ifndef NDEBUG
inline std::ostream &operator<<(std::ostream &out, Energy &en) {
  // Save existing settings for the ostream
  std::streamsize ss = out.precision();
  std::ios_base::fmtflags ff = out.flags();

  en.Total();
  en.TotalElect();

  out << std::setprecision(6) << std::fixed;
  out << "\tTotal: " << en.total << "  IntraB: " << en.intraBond
      << "  IntraNB: " << en.intraNonbond << "  Inter: " << en.inter
      << "  LRC: " << en.tailCorrection;
  if (en.totalElect != 0.0) {
    out << std::endl
        << "\tTotal Electric: " << en.totalElect << "  Real: " << en.real
        << "  Recip: " << en.recip << "  Self: " << en.self
        << "  Correction: " << en.correction;
  }

  // Restore ostream settings to prior value
  out << std::setprecision(ss);
  out.flags(ff);
  return out;
}
#endif

#endif
