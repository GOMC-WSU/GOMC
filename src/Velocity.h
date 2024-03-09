/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef VELOCITY_H
#define VELOCITY_H

#include "ConfigSetup.h"
#include "Forcefield.h"
#include "MoleculeLookup.h"
#include "Molecules.h"
#include "PDBSetup.h"
#include "PRNG.h"
#include "UnitConst.h"
#include "XYZArray.h"

namespace config_setup {
struct Input;
}

class Velocity : public XYZArray {
public:
  Velocity(Forcefield &ff, MoleculeLookup &molLook, Molecules const &mol,
           PRNG &prng)
      : prngRef(prng), molRef(mol), molLookRef(molLook), hasVelocity(false),
        temperature(ff.T_in_K) {}

  void Init(pdb_setup::Atoms const &atoms,
            config_setup::Input const &inputFile) {
    if (inputFile.restart.restartFromBinaryVelFile) {
      hasVelocity = true;
      // Allocate master array
      XYZArray::Init(atoms.beta.size());
      for (int b = 0; b < BOX_TOTAL; ++b) {
        updateBoxVelocity[b] = false;
      }
    } else {
      hasVelocity = false;
    }
    // We dont set the velocity values. We read it from restart
    // binary velocity outputted by MD simulation
  }

  void UpdateBoxVelocity(const uint &box) { updateBoxVelocity[box] = true; }

  void UpdateMolVelocity(const uint &mol, const uint &box) {
    // Do not do anything if we did not read velocity
    // If we performed a move that requires to update all atom's
    // velocity, we do not waste time and return
    if (!hasVelocity || updateBoxVelocity[box]) {
      return;
    }

    molRef.GetRangeStartStop(pStart, pEnd, mol);
    for (uint p = pStart; p < pEnd; ++p) {
      mass = molRef.GetKind(mol).atomMass[p - pStart];
      this->Set(p, Random_velocity(mass));
    }
  }

  void UpdateVelocityInBox(const int box) {
    // Do not do anything if we did not read velocity
    // If we Did not performed a move that requires to update all atom's
    // velocity, we should not update the velocity
    if (!updateBoxVelocity[box] || !hasVelocity) {
      return;
    }

    MoleculeLookup::box_iterator m = molLookRef.BoxBegin(box),
                                 end = molLookRef.BoxEnd(box);
    // loop thorugh all molecule in box
    while (m != end) {
      molRef.GetRangeStartStop(pStart, pEnd, *m);
      for (uint p = pStart; p < pEnd; ++p) {
        mass = molRef.GetKind(*m).atomMass[p - pStart];
        this->Set(p, Random_velocity(mass));
      }
      ++m;
    }
  }

  Velocity &operator=(Velocity const &rhs) {
    this->XYZArray::operator=(rhs);
    for (int b = 0; b < BOX_TOTAL; ++b) {
      updateBoxVelocity[b] = rhs.updateBoxVelocity[b];
    }
    temperature = rhs.temperature;
    return *this;
  }

private:
  //  The following comment was taken from X-PLOR where
  //  the following section of code was adapted from.
  //  NAMD uses same function
  //  This section generates a Gaussian random
  //  deviate of 0.0 mean and standard deviation RFD for
  //  each of the three spatial dimensions.
  //  The algorithm is a "sum of uniform deviates algorithm"
  //  which may be found in Abramowitz and Stegun,
  //  "Handbook of Mathematical Functions", pg 952.
  XYZ Random_velocity(const double &mass) {
    XYZ vel;
    int i;
    if (mass <= 0.0) {
      kbToverM = 0.0;
    } else {
      kbToverM = sqrt(temperature * unit::BOLTZMANN / mass);
    }

    for (randnum = 0.0, i = 0; i < 12; ++i) {
      randnum += prngRef();
    }

    randnum -= 6.0;
    vel.x = randnum * kbToverM;

    for (randnum = 0.0, i = 0; i < 12; ++i) {
      randnum += prngRef();
    }

    randnum -= 6.0;
    vel.y = randnum * kbToverM;

    for (randnum = 0.0, i = 0; i < 12; ++i) {
      randnum += prngRef();
    }

    randnum -= 6.0;
    vel.z = randnum * kbToverM;

    return vel;
  }

  PRNG &prngRef;
  Molecules const &molRef;
  MoleculeLookup &molLookRef;
  uint pStart, pEnd; // start and end of atom index
  double randnum;    // Random number from -6.0 to 6.0
  double kbToverM;   // sqrt(Kb*Temp/Mass)
  double mass;       // Mass of particle

  bool hasVelocity;
  bool updateBoxVelocity[BOX_TOTAL]; // update all atom's velocities
  double &temperature;               // system temperature
};

#endif