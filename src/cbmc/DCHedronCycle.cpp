/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#define _USE_MATH_DEFINES
#include "DCHedronCycle.h"

#include <cassert>
#include <cmath>
#include <numeric>

#include "DCData.h"
#include "Forcefield.h"
#include "Geometry.h"
#include "MolSetup.h"
#include "NumLib.h"
#include "PRNG.h"
#include "TrialMol.h"

namespace {
// Wish I could use lambdas..
struct FindA1 {
  FindA1(uint x) : x(x){};
  bool operator()(const mol_setup::Bond &b) { return (b.a1 == x); }
  uint x;
};

struct FindAngle {
  FindAngle(uint x, uint y) : x(x), y(y) {}
  uint x, y;
  bool operator()(const mol_setup::Angle &a) {
    return (a.a0 == x && a.a2 == y) || (a.a0 == y && a.a2 == x);
  }
};

// Check to see if angle a is in the ring or not
bool IsInRing(const std::vector<int> &cycAtoms, const mol_setup::Angle &a) {
  bool res = true;
  if (std::find(cycAtoms.begin(), cycAtoms.end(), a.a0) == cycAtoms.end()) {
    res = false;
  } else if (std::find(cycAtoms.begin(), cycAtoms.end(), a.a1) ==
             cycAtoms.end()) {
    res = false;
  } else if (std::find(cycAtoms.begin(), cycAtoms.end(), a.a2) ==
             cycAtoms.end()) {
    res = false;
  }
  return res;
}

} // namespace

namespace cbmc {

DCHedronCycle::DCHedronCycle(DCData *data, const mol_setup::MolKind &kind,
                             const std::vector<int> &cycAtoms, uint focus,
                             uint prev)
    : data(data), focus(focus), prev(prev) {
  using namespace mol_setup;
  std::vector<Bond> onFocus = AtomBonds(kind, focus);
  onFocus.erase(remove_if(onFocus.begin(), onFocus.end(), FindA1(prev)),
                onFocus.end());
  nBonds = onFocus.size();

  for (uint i = 0; i < nBonds; ++i) {
    bonded[i] = onFocus[i].a1;
  }

  std::vector<Angle> angles = AtomMidAngles(kind, focus);
  double sumAngle = 0.0;
  for (uint a = 0; a < angles.size(); a++) {
    sumAngle += data->ff.angles->Angle(angles[a].kind);
  }
  // If sum of angles (sumAngle +/- 10) ~ 2*pi = 6.283, it means they are in a
  // plane To avoid geometric conflict for flexible angle, we consider it fixed
  // and let crankshaft sample it. 0.1745 ~= 10 degree tolerance
  bool angleInPlane = (std::abs(2.0 * M_PI - sumAngle) < 0.1745);
  bool constrainAngInRing = false;

  for (uint i = 0; i < nBonds; ++i) {
    typedef std::vector<Angle>::const_iterator Aiter;
    Aiter free =
        find_if(angles.begin(), angles.end(), FindAngle(prev, bonded[i]));
    assert(free != angles.end());
    angleKinds[i][i] = free->kind;
    angleInRing[i][i] = IsInRing(cycAtoms, *free);

    for (uint j = i + 1; j < nBonds; ++j) {
      Aiter pair = find_if(angles.begin(), angles.end(),
                           FindAngle(bonded[i], bonded[j]));
      angleKinds[i][j] = pair->kind;
      angleKinds[j][i] = pair->kind;
      angleInRing[i][j] = IsInRing(cycAtoms, *pair);
      angleInRing[j][i] = angleInRing[i][j];
      constrainAngInRing |= angleInRing[i][j];
    }
  }

  // If one of the constrained angles belongs to a ring and the angle forms a
  // plane
  // we fix the free and constrained angles
  if (constrainAngInRing && angleInPlane) {
    for (uint i = 0; i < nBonds; ++i) {
      angleInRing[i][i] = true;
      for (uint j = i + 1; j < nBonds; ++j) {
        angleInRing[i][j] = true;
        angleInRing[j][i] = true;
      }
    }
  }

  phi[0] = 0.0;
  phiWeight[0] = 1.0;

  if (data->nAngleTrials < 1) {
    std::cout << "Error: CBMC angle trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

void DCHedronCycle::SetBondNew(double const *bondLen, double const &anchBond) {
  for (uint i = 0; i < nBonds; ++i) {
    bondLength[i] = bondLen[i];
  }
  anchorBond = anchBond;
}

void DCHedronCycle::SetBondOld(double const *bondLen, double const &anchBond) {
  for (uint i = 0; i < nBonds; ++i) {
    bondLengthOld[i] = bondLen[i];
  }
  anchorBondOld = anchBond;
}

double DCHedronCycle::GetWeight() {
  double result = 1;
  for (uint i = 0; i < nBonds; ++i) {
    result *= thetaWeight[i];
    result *= phiWeight[i];
  }
  return result;
}

double DCHedronCycle::CalcTheta(TrialMol &mol, const uint a0, const uint a1,
                                const uint a2) {
  double theta = 0.0;
  if (mol.AtomExists(a0) && mol.AtomExists(a2)) {
    // Calculate theta using tCoords
    theta = mol.GetTheta(a0, a1, a2);
  } else {
    // Calculate theta using bCoords
    const XYZArray &coords = mol.GetBCoords();
    // Since data in bCoords in unwrap, no need to calc minImage
    theta = geom::Theta(coords.Difference(a0, a1), coords.Difference(a2, a1));
  }
  return theta;
}

// Randomly generate nTrials angles and save energy and weight
void DCHedronCycle::GenerateAnglesNew(TrialMol &newMol, uint molIndex,
                                      uint kind, uint nTrials, uint bType) {
  double *nonbonded_1_3 = data->nonbonded_1_3;
  bool angleFix = data->ff.angles->AngleFixed(kind);
  std::fill_n(nonbonded_1_3, nTrials, 0.0);
  // use backup coordinate to find theta and phi of the ring
  if (angleInRing[bType][bType] || angleFix) {
    double th = CalcTheta(newMol, bonded[bType], focus, prev);
    std::fill_n(data->angles, nTrials, th);
    double en = data->ff.angles->Calc(kind, th);
    std::fill_n(data->angleEnergy, nTrials, en);
    double distSq = newMol.AngleDist(anchorBond, bondLength[bType], th);
    double enNB =
        data->calc.IntraEnergy_1_3(distSq, prev, bonded[bType], molIndex);
    std::fill_n(nonbonded_1_3, nTrials, enNB);
    double w = exp((en + enNB) * -data->ff.beta);
    std::fill_n(data->angleWeights, nTrials, w);
    return;
  }

  for (int i = 0; i < (int)nTrials; ++i) {
    data->angles[i] = data->prng.rand(M_PI);
  }

#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(bType, kind, molIndex, newMol, nonbonded_1_3, nTrials)
#endif
  for (int i = 0; i < (int)nTrials; ++i) {
    data->angleEnergy[i] = data->ff.angles->Calc(kind, data->angles[i]);
    double distSq =
        newMol.AngleDist(anchorBond, bondLength[bType], data->angles[i]);
    nonbonded_1_3[i] =
        data->calc.IntraEnergy_1_3(distSq, prev, bonded[bType], molIndex);
    data->angleWeights[i] =
        exp((data->angleEnergy[i] + nonbonded_1_3[i]) * -data->ff.beta);
  }
}

void DCHedronCycle::GenerateAnglesOld(TrialMol &oldMol, uint molIndex,
                                      uint kind, uint nTrials, uint bType) {
  double *nonbonded_1_3 = data->nonbonded_1_3;
  bool angleFix = data->ff.angles->AngleFixed(kind);
  std::fill_n(nonbonded_1_3, nTrials, 0.0);
  // use backup coordinate to find theta and phi of the ring
  if (angleInRing[bType][bType] || angleFix) {
    double th = CalcTheta(oldMol, bonded[bType], focus, prev);
    std::fill_n(data->angles, nTrials, th);
    double en = data->ff.angles->Calc(kind, th);
    std::fill_n(data->angleEnergy, nTrials, en);
    double distSq = oldMol.AngleDist(anchorBondOld, bondLengthOld[bType], th);
    double enNB =
        data->calc.IntraEnergy_1_3(distSq, prev, bonded[bType], molIndex);
    std::fill_n(nonbonded_1_3, nTrials, enNB);
    double w = exp((en + enNB) * -data->ff.beta);
    std::fill_n(data->angleWeights, nTrials, w);
    return;
  }

  for (int i = 0; i < (int)nTrials; ++i) {
    data->angles[i] = data->prng.rand(M_PI);
  }

#ifdef _OPENMP
#pragma omp parallel for default(none)                                         \
    shared(bType, kind, molIndex, nonbonded_1_3, nTrials, oldMol)
#endif
  for (int i = 0; i < (int)nTrials; ++i) {
    data->angleEnergy[i] = data->ff.angles->Calc(kind, data->angles[i]);

    double distSq =
        oldMol.AngleDist(anchorBondOld, bondLengthOld[bType], data->angles[i]);
    nonbonded_1_3[i] =
        data->calc.IntraEnergy_1_3(distSq, prev, bonded[bType], molIndex);

    data->angleWeights[i] =
        exp((data->angleEnergy[i] + nonbonded_1_3[i]) * -data->ff.beta);
  }
}

void DCHedronCycle::FreeAnglesNew(TrialMol &newMol, uint molIndex,
                                  uint nTrials) {
  for (uint i = 0; i < nBonds; ++i) {
    GenerateAnglesNew(newMol, molIndex, angleKinds[i][i], nTrials, i);
    double stepWeight =
        std::accumulate(data->angleWeights, data->angleWeights + nTrials, 0.0);
    uint winner =
        data->prng.PickWeighted(data->angleWeights, nTrials, stepWeight);
    theta[i] = data->angles[winner];
    bendEnergy += data->angleEnergy[winner];
    oneThree += data->nonbonded_1_3[winner];
    thetaWeight[i] = stepWeight;
  }
}

void DCHedronCycle::FreeAnglesOld(TrialMol &oldMol, uint molIndex,
                                  uint nTrials) {
  for (uint i = 0; i < nBonds; ++i) {
    GenerateAnglesOld(oldMol, molIndex, angleKinds[i][i], nTrials, i);
    double stepWeight =
        std::accumulate(data->angleWeights, data->angleWeights + nTrials, 0.0);
    thetaWeight[i] = stepWeight;
  }
}

void DCHedronCycle::PrepareNew(TrialMol &newMol, uint molIndex) {
  bendEnergy = 0.0;
  oneThree = 0.0;
  phi[0] = 0.0;
  FreeAnglesNew(newMol, molIndex, data->nAngleTrials);
  ConstrainedAngles(newMol, molIndex, data->nAngleTrials);
}

void DCHedronCycle::PrepareOld(TrialMol &oldMol, uint molIndex) {
  oneThree = 0.0;
  bendEnergy = 0.0;
  FreeAnglesOld(oldMol, molIndex, data->nAngleTrials - 1);
}

void DCHedronCycle::IncorporateOld(TrialMol &oldMol, uint molIndex) {
  bendEnergy = 0.0;
  oneThree = 0.0;
  const Forcefield &ff = data->ff;
  for (uint b = 0; b < nBonds; ++b) {
    oldMol.OldThetaAndPhi(bonded[b], focus, theta[b], phi[b]);
    double thetaEnergy = data->ff.angles->Calc(angleKinds[b][b], theta[b]);
    double distSq = oldMol.OldDistSq(prev, bonded[b]);
    double nonbondedEn =
        data->calc.IntraEnergy_1_3(distSq, prev, bonded[b], molIndex);

    thetaWeight[b] += exp(-1 * data->ff.beta * (thetaEnergy + nonbondedEn));
    bendEnergy += thetaEnergy;
    oneThree += nonbondedEn;

    if (b != 0) {
      double phiEnergy = 0.0;
      nonbondedEn = 0.0;
      phiWeight[b] = 0.0;
      for (uint c = 0; c < b; ++c) {
        double cosTerm = cos(theta[b]) * cos(theta[c]);
        double sinTerm = sin(theta[b]) * sin(theta[c]);
        double bfcTheta = acos(sinTerm * cos(phi[b] - phi[c]) + cosTerm);

        double distSq = oldMol.OldDistSq(bonded[c], bonded[b]);
        nonbondedEn +=
            data->calc.IntraEnergy_1_3(distSq, bonded[c], bonded[b], molIndex);

        phiEnergy += ff.angles->Calc(angleKinds[b][c], bfcTheta);
      }
      phiWeight[b] = exp(-ff.beta * (phiEnergy + nonbondedEn));
      bendEnergy += phiEnergy;
      oneThree += nonbondedEn;
    }
  }
}

void DCHedronCycle::ConstrainedAngles(TrialMol &newMol, uint molIndex,
                                      uint nTrials) {
  double *angles = data->angles;
  double *energies = data->angleEnergy;
  double *weights = data->angleWeights;
  double *nonbonded_1_3 = data->nonbonded_1_3;

  SetBasis(newMol, focus, prev);
  phi[0] = 0.0;
  // reference phi
  double refPhi = CalcOldPhi(newMol, bonded[0], focus);

  for (uint b = 1; b < nBonds; ++b) {
    std::fill_n(energies, nTrials, 0.0);
    std::fill_n(nonbonded_1_3, nTrials, 0.0);
    // pick "twist" angles
    for (uint i = 0; i < nTrials; ++i) {
      angles[i] = data->prng.rand(2.0 * M_PI);
    }

    // tan2 output is [-pi, pi], acos output is [0, pi]
    // Need to determine if we want to use ang or -ang
    // Use phi[0] to compare it
    double phiDiff = CalcOldPhi(newMol, bonded[b], focus) - refPhi;
    phiDiff += (phiDiff < 0.0 ? M_PI : -M_PI);
    bool flip = (phiDiff > 0.0 ? true : false);

    // modify the twist angle if it was fixed or was part of ring
    for (uint c = 0; c < b; ++c) {
      bool angleFix = data->ff.angles->AngleFixed(angleKinds[b][c]);

      if (angleInRing[b][c] || angleFix) {
        double cosTerm = cos(theta[b]) * cos(theta[c]);
        double sinTerm = sin(theta[b]) * sin(theta[c]);
        double bfcRing = CalcTheta(newMol, bonded[b], focus, bonded[c]);
        double var = (cos(bfcRing) - cosTerm) / sinTerm;
        // To fix the numerical problem for flat molecule
        var = (var > 1.0 && var < 1.1 ? 1.0 : var);
        var = (var < -1.0 && var > -1.1 ? -1.0 : var);
        double ang = acos(var);
        ang = (flip ? 2.0 * M_PI - ang : ang);
        ang += phi[c];
        std::fill_n(angles, nTrials, ang);
        if (std::isnan(ang)) {
          std::cout << "Error: Cannot Construct Angle "
                    << newMol.GetKind().atomTypeNames[bonded[b]] << " "
                    << newMol.GetKind().atomTypeNames[focus] << " "
                    << newMol.GetKind().atomTypeNames[bonded[c]] << " !\n";
          if (angleFix) {
            std::cout << "Note: This issue might happened due to defining "
                      << "fix angle.\n";
          }
          exit(EXIT_FAILURE);
        }
        break;
      }
    }

    // compare to angles determined in previous iterations
    for (uint c = 0; c < b; ++c) {
      double cosTerm = cos(theta[b]) * cos(theta[c]);
      double sinTerm = sin(theta[b]) * sin(theta[c]);

      for (uint i = 0; i < nTrials; ++i) {
        double bfcTheta = acos(sinTerm * cos(angles[i] - phi[c]) + cosTerm);
        double distSq =
            newMol.AngleDist(bondLength[b], bondLength[c], bfcTheta);
        double tempEn =
            data->calc.IntraEnergy_1_3(distSq, bonded[b], bonded[c], molIndex);
        nonbonded_1_3[i] += tempEn;
        energies[i] += data->ff.angles->Calc(angleKinds[b][c], bfcTheta);
      }
    }

    // calculate weights from combined energy
    double stepWeight = 0.0;
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(energies, nonbonded_1_3, nTrials, weights) reduction(+:stepWeight)
#endif
    for (int i = 0; i < (int)nTrials; ++i) {
      weights[i] = exp(-1 * data->ff.beta * (energies[i] + nonbonded_1_3[i]));
      stepWeight += weights[i];
    }

    uint winner = data->prng.PickWeighted(weights, nTrials, stepWeight);
    phi[b] = angles[winner];
    bendEnergy += energies[winner];
    oneThree += nonbonded_1_3[winner];
    phiWeight[b] = stepWeight;
  }
}

// Calculate OldMol Bond Energy &
// Calculate phi weight for nTrials using actual theta of OldMol
void DCHedronCycle::ConstrainedAnglesOld(uint nTrials, TrialMol &oldMol,
                                         uint molIndex) {
  double *angles = data->angles;
  IncorporateOld(oldMol, molIndex);

  for (uint b = 1; b < nBonds; ++b) {
    double stepWeight = 0.0;
    // pick "twist" angles
    for (uint i = 0; i < nTrials; ++i) {
      angles[i] = data->prng.rand(2.0 * M_PI);
    }

    // tan2 output is [-pi, pi], acos output is [0, pi]
    // Need to determine if we want to use ang or -ang
    // Use phi[0] to compare it
    double phiDiff = phi[b] - phi[0];
    phiDiff += (phiDiff < 0.0 ? M_PI : -M_PI);
    bool flip = (phiDiff > 0.0 ? true : false);

    // modify the twist angle if it was fixed or was part of ring
    for (uint c = 0; c < b; ++c) {
      bool angleFix = data->ff.angles->AngleFixed(angleKinds[b][c]);

      if (angleInRing[b][c] || angleFix) {
        double cosTerm = cos(theta[b]) * cos(theta[c]);
        double sinTerm = sin(theta[b]) * sin(theta[c]);
        double bfcRing = CalcTheta(oldMol, bonded[b], focus, bonded[c]);
        double var = (cos(bfcRing) - cosTerm) / sinTerm;
        // To fix the numerical problem for flat molecule
        var = (var > 1.0 && var < 1.1 ? 1.0 : var);
        var = (var < -1.0 && var > -1.1 ? -1.0 : var);
        double ang = acos(var);
        ang = (flip ? 2.0 * M_PI - ang : ang);
        ang += phi[c];
        std::fill_n(angles, nTrials, ang);
        if (std::isnan(ang)) {
          std::cout << "Error: Cannot Construct Angle"
                    << oldMol.GetKind().atomTypeNames[bonded[b]] << " "
                    << oldMol.GetKind().atomTypeNames[focus] << " "
                    << oldMol.GetKind().atomTypeNames[bonded[c]] << " !\n";
          if (angleFix) {
            std::cout << "Note: This issue might happened due to defining "
                      << "fix angle.\n";
          }
          exit(EXIT_FAILURE);
        }
        break;
      }
    }

    for (uint i = 0; i < nTrials; ++i) {
      double energies = 0.0;
      double nonbondedEng = 0.0;
      // compare to angles determined in previous iterations
      for (uint c = 0; c < b; ++c) {
        double cosTerm = cos(theta[b]) * cos(theta[c]);
        double sinTerm = sin(theta[b]) * sin(theta[c]);

        double bfcTheta = acos(sinTerm * cos(angles[i] - phi[c]) + cosTerm);
        double distSq =
            oldMol.AngleDist(bondLengthOld[b], bondLengthOld[c], bfcTheta);
        nonbondedEng +=
            data->calc.IntraEnergy_1_3(distSq, bonded[b], bonded[c], molIndex);
        energies += data->ff.angles->Calc(angleKinds[b][c], bfcTheta);
      }
      // calculate weights from combined energy
      double weights = exp(-1 * data->ff.beta * (energies + nonbondedEng));
      stepWeight += weights;
    }
    phiWeight[b] += stepWeight;
  }
}

void DCHedronCycle::SetBasis(TrialMol &mol, uint p1, uint p2) {
  // Calculate the dihedral using bCoords
  const XYZArray &coords = mol.GetBCoords();
  // W is unit vec of p1->p2
  XYZ wVec = coords.Difference(p2, p1);
  wVec.Normalize();
  XYZ uVec;
  // check to make sure our W isn't in line with the standard X Axis
  if (std::abs(wVec.x) < 0.8) {
    // V will be W x the standard X unit vec
    uVec = XYZ(1.0, 0.0, 0.0);
  } else {
    // V will be W x the standard Y unit vec
    uVec = XYZ(0.0, 1.0, 0.0);
  }
  XYZ vVec = geom::Cross(wVec, uVec);
  vVec.Normalize();
  // U is unit vec perpendicular to both V and W
  uVec = geom::Cross(vVec, wVec);
  growthToWorld.BasisRotation(uVec, vVec, wVec);
  worldToGrowth = growthToWorld.Inverse();
}

double DCHedronCycle::CalcOldPhi(TrialMol &mol, uint atom,
                                 uint lastAtom) const {
  // Calculate the dihedral using bCoords
  const XYZArray &coords = mol.GetBCoords();
  XYZ diff = coords.Difference(atom, lastAtom);
  XYZ growthCoords = worldToGrowth.Apply(diff);
  return atan2(growthCoords.y, growthCoords.x);
}

} // namespace cbmc
