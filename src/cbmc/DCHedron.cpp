/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "DCHedron.h"

#include <cassert>
#include <numeric>

#include "DCData.h"
#include "Forcefield.h"
#include "MolSetup.h"
#include "NumLib.h"
#include "PRNG.h"
#include "TrialMol.h"

namespace {
// Wish I could use lambdas..
struct FindA1 {
  FindA1(uint x) : x(x) {};
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

} // namespace

namespace cbmc {

DCHedron::DCHedron(DCData *data, const mol_setup::MolKind &kind, uint focus,
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
  for (uint i = 0; i < nBonds; ++i) {
    typedef std::vector<Angle>::const_iterator Aiter;
    Aiter free =
        find_if(angles.begin(), angles.end(), FindAngle(prev, bonded[i]));
    assert(free != angles.end());
    angleKinds[i][i] = free->kind;

    for (uint j = i + 1; j < nBonds; ++j) {
      Aiter pair = find_if(angles.begin(), angles.end(),
                           FindAngle(bonded[i], bonded[j]));
      angleKinds[i][j] = pair->kind;
      angleKinds[j][i] = pair->kind;
    }
  }

  phi[0] = 0.0;
  phiWeight[0] = 1.0;

  if (data->nAngleTrials < 1) {
    std::cout << "Error: CBMC angle trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

void DCHedron::SetBondNew(double const *bondLen, double const &anchBond) {
  for (uint i = 0; i < nBonds; ++i) {
    bondLength[i] = bondLen[i];
  }
  anchorBond = anchBond;
}

void DCHedron::SetBondOld(double const *bondLen, double const &anchBond) {
  for (uint i = 0; i < nBonds; ++i) {
    bondLengthOld[i] = bondLen[i];
  }
  anchorBondOld = anchBond;
}

double DCHedron::GetWeight() const {
  double result = 1;
  for (uint i = 0; i < nBonds; ++i) {
    result *= thetaWeight[i];
    result *= phiWeight[i];
  }
  return result;
}

// Randomly generate nTrials angles and save energy and weight
void DCHedron::GenerateAnglesNew(TrialMol &newMol, uint molIndex, uint kind,
                                 uint nTrials, uint bType) {
  double *nonbonded_1_3 = data->nonbonded_1_3;
  double thetaFix;
  bool angleFix = false;
  std::fill_n(nonbonded_1_3, nTrials, 0.0);

  if (data->ff.angles->AngleFixed(kind)) {
    angleFix = true;
    thetaFix = data->ff.angles->Angle(kind);
  }

  for (int i = 0; i < (int)nTrials; ++i) {
    if (angleFix)
      data->angles[i] = thetaFix;
    else
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

void DCHedron::GenerateAnglesOld(TrialMol &oldMol, uint molIndex, uint kind,
                                 uint nTrials, uint bType) {
  double *nonbonded_1_3 = data->nonbonded_1_3;
  double thetaFix;
  bool angleFix = false;
  std::fill_n(nonbonded_1_3, nTrials, 0.0);

  if (data->ff.angles->AngleFixed(kind)) {
    angleFix = true;
    thetaFix = data->ff.angles->Angle(kind);
  }

  for (int i = 0; i < (int)nTrials; ++i) {
    if (angleFix)
      data->angles[i] = thetaFix;
    else
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

void DCHedron::FreeAnglesNew(TrialMol &newMol, uint molIndex, uint nTrials) {
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

void DCHedron::FreeAnglesOld(TrialMol &oldMol, uint molIndex, uint nTrials) {
  for (uint i = 0; i < nBonds; ++i) {
    GenerateAnglesOld(oldMol, molIndex, angleKinds[i][i], nTrials, i);
    double stepWeight =
        std::accumulate(data->angleWeights, data->angleWeights + nTrials, 0.0);
    thetaWeight[i] = stepWeight;
  }
}

void DCHedron::PrepareNew(TrialMol &newMol, uint molIndex) {
  bendEnergy = 0.0;
  oneThree = 0.0;
  FreeAnglesNew(newMol, molIndex, data->nAngleTrials);
  ConstrainedAngles(newMol, molIndex, data->nAngleTrials);
}

void DCHedron::PrepareOld(TrialMol &oldMol, uint molIndex) {
  oneThree = 0.0;
  bendEnergy = 0.0;
  FreeAnglesOld(oldMol, molIndex, data->nAngleTrials - 1);
}

void DCHedron::IncorporateOld(TrialMol &oldMol, uint molIndex) {
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

void DCHedron::ConstrainedAngles(TrialMol &newMol, uint molIndex,
                                 uint nTrials) {
  double *angles = data->angles;
  double *energies = data->angleEnergy;
  double *weights = data->angleWeights;
  double *nonbonded_1_3 = data->nonbonded_1_3;
  phi[0] = 0.0;

  for (uint b = 1; b < nBonds; ++b) {
    std::fill_n(energies, nTrials, 0.0);
    std::fill_n(nonbonded_1_3, nTrials, 0.0);
    // pick "twist" angles
    for (uint i = 0; i < nTrials; ++i) {
      angles[i] = data->prng.rand(2.0 * M_PI);
    }

    // modify the twist angle if it was fixed
    for (uint c = 0; c < b; ++c) {
      double cosTerm = cos(theta[b]) * cos(theta[c]);
      double sinTerm = sin(theta[b]) * sin(theta[c]);

      if (data->ff.angles->AngleFixed(angleKinds[b][c])) {
        double bfcTheta = data->ff.angles->Angle(angleKinds[b][c]);
        double var = (cos(bfcTheta) - cosTerm) / sinTerm;
        // To fix the numerical problem for flat molecule
        var = (var > 1.0 && var < 1.1 ? 1.0 : var);
        var = (var < -1.0 && var > -1.1 ? -1.0 : var);
        double ang = acos(var) + phi[c];
        std::fill_n(angles, nTrials, ang);
        if (std::isnan(ang)) {
          // printf("Val: %2.10f, angle: %2.5f \n", var, ang);
          std::cout << "Error: Cannot constrain fix angle for "
                    << newMol.GetKind().atomTypeNames[bonded[b]] << " "
                    << newMol.GetKind().atomTypeNames[focus] << " "
                    << newMol.GetKind().atomTypeNames[bonded[c]] << " !\n";
          exit(EXIT_FAILURE);
        }
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
#pragma omp parallel for default(none)                                         \
    shared(energies, nonbonded_1_3, nTrials, weights)                          \
    reduction(+ : stepWeight)
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
void DCHedron::ConstrainedAnglesOld(uint nTrials, TrialMol &oldMol,
                                    uint molIndex) {
  double *angles = data->angles;
  IncorporateOld(oldMol, molIndex);

  for (uint b = 1; b < nBonds; ++b) {
    double stepWeight = 0.0;

    // pick "twist" angles
    for (uint i = 0; i < nTrials; ++i) {
      angles[i] = data->prng.rand(2.0 * M_PI);
    }

    // modify the twist angle if it was fixed
    for (uint c = 0; c < b; ++c) {
      double cosTerm = cos(theta[b]) * cos(theta[c]);
      double sinTerm = sin(theta[b]) * sin(theta[c]);

      if (data->ff.angles->AngleFixed(angleKinds[b][c])) {
        double bfcTheta = data->ff.angles->Angle(angleKinds[b][c]);
        double var = (cos(bfcTheta) - cosTerm) / sinTerm;
        // To fix the numerical problem for flat molecule
        var = (var > 1.0 && var < 1.1 ? 1.0 : var);
        var = (var < -1.0 && var > -1.1 ? -1.0 : var);
        double ang = acos(var) + phi[c];
        std::fill_n(angles, nTrials, ang);
        if (std::isnan(ang)) {
          // printf("Val: %2.10f, angle: %2.5f \n", var, ang);
          std::cout << "Error: Cannot constrain fix angle for "
                    << oldMol.GetKind().atomTypeNames[bonded[b]] << " "
                    << oldMol.GetKind().atomTypeNames[focus] << " "
                    << oldMol.GetKind().atomTypeNames[bonded[c]] << " !\n";
          exit(EXIT_FAILURE);
        }
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
} // namespace cbmc
