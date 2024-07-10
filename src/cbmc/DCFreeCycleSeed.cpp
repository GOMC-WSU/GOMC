/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#define _USE_MATH_DEFINES
#include "DCFreeCycleSeed.h"

#include <cmath>

#include "DCData.h"
#include "Forcefield.h"
#include "MolSetup.h"
#include "NumLib.h"
#include "PRNG.h"
#include "TrialMol.h"

namespace cbmc {

struct FindA1 {
  FindA1(uint x) : x(x){};
  bool operator()(const mol_setup::Bond &b) { return (b.a1 == x); }
  uint x;
};

DCFreeCycleSeed::DCFreeCycleSeed(DCData *data, const mol_setup::MolKind &kind,
                                 const std::vector<int> &cycAtoms, uint focus,
                                 uint prev)
    : data(data), hed(data, kind, cycAtoms, focus, prev) {
  using namespace mol_setup;
  std::fill_n(bondedInRing, MAX_BONDS, false);
  std::vector<Bond> onFocus = AtomBonds(kind, hed.Focus());
  for (uint i = 0; i < onFocus.size(); ++i) {
    if (onFocus[i].a1 == prev) {
      anchorKind = onFocus[i].kind;
      if (std::find(cycAtoms.begin(), cycAtoms.end(), onFocus[i].a1) !=
          cycAtoms.end()) {
        bondedInRing[hed.NumBond()] = true;
      }
      break;
    }
  }

  onFocus.erase(remove_if(onFocus.begin(), onFocus.end(), FindA1(prev)),
                onFocus.end());
  // Find the atoms bonded to focus, except prev
  for (uint i = 0; i < onFocus.size(); ++i) {
    bondKinds[i] = onFocus[i].kind;
    if (std::find(cycAtoms.begin(), cycAtoms.end(), onFocus[i].a1) !=
        cycAtoms.end()) {
      bondedInRing[i] = true;
    }
  }

  if (data->nLJTrialsNth < 1) {
    std::cout << "Error: CBMC secondary atom trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

void DCFreeCycleSeed::PrepareNew(TrialMol &newMol, uint molIndex) {
  // Get new bond information
  SetBondLengthNew(newMol);
  hed.SetBondNew(bondLength, anchorBond);
  hed.PrepareNew(newMol, molIndex);
  bondEnergy = 0.0;
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondEnergy += data->ff.bonds.Calc(bondKinds[i], bondLength[i]);
  }
  bondEnergy += data->ff.bonds.Calc(anchorKind, anchorBond);
}

void DCFreeCycleSeed::PrepareOld(TrialMol &oldMol, uint molIndex) {
  // Get old bond information
  SetBondLengthOld(oldMol);
  hed.SetBondOld(bondLengthOld, anchorBondOld);
  hed.PrepareOld(oldMol, molIndex);
  bondEnergy = 0.0;
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondEnergy += data->ff.bonds.Calc(bondKinds[i], bondLengthOld[i]);
  }
  bondEnergy += data->ff.bonds.Calc(anchorKind, anchorBondOld);
}

void DCFreeCycleSeed::SetBondLengthNew(TrialMol &newMol) {
  for (uint i = 0; i < hed.NumBond(); ++i) {
    // use bondLength from bCoords if it is in the ring body
    if (bondedInRing[i]) {
      bondLength[i] =
          newMol.GetBCoords().Difference(hed.Bonded(i), hed.Focus()).Length();
    } else {
      bondLength[i] = data->ff.bonds.Length(bondKinds[i]);
    }
  }
  if (bondedInRing[hed.NumBond()]) {
    anchorBond =
        newMol.GetBCoords().Difference(hed.Focus(), hed.Prev()).Length();
  } else {
    anchorBond = data->ff.bonds.Length(anchorKind);
  }
}

void DCFreeCycleSeed::SetBondLengthOld(TrialMol &oldMol) {
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondLengthOld[i] = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Bonded(i)));
  }
  anchorBondOld = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Prev()));
}

void DCFreeCycleSeed::BuildNew(TrialMol &newMol, uint molIndex) {
  PRNG &prng = data->prng;
  const CalculateEnergy &calc = data->calc;
  // const Ewald *calcEwald = data->calcEwald;
  const Forcefield &ff = data->ff;
  uint nLJTrials = data->nLJTrialsNth;
  double *ljWeights = data->ljWeights;
  double *inter = data->inter;
  double *real = data->real;
  bool *overlap = data->overlap;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  // get info about existing geometry
  newMol.ShiftBasis(hed.Focus());
  const XYZ center = newMol.AtomPosition(hed.Focus());
  XYZArray *positions = data->multiPositions;

  for (uint i = 0; i < hed.NumBond(); ++i) {
    positions[i].Set(
        0, newMol.RawRectCoords(bondLength[i], hed.Theta(i), hed.Phi(i)));
  }
  // add anchor atom
  positions[hed.NumBond()].Set(0, newMol.RawRectCoords(anchorBond, 0, 0));

  // counting backward to preserve prototype
  double u1, u2, u3;
  for (uint lj = nLJTrials; lj-- > 0;) {
    // convert chosen torsion to 3D positions
    u1 = prng();
    u2 = prng();
    u3 = prng();
    RotationMatrix spin = RotationMatrix::UniformRandom(u1, u2, u3);
    // RotationMatrix spin = RotationMatrix::UniformRandom(prng(), prng(), prng());
    for (uint b = 0; b < hed.NumBond() + 1; ++b) {
      // find positions
      positions[b].Set(lj, spin.Apply(positions[b][0]));
      positions[b].Add(lj, center);
    }
  }

  for (uint b = 0; b < hed.NumBond() + 1; ++b) {
    data->axes.WrapPBC(positions[b], newMol.GetBox());
  }

  for (uint b = 0; b < hed.NumBond(); ++b) {
    calc.ParticleInter(inter, real, positions[b], overlap, hed.Bonded(b),
                       molIndex, newMol.GetBox(), nLJTrials);
  }
  calc.ParticleInter(inter, real, positions[hed.NumBond()], overlap, hed.Prev(),
                     molIndex, newMol.GetBox(), nLJTrials);

  double stepWeight = 0;
  for (uint lj = 0; lj < nLJTrials; ++lj) {
    ljWeights[lj] = exp(-ff.beta * (inter[lj] + real[lj]));
    stepWeight += ljWeights[lj];
  }

  uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  for (uint b = 0; b < hed.NumBond(); ++b) {
    newMol.AddAtom(hed.Bonded(b), positions[b][winner]);
    newMol.AddBonds(hed.Bonded(b), hed.Focus());
  }

  newMol.AddAtom(hed.Prev(), positions[hed.NumBond()][winner]);
  newMol.AddBonds(hed.Prev(), hed.Focus());
  newMol.UpdateOverlap(overlap[winner]);
  newMol.AddEnergy(Energy(hed.GetEnergy() + bondEnergy, hed.GetNonBondedEn(),
                          inter[winner], real[winner], 0.0, 0.0, 0.0));
  newMol.MultWeight(hed.GetWeight());
  newMol.MultWeight(stepWeight / nLJTrials);
}

void DCFreeCycleSeed::BuildOld(TrialMol &oldMol, uint molIndex) {
  PRNG &prng = data->prng;
  const CalculateEnergy &calc = data->calc;
  // const Ewald * calcEwald = data->calcEwald;
  const Forcefield &ff = data->ff;
  uint nLJTrials = data->nLJTrialsNth;
  double *ljWeights = data->ljWeights;
  double *inter = data->inter;
  double *real = data->real;
  bool *overlap = data->overlap;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  // get info about existing geometry
  oldMol.SetBasis(hed.Focus(), hed.Prev());
  // Calculate OldMol Bond Energy &
  // Calculate phi weight for nTrials using actual theta of OldMol
  hed.ConstrainedAnglesOld(data->nAngleTrials - 1, oldMol, molIndex);
  const XYZ center = oldMol.AtomPosition(hed.Focus());
  XYZArray *positions = data->multiPositions;
  for (uint i = 0; i < hed.NumBond(); ++i) {
    // get position and shift to origin
    positions[i].Set(0, oldMol.AtomPosition(hed.Bonded(i)));
    data->axes.UnwrapPBC(positions[i], 0, 1, oldMol.GetBox(), center);
    positions[i].Add(0, -center);
  }

  // add anchor atom
  positions[hed.NumBond()].Set(0, oldMol.AtomPosition(hed.Prev()));
  data->axes.UnwrapPBC(positions[hed.NumBond()], 0, 1, oldMol.GetBox(), center);
  positions[hed.NumBond()].Add(0, -center);

  // counting backward to preserve prototype
  double u1, u2, u3;
  for (uint lj = nLJTrials; lj-- > 1;) {
    // convert chosen torsion to 3D positions
    u1 = prng();
    u2 = prng();
    u3 = prng();
    RotationMatrix spin = RotationMatrix::UniformRandom(u1, u2, u3);
    // RotationMatrix spin = RotationMatrix::UniformRandom(prng(), prng(), prng());
    for (uint b = 0; b < hed.NumBond() + 1; ++b) {
      // find positions
      positions[b].Set(lj, spin.Apply(positions[b][0]));
      positions[b].Add(lj, center);
    }
  }

  for (uint b = 0; b < hed.NumBond() + 1; ++b) {
    positions[b].Add(0, center);
    data->axes.WrapPBC(positions[b], oldMol.GetBox());
  }

  for (uint b = 0; b < hed.NumBond(); ++b) {
    calc.ParticleInter(inter, real, positions[b], overlap, hed.Bonded(b),
                       molIndex, oldMol.GetBox(), nLJTrials);
  }
  double stepWeight = 0;
  calc.ParticleInter(inter, real, positions[hed.NumBond()], overlap, hed.Prev(),
                     molIndex, oldMol.GetBox(), nLJTrials);

  for (uint lj = 0; lj < nLJTrials; ++lj) {
    stepWeight += exp(-ff.beta * (inter[lj] + real[lj]));
  }

  for (uint b = 0; b < hed.NumBond(); ++b) {
    oldMol.ConfirmOldAtom(hed.Bonded(b));
    oldMol.AddBonds(hed.Bonded(b), hed.Focus());
  }

  oldMol.ConfirmOldAtom(hed.Prev());
  oldMol.AddBonds(hed.Prev(), hed.Focus());
  oldMol.UpdateOverlap(overlap[0]);
  oldMol.AddEnergy(Energy(hed.GetEnergy() + bondEnergy, hed.GetNonBondedEn(),
                          inter[0], real[0], 0.0, 0.0, 0.0));
  oldMol.MultWeight(hed.GetWeight());
  oldMol.MultWeight(stepWeight / nLJTrials);
}

} // namespace cbmc
