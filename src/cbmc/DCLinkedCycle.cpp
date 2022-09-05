/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#define _USE_MATH_DEFINES
#include "DCLinkedCycle.h"

#include <cassert>
#include <cmath>
#include <numeric>

#include "DCData.h"
#include "Forcefield.h"
#include "MolSetup.h"
#include "NumLib.h"
#include "PRNG.h"
#include "TrialMol.h"

namespace {
struct FindA1 {
  FindA1(uint x) : x(x){};
  bool operator()(const mol_setup::Bond &b) { return (b.a1 == x); }
  uint x;
};

struct FindDih {
  FindDih(uint x, uint y) : x(x), y(y) {}
  uint x, y;
  bool operator()(const mol_setup::Dihedral d) {
    return (d.a0 == x && d.a3 == y) || (d.a0 == y && d.a3 == x);
  }
};

} // namespace

namespace cbmc {
DCLinkedCycle::DCLinkedCycle(DCData *data, const mol_setup::MolKind &kind,
                             std::vector<int> cycAtoms, uint focus, uint prev)
    : data(data), hed(data, kind, cycAtoms, focus, prev) {
  using namespace mol_setup;
  std::fill_n(bondedInRing, MAX_BONDS, false);
  std::vector<Bond> onFocus = AtomBonds(kind, hed.Focus());
  onFocus.erase(remove_if(onFocus.begin(), onFocus.end(), FindA1(prev)),
                onFocus.end());
  // Find the atoms bonded to focus, except prev
  focBondedRing = -1;
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondKinds[i] = onFocus[i].kind;
    // store the dihedral that boded[i] and focus are in the middle
    bondedFocusDih.push_back(DihsOnBond(kind, onFocus[i].a1, focus));
    if (std::find(cycAtoms.begin(), cycAtoms.end(), onFocus[i].a1) !=
        cycAtoms.end()) {
      focBondedRing = i;
      bondedInRing[i] = true;
    }
  }
  assert(focBondedRing != -1);

  std::vector<Bond> onPrev = AtomBonds(kind, hed.Prev());
  onPrev.erase(remove_if(onPrev.begin(), onPrev.end(), FindA1(hed.Focus())),
               onPrev.end());
  nPrevBonds = onPrev.size();
  prevBondedRing = -1;
  for (uint i = 0; i < nPrevBonds; ++i) {
    prevBonded[i] = onPrev[i].a1;
    if (std::find(cycAtoms.begin(), cycAtoms.end(), onPrev[i].a1) !=
        cycAtoms.end()) {
      prevBondedRing = i;
    }
  }

  std::vector<Dihedral> dihs = DihsOnBond(kind, hed.Focus(), hed.Prev());
  for (uint i = 0; i < hed.NumBond(); ++i) {
    for (uint j = 0; j < nPrevBonds; ++j) {
      std::vector<Dihedral>::const_iterator match = find_if(
          dihs.begin(), dihs.end(), FindDih(hed.Bonded(i), prevBonded[j]));
      assert(match != dihs.end());
      dihKinds[i][j] = match->kind;
    }
  }

  if (data->nLJTrialsNth < 1) {
    std::cout << "Error: CBMC secondary atom trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }

  if (data->nDihTrials < 1) {
    std::cout << "Error: CBMC dihedral trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

double DCLinkedCycle::CalcDih(TrialMol &mol, uint a0, uint a1, uint a2,
                              uint a3) {
  double phi = 0.0;
  if (mol.AtomExists(a0)) {
    // Calculate theta using tCoords
    phi = mol.GetPhi(a0, a1, a2, a3);
  } else {
    // Calculate the dihedral using bCoords
    const XYZArray &coords = mol.GetBCoords();
    phi = geom::Phi(coords.Difference(a1, a0), coords.Difference(a2, a1),
                    coords.Difference(a3, a2));
  }
  return phi;
}

void DCLinkedCycle::PrepareNew(TrialMol &newMol, uint molIndex) {
  // Get new bond information
  SetBondLengthNew(newMol);
  hed.SetBondNew(bondLength, anchorBond);
  hed.PrepareNew(newMol, molIndex);
  bondEnergy = 0.0;
  bExist.resize(hed.NumBond(), false);
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bExist[i] = newMol.BondsExist(hed.Bonded(i), hed.Focus());
    if (!bExist[i]) {
      bondEnergy += data->ff.bonds.Calc(bondKinds[i], bondLength[i]);
    }
  }
}

void DCLinkedCycle::PrepareOld(TrialMol &oldMol, uint molIndex) {
  // Get old bond information
  SetBondLengthOld(oldMol);
  hed.SetBondOld(bondLengthOld, anchorBondOld);
  hed.PrepareOld(oldMol, molIndex);
  bondEnergy = 0.0;
  bExist.resize(hed.NumBond(), false);
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bExist[i] = oldMol.BondsExist(hed.Bonded(i), hed.Focus());
    if (!bExist[i]) {
      bondEnergy += data->ff.bonds.Calc(bondKinds[i], bondLengthOld[i]);
    }
  }
}

void DCLinkedCycle::SetBondLengthNew(TrialMol &newMol) {
  for (uint i = 0; i < hed.NumBond(); ++i) {
    // use bondLength from bCoords if it is in the ring body
    if (bondedInRing[i]) {
      bondLength[i] =
          newMol.GetBCoords().Difference(hed.Bonded(i), hed.Focus()).Length();
    } else {
      bondLength[i] = data->ff.bonds.Length(bondKinds[i]);
    }
  }
  // anchorBond is built, we need the actual length
  anchorBond = sqrt(newMol.OldDistSq(hed.Focus(), hed.Prev()));
}

void DCLinkedCycle::SetBondLengthOld(TrialMol &oldMol) {
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondLengthOld[i] = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Bonded(i)));
  }
  anchorBondOld = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Prev()));
}

void DCLinkedCycle::BuildNew(TrialMol &newMol, uint molIndex) {
  PRNG &prng = data->prng;
  // const CalculateEnergy& calc = data->calc;
  // const Forcefield& ff = data->ff;
  double *torsion = data->angles;
  double *torWeights = data->angleWeights;
  double *torEnergy = data->angleEnergy;
  double *ljWeights = data->ljWeights;
  double *bondedEn = data->bonded;
  double *inter = data->inter;
  double *nonbonded = data->nonbonded;
  double *nonbonded_1_4 = data->nonbonded_1_4;
  double *real = data->real;
  double *oneFour = data->oneFour;
  bool *overlap = data->overlap;
  uint nDihTrials = data->nDihTrials;
  uint nLJTrials;
  if (prevBondedRing == -1) {
    // we can perform trial
    nLJTrials = data->nLJTrialsNth;
  } else {
    // it is in the body of ring, no trial
    nLJTrials = 1;
  }

  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(bondedEn, nLJTrials, 0.0);
  std::fill_n(oneFour, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  // get info about existing geometry
  newMol.SetBasis(hed.Focus(), hed.Prev());
  const XYZ center = newMol.AtomPosition(hed.Focus());
  XYZArray *positions = data->multiPositions;
  double prevPhi[MAX_BONDS];
  for (uint i = 0; i < hed.NumBond(); ++i) {
    // get position and shift to origin
    positions[i].Set(
        0, newMol.RawRectCoords(bondLength[i], hed.Theta(i), hed.Phi(i)));
  }
  for (uint i = 0; i < nPrevBonds; ++i) {
    double th;
    // not using theta, so this is a wasted cos and sqrt
    newMol.OldThetaAndPhi(prevBonded[i], hed.Prev(), th, prevPhi[i]);
  }
  if (prevBondedRing != -1) {
    // calculate dihedral in ring using bCoords
    double phi = CalcDih(newMol, hed.Bonded(focBondedRing), hed.Focus(),
                         hed.Prev(), prevBonded[prevBondedRing]);
    // find the torsion that give the same dihedral in the ring
    torDiff = phi - (hed.Phi(focBondedRing) - prevPhi[prevBondedRing]);
  }
  XYZ rotationAxis =
      newMol.AtomPosition(hed.Focus()) - newMol.AtomPosition(hed.Prev());
  rotationAxis = data->axes.MinImage(rotationAxis, newMol.GetBox());
  rotationAxis *= (1 / rotationAxis.Length());
  RotationMatrix cross = RotationMatrix::CrossProduct(rotationAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotationAxis);

  // counting backward to preserve prototype
  for (uint lj = nLJTrials; lj-- > 0;) {
    ChooseTorsion(newMol, molIndex, prevPhi, cross, tensor);
    ljWeights[lj] = std::accumulate(torWeights, torWeights + nDihTrials, 0.0);
    uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
    bondedEn[lj] = torEnergy[winner];
    oneFour[lj] = nonbonded_1_4[winner];
    // convert chosen torsion to 3D positions
    RotationMatrix spin =
        RotationMatrix::FromAxisAngle(-torsion[winner], cross, tensor);
    for (uint b = 0; b < hed.NumBond(); ++b) {
      // find positions
      positions[b].Set(lj, spin.Apply(positions[b][0]));
      positions[b].Add(lj, center);
    }
  }

  for (uint b = 0; b < hed.NumBond(); ++b) {
    data->axes.WrapPBC(positions[b], newMol.GetBox());
  }

  double stepWeight = EvalLJ(newMol, molIndex);
  uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);

  for (uint b = 0; b < hed.NumBond(); ++b) {
    if (!newMol.AtomExists(hed.Bonded(b)))
      newMol.AddAtom(hed.Bonded(b), positions[b][winner]);
    newMol.AddBonds(hed.Bonded(b), hed.Focus());
  }

  // In ring, we will miss the dihedral for existed bonds
  // Calculate the dihedral and 1-4 interaction for existed bonds
  for (uint b = 0; b < hed.NumBond(); ++b) {
    if (bExist[b]) {
      CaclIntraEnergy(newMol, b, molIndex);
    }
  }

  newMol.UpdateOverlap(overlap[winner]);
  newMol.AddEnergy(
      Energy(bondedEn[winner] + hed.GetEnergy() + bondEnergy,
             nonbonded[winner] + hed.GetNonBondedEn() + oneFour[winner],
             inter[winner], real[winner], 0.0, 0.0, 0.0));
  newMol.MultWeight(hed.GetWeight());
  newMol.MultWeight(stepWeight / nLJTrials);
}

void DCLinkedCycle::BuildOld(TrialMol &oldMol, uint molIndex) {
  PRNG &prng = data->prng;
  // const CalculateEnergy& calc = data->calc;
  const Forcefield &ff = data->ff;
  double *torsion = data->angles;
  double *torWeights = data->angleWeights;
  double *torEnergy = data->angleEnergy;
  double *ljWeights = data->ljWeights;
  double *bondedEn = data->bonded;
  double *inter = data->inter;
  double *nonbonded = data->nonbonded;
  double *nonbonded_1_4 = data->nonbonded_1_4;
  double *real = data->real;
  double *oneFour = data->oneFour;
  bool *overlap = data->overlap;
  uint nDihTrials = data->nDihTrials;
  uint nLJTrials;
  if (prevBondedRing == -1) {
    // we can perform trial
    nLJTrials = data->nLJTrialsNth;
  } else {
    // it is in the body of ring, no trial
    nLJTrials = 1;
  }

  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(bondedEn, nLJTrials, 0.0);
  std::fill_n(oneFour, nLJTrials, 0.0);
  std::fill_n(oneFour, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  // get info about existing geometry
  oldMol.SetBasis(hed.Focus(), hed.Prev());
  // Calculate OldMol Bond Energy &
  // Calculate phi weight for nTrials using actual theta of OldMol
  hed.ConstrainedAnglesOld(data->nAngleTrials - 1, oldMol, molIndex);
  const XYZ center = oldMol.AtomPosition(hed.Focus());
  XYZArray *positions = data->multiPositions;
  double prevPhi[MAX_BONDS];
  for (uint i = 0; i < hed.NumBond(); ++i) {
    // get position and shift to origin
    positions[i].Set(0, oldMol.AtomPosition(hed.Bonded(i)));
    data->axes.UnwrapPBC(positions[i], 0, 1, oldMol.GetBox(), center);
    positions[i].Add(0, -center);
  }
  for (uint i = 0; i < nPrevBonds; ++i) {
    double t;
    // not using theta, so this is a wasted cos and sqrt
    oldMol.OldThetaAndPhi(prevBonded[i], hed.Prev(), t, prevPhi[i]);
  }
  if (prevBondedRing != -1) {
    // calculate dihedral in ring using bCoords
    double phi = CalcDih(oldMol, hed.Bonded(focBondedRing), hed.Focus(),
                         hed.Prev(), prevBonded[prevBondedRing]);
    // find the torsion that give the same dihedral in the ring
    torDiff = phi - (hed.Phi(focBondedRing) - prevPhi[prevBondedRing]);
  }
  XYZ rotationAxis =
      oldMol.AtomPosition(hed.Focus()) - oldMol.AtomPosition(hed.Prev());
  rotationAxis = data->axes.MinImage(rotationAxis, oldMol.GetBox());
  rotationAxis *= (1 / rotationAxis.Length());
  RotationMatrix cross = RotationMatrix::CrossProduct(rotationAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotationAxis);

  // counting backward to preserve prototype
  for (uint lj = nLJTrials; lj-- > 1;) {
    ChooseTorsion(oldMol, molIndex, prevPhi, cross, tensor);
    ljWeights[lj] = std::accumulate(torWeights, torWeights + nDihTrials, 0.0);
    uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
    bondedEn[lj] = torEnergy[winner];
    oneFour[lj] = nonbonded_1_4[winner];
    // convert chosen torsion to 3D positions
    RotationMatrix spin =
        RotationMatrix::FromAxisAngle(-torsion[winner], cross, tensor);
    for (uint b = 0; b < hed.NumBond(); ++b) {
      // find positions
      positions[b].Set(lj, spin.Apply(positions[b][0]));
      positions[b].Add(lj, center);
    }
  }
  // for actual atom position, we perform nDihTrials - 1 dihedral trial
  ljWeights[0] = 0.0;
  for (uint tor = 0; tor < nDihTrials; ++tor) {
    // No trial torsion if it is not free end
    if (prevBondedRing == -1) {
      torsion[tor] = (tor == 0) ? 0.0 : data->prng.rand(2.0 * M_PI);
    } else {
      torsion[tor] = 0.0;
    }
    torEnergy[tor] = 0.0;
    nonbonded_1_4[tor] = 0.0;
    for (uint b = 0; b < hed.NumBond(); ++b) {
      double trialPhi = hed.Phi(b) + torsion[tor];
      XYZ bondedC;
      if (oldMol.OneFour()) {
        // convert chosen torsion to 3D positions for bonded atoms to focus
        RotationMatrix spin =
            RotationMatrix::FromAxisAngle(-torsion[tor], cross, tensor);
        bondedC = spin.Apply(positions[b][0]) + center;
      }

      for (uint p = 0; p < nPrevBonds; ++p) {
        if (oldMol.OneFour()) {
          double distSq =
              oldMol.DistSq(bondedC, oldMol.AtomPosition(prevBonded[p]));
          nonbonded_1_4[tor] += data->calc.IntraEnergy_1_4(
              distSq, prevBonded[p], hed.Bonded(b), molIndex);
          if (std::isnan(nonbonded_1_4[tor]))
            nonbonded_1_4[tor] = num::BIGNUM;
        }
        torEnergy[tor] +=
            ff.dihedrals.Calc(dihKinds[b][p], trialPhi - prevPhi[p]);
      }
    }
    ljWeights[0] += exp(-ff.beta * (torEnergy[tor] + nonbonded_1_4[tor]));
  }
  bondedEn[0] = torEnergy[0];
  oneFour[0] = nonbonded_1_4[0];

  for (uint b = 0; b < hed.NumBond(); ++b) {
    positions[b].Add(0, center);
    data->axes.WrapPBC(positions[b], oldMol.GetBox());
  }
  double stepWeight = EvalLJ(oldMol, molIndex);

  for (uint b = 0; b < hed.NumBond(); ++b) {
    oldMol.ConfirmOldAtom(hed.Bonded(b));
    oldMol.AddBonds(hed.Bonded(b), hed.Focus());
  }

  // In ring, we will miss the dihedral for existed bonds
  // Calculate the dihedral and 1-4 interaction for existed bonds
  for (uint b = 0; b < hed.NumBond(); ++b) {
    if (bExist[b]) {
      CaclIntraEnergy(oldMol, b, molIndex);
    }
  }

  oldMol.UpdateOverlap(overlap[0]);
  oldMol.AddEnergy(Energy(bondedEn[0] + hed.GetEnergy() + bondEnergy,
                          nonbonded[0] + hed.GetNonBondedEn() + oneFour[0],
                          inter[0], real[0], 0.0, 0.0, 0.0));

  oldMol.MultWeight(hed.GetWeight());
  oldMol.MultWeight(stepWeight / nLJTrials);
}

double DCLinkedCycle::EvalLJ(TrialMol &mol, uint molIndex) {
  double *inter = data->inter;
  double *nonbonded = data->nonbonded;
  double *real = data->real;
  bool *overlap = data->overlap;
  XYZArray *positions = data->multiPositions;
  uint nLJTrials;
  if (prevBondedRing == -1) {
    // we can perform trial
    nLJTrials = data->nLJTrialsNth;
  } else {
    // it is in the body of ring, no trial
    nLJTrials = 1;
  }

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(nonbonded, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);

  for (uint b = 0; b < hed.NumBond(); ++b) {
    // Avoid double calculating energy for existed atom
    if (mol.AtomExists(hed.Bonded(b))) {
      continue;
    }
    data->calc.ParticleInter(inter, real, positions[b], overlap, hed.Bonded(b),
                             molIndex, mol.GetBox(), nLJTrials);

    data->calc.ParticleNonbonded(nonbonded, mol, positions[b], hed.Bonded(b),
                                 mol.GetBox(), nLJTrials);
  }
  double stepWeight = 0;
  for (uint lj = 0; lj < nLJTrials; ++lj) {
    data->ljWeights[lj] *=
        exp(-data->ff.beta * (inter[lj] + nonbonded[lj] + real[lj]));
    stepWeight += data->ljWeights[lj];
  }
  return stepWeight;
}

void DCLinkedCycle::ChooseTorsion(TrialMol &mol, uint molIndex,
                                  double prevPhi[], RotationMatrix &cross,
                                  RotationMatrix &tensor) {
  double *torsion = data->angles;
  double *torEnergy = data->angleEnergy;
  double *torWeights = data->angleWeights;
  double *nonbonded_1_4 = data->nonbonded_1_4;
  uint nDihTrials = data->nDihTrials;
  const Forcefield &ff = data->ff;
  // To get the information if initial posotion before applying torsion
  XYZArray *positions = data->multiPositions;

  std::fill_n(torsion, data->nDihTrials, 0.0);
  std::fill_n(torWeights, data->nDihTrials, 0.0);
  std::fill_n(torEnergy, data->nDihTrials, 0.0);
  std::fill_n(nonbonded_1_4, data->nDihTrials, 0.0);

  const XYZ center = mol.AtomPosition(hed.Focus());
  // select torsion based on all dihedral angles
  for (uint tor = 0; tor < nDihTrials; ++tor) {
    if (prevBondedRing != -1) {
      torsion[tor] = torDiff;
    } else {
      torsion[tor] = data->prng.rand(2.0 * M_PI);
    }
    torEnergy[tor] = 0.0;
    nonbonded_1_4[tor] = 0.0;
    for (uint b = 0; b < hed.NumBond(); ++b) {
      double trialPhi = hed.Phi(b) + torsion[tor];
      XYZ bondedC;
      if (mol.OneFour()) {
        // convert chosen torsion to 3D positions for bonded atoms to focus
        RotationMatrix spin =
            RotationMatrix::FromAxisAngle(-torsion[tor], cross, tensor);
        bondedC = spin.Apply(positions[b][0]) + center;
      }

      for (uint p = 0; p < nPrevBonds; ++p) {
        if (mol.OneFour()) {
          double distSq = mol.DistSq(bondedC, mol.AtomPosition(prevBonded[p]));
          nonbonded_1_4[tor] += data->calc.IntraEnergy_1_4(
              distSq, prevBonded[p], hed.Bonded(b), molIndex);
          if (std::isnan(nonbonded_1_4[tor]))
            nonbonded_1_4[tor] = num::BIGNUM;
        }

        torEnergy[tor] +=
            ff.dihedrals.Calc(dihKinds[b][p], trialPhi - prevPhi[p]);
      }
    }
    torWeights[tor] = exp(-ff.beta * (torEnergy[tor] + nonbonded_1_4[tor]));
  }
}

void DCLinkedCycle::CaclIntraEnergy(TrialMol &mol, const uint bIdx,
                                    const uint molIndex) {
  double dihEnergy = 0.0;
  double nonBondedEn = 0.0;
  uint size = bondedFocusDih[bIdx].size();
  const std::vector<mol_setup::Dihedral> &dih = bondedFocusDih[bIdx];
  for (uint d = 0; d < size; d++) {
    // We know that a1(bonded[i]) and a2(focus) already exist
    if (mol.AtomExists(dih[d].a0) && mol.AtomExists(dih[d].a3)) {
      double phi = mol.GetPhi(dih[d].a0, dih[d].a1, dih[d].a2, dih[d].a3);
      dihEnergy += data->ff.dihedrals.Calc(dih[d].kind, phi);
      if (mol.OneFour()) {
        double distSq = mol.DistSq(mol.AtomPosition(dih[d].a0),
                                   mol.AtomPosition(dih[d].a3));
        nonBondedEn +=
            data->calc.IntraEnergy_1_4(distSq, dih[d].a0, dih[d].a3, molIndex);
      }
    }
  }
  double weight = exp(-1.0 * data->ff.beta * (dihEnergy + nonBondedEn));
  mol.MultWeight(weight);
  mol.AddEnergy(Energy(dihEnergy, nonBondedEn, 0.0, 0.0, 0.0, 0.0, 0.0));
}

} // namespace cbmc
