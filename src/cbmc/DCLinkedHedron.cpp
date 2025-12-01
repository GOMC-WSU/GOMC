/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "DCLinkedHedron.h"

#include <cassert>
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
DCLinkedHedron::DCLinkedHedron(DCData *data, const mol_setup::MolKind &kind,
                               uint focus, uint prev)
    : data(data), hed(data, kind, focus, prev) {
  using namespace mol_setup;
  std::vector<Bond> onFocus = AtomBonds(kind, hed.Focus());
  onFocus.erase(remove_if(onFocus.begin(), onFocus.end(), FindA1(prev)),
                onFocus.end());
  // Find the atoms bonded to focus, except prev
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondKinds[i] = onFocus[i].kind;
  }

  std::vector<Bond> onPrev = AtomBonds(kind, hed.Prev());
  onPrev.erase(remove_if(onPrev.begin(), onPrev.end(), FindA1(hed.Focus())),
               onPrev.end());
  nPrevBonds = onPrev.size();

  for (uint i = 0; i < nPrevBonds; ++i) {
    prevBonded[i] = onPrev[i].a1;
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

void DCLinkedHedron::PrepareNew(TrialMol &newMol, uint molIndex) {
  // Get new bond information
  SetBondLengthNew(newMol);
  hed.SetBondNew(bondLength, anchorBond);
  hed.PrepareNew(newMol, molIndex);
  bondEnergy = 0.0;
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondEnergy += data->ff.bonds.Calc(bondKinds[i], bondLength[i]);
  }
}

void DCLinkedHedron::PrepareOld(TrialMol &oldMol, uint molIndex) {
  // Get old bond information
  SetBondLengthOld(oldMol);
  hed.SetBondOld(bondLengthOld, anchorBondOld);
  hed.PrepareOld(oldMol, molIndex);
  bondEnergy = 0.0;
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondEnergy += data->ff.bonds.Calc(bondKinds[i], bondLengthOld[i]);
  }
}

void DCLinkedHedron::SetBondLengthNew(TrialMol &newMol) {
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondLength[i] = data->ff.bonds.Length(bondKinds[i]);
  }
  // anchorBond is built, we need the actual length
  anchorBond = sqrt(newMol.OldDistSq(hed.Focus(), hed.Prev()));
}

void DCLinkedHedron::SetBondLengthOld(TrialMol &oldMol) {
  for (uint i = 0; i < hed.NumBond(); ++i) {
    bondLengthOld[i] = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Bonded(i)));
  }
  anchorBondOld = sqrt(oldMol.OldDistSq(hed.Focus(), hed.Prev()));
}

void DCLinkedHedron::BuildNew(TrialMol &newMol, uint molIndex) {
  PRNG &prng = data->prng;
  // const CalculateEnergy& calc = data->calc;
  // const Forcefield& ff = data->ff;
  uint nLJTrials = data->nLJTrialsNth;
  uint nDihTrials = data->nDihTrials;
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
    newMol.AddAtom(hed.Bonded(b), positions[b][winner]);
    newMol.AddBonds(hed.Bonded(b), hed.Focus());
  }
  newMol.UpdateOverlap(overlap[winner]);
  newMol.AddEnergy(
      Energy(bondedEn[winner] + hed.GetEnergy() + bondEnergy,
             nonbonded[winner] + hed.GetNonBondedEn() + oneFour[winner],
             inter[winner], real[winner], 0.0, 0.0, 0.0));
  newMol.MultWeight(hed.GetWeight());
  newMol.MultWeight(stepWeight / nLJTrials);
}

void DCLinkedHedron::BuildOld(TrialMol &oldMol, uint molIndex) {
  PRNG &prng = data->prng;
  // const CalculateEnergy& calc = data->calc;
  const Forcefield &ff = data->ff;
  uint nLJTrials = data->nLJTrialsNth;
  uint nDihTrials = data->nDihTrials;
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
  ljWeights[0] = 0.0;
  for (uint tor = 0; tor < nDihTrials; ++tor) {
    torsion[tor] = (tor == 0) ? 0.0 : data->prng.rand(2.0 * M_PI);
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
  oldMol.UpdateOverlap(overlap[0]);
  oldMol.AddEnergy(Energy(bondedEn[0] + hed.GetEnergy() + bondEnergy,
                          nonbonded[0] + hed.GetNonBondedEn() + oneFour[0],
                          inter[0], real[0], 0.0, 0.0, 0.0));

  oldMol.MultWeight(hed.GetWeight());
  oldMol.MultWeight(stepWeight / nLJTrials);
}

double DCLinkedHedron::EvalLJ(TrialMol &mol, uint molIndex) {
  uint nLJTrials = data->nLJTrialsNth;
  double *inter = data->inter;
  double *nonbonded = data->nonbonded;
  double *real = data->real;
  bool *overlap = data->overlap;
  XYZArray *positions = data->multiPositions;

  std::fill_n(data->inter, nLJTrials, 0.0);
  std::fill_n(data->nonbonded, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);

  for (uint b = 0; b < hed.NumBond(); ++b) {
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

void DCLinkedHedron::ChooseTorsion(TrialMol &mol, uint molIndex,
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
    torsion[tor] = data->prng.rand(2.0 * M_PI);
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

} // namespace cbmc
