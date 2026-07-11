/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "DCOnSphere.h"

#include "DCData.h"
#include "Forcefield.h"
#include "MolSetup.h"
#include "PRNG.h"
#include "TrialMol.h"
#include "XYZArray.h"

namespace cbmc {
DCOnSphere::DCOnSphere(DCData *data, const mol_setup::MolKind kind, uint atom,
                       uint focus)
    : data(data), atom(atom), focus(focus) {
  using namespace mol_setup;
  std::vector<Bond> bonds = AtomBonds(kind, atom);
  for (uint i = 0; i < bonds.size(); ++i) {
    if (bonds[i].a0 == focus || bonds[i].a1 == focus) {
      bondKind = bonds[i].kind;
      break;
    }
  }

  if (data->nLJTrialsNth < 1) {
    std::cout << "Error: CBMC secondary atom trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

void DCOnSphere::SetBondLengthNew(TrialMol &newMol) {
  bondLength = data->ff.bonds.Length(bondKind);
}

void DCOnSphere::SetBondLengthOld(TrialMol &oldMol) {
  bondLengthOld = sqrt(oldMol.OldDistSq(focus, atom));
}

double DCOnSphere::BondEnergyNew(TrialMol &newMol) {
  return data->ff.bonds.Calc(bondKind, bondLength);
}

double DCOnSphere::BondEnergyOld(TrialMol &oldMol) {
  return data->ff.bonds.Calc(bondKind, bondLengthOld);
}

void DCOnSphere::BuildOld(TrialMol &oldMol, uint molIndex) {
  XYZArray &positions = data->positions;
  uint nLJTrials = data->nLJTrialsNth;
  double *inter = data->inter;
  double *real = data->real;
  bool *overlap = data->overlap;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);
  // calculate bond energy for old molecule.
  SetBondLengthOld(oldMol);
  bondEnergy = BondEnergyOld(oldMol);

  data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
                                    oldMol.AtomPosition(focus));
  positions.Set(0, oldMol.AtomPosition(atom));
  data->axes.WrapPBC(positions, oldMol.GetBox());

  data->calc.ParticleInter(inter, real, positions, overlap, atom, molIndex,
                           oldMol.GetBox(), nLJTrials);

  for (uint trial = 0; trial < nLJTrials; trial++) {
    stepWeight += exp(-1 * data->ff.beta * (inter[trial] + real[trial]));
  }
  oldMol.UpdateOverlap(overlap[0]);
  oldMol.MultWeight(stepWeight / nLJTrials);
  oldMol.AddEnergy(Energy(bondEnergy, 0.0, inter[0], real[0], 0.0, 0.0, 0.0));
  oldMol.ConfirmOldAtom(atom);
  oldMol.AddBonds(atom, focus);
}

void DCOnSphere::BuildNew(TrialMol &newMol, uint molIndex) {
  XYZArray &positions = data->positions;
  uint nLJTrials = data->nLJTrialsNth;
  double *inter = data->inter;
  double *real = data->real;
  double *ljWeights = data->ljWeights;
  bool *overlap = data->overlap;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);
  // calculate bond energy for old molecule.
  SetBondLengthNew(newMol);
  bondEnergy = BondEnergyNew(newMol);

  data->prng.FillWithRandomOnSphere(positions, nLJTrials, bondLength,
                                    newMol.AtomPosition(focus));
  data->axes.WrapPBC(positions, newMol.GetBox());

  data->calc.ParticleInter(inter, real, positions, overlap, atom, molIndex,
                           newMol.GetBox(), nLJTrials);

  for (uint trial = 0; trial < nLJTrials; trial++) {
    ljWeights[trial] = exp(-1 * data->ff.beta * (inter[trial] + real[trial]));
    stepWeight += ljWeights[trial];
  }
  uint winner = data->prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  newMol.UpdateOverlap(overlap[winner]);
  newMol.MultWeight(stepWeight / nLJTrials);
  newMol.AddEnergy(
      Energy(bondEnergy, 0, inter[winner], real[winner], 0.0, 0.0, 0.0));
  newMol.AddAtom(atom, positions[winner]);
  newMol.AddBonds(atom, focus);
}
} // namespace cbmc
