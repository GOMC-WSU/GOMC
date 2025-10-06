/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "DCRotateOnAtom.h"

#include <cassert>
#include <numeric>

#include "CalculateEnergy.h"
#include "Forcefield.h"
#include "Geometry.h"
#include "NumLib.h"
#include "PRNG.h"
#include "TrialMol.h"
#include "XYZArray.h"

namespace cbmc {
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

DCRotateOnAtom::DCRotateOnAtom(DCData *data, const mol_setup::MolKind &kind,
                               uint a0, uint a1, uint a2)
    : data(data), a0(a0), a1(a1), a2(a2) {
  using namespace mol_setup;
  std::vector<bool> visited(kind.atoms.size(), false);
  totAtoms = kind.atoms.size();
  // Find all the atoms that bonds with atoms a1
  std::vector<Bond> bonds = AtomBonds(kind, a1);
  // Remove the a0-a1 and a1-a2 bond
  bonds.erase(remove_if(bonds.begin(), bonds.end(), FindA1(a0)), bonds.end());
  bonds.erase(remove_if(bonds.begin(), bonds.end(), FindA1(a2)), bonds.end());
  visited[a1] = true;

  // Loop through other atoms that are bonded to a1
  for (uint b = 0; b < bonds.size(); b++) {
    // Store the atom index if it doesnot exist
    if (!visited[bonds[b].a0]) {
      atoms.push_back(bonds[b].a0);
      visited[bonds[b].a0] = true;
    }

    std::vector<Bond> temp = AtomBonds(kind, bonds[b].a1);
    for (uint i = 0; i < temp.size(); i++) {
      if (!visited[temp[i].a0]) {
        bonds.push_back(temp[i]);
      }
    }
  }

  numAtom = atoms.size();

  multiPosRotions = new XYZArray[numAtom];
  for (uint i = 0; i < numAtom; ++i) {
    multiPosRotions[i] = XYZArray(data->nLJTrialsNth);
  }

  if (data->nLJTrialsNth < 1) {
    std::cout << "Error: CBMC secondary atom trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }

  if (data->nDihTrials < 1) {
    std::cout << "Error: CBMC dihedral trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }

  // Find the angles affected by rotation
  // First find the angle x-a1-x
  ang = AtomMidAngles(kind, a1);
  // Remove the a0-a1-a2 angle
  ang.erase(remove_if(ang.begin(), ang.end(), FindAngle(a0, a2)), ang.end());

  // Find the dihedral affected by rotation
  // First find the dihedral with x-a0-a1-x in the middle
  std::vector<Dihedral> tempDih = DihsOnBond(kind, a1, a0);
  for (uint i = 0; i < tempDih.size(); i++) {
    // Make sure that the dihedral atoms are in the list since they are
    // constant.
    if (std::find(atoms.begin(), atoms.end(), tempDih[i].a0) != atoms.end()) {
      dih.push_back(tempDih[i]);
    }
  }
  // Add dihedral with x-a1-a2-x
  tempDih = DihsOnBond(kind, a1, a2);
  for (uint i = 0; i < tempDih.size(); i++) {
    // Make sure that the dihedral atoms are in the list since they are
    // constant.
    if (std::find(atoms.begin(), atoms.end(), tempDih[i].a0) != atoms.end()) {
      dih.push_back(tempDih[i]);
    }
  }

  // Add dihedral with atom a0 in one end: x-x-x-a0
  tempDih = AtomEndDihs(kind, a0);
  for (uint i = 0; i < tempDih.size(); i++) {
    // Make sure that the dihedral atoms are in the list since they are
    // constant.
    if (std::find(atoms.begin(), atoms.end(), tempDih[i].a3) != atoms.end()) {
      dih.push_back(tempDih[i]);
    }
  }
  // Add dihedral with atom a2 in one end: x-x-x-a2
  tempDih = AtomEndDihs(kind, a2);
  for (uint i = 0; i < tempDih.size(); i++) {
    // Make sure that the dihedral atoms are in the list since they are
    // constant.
    if (std::find(atoms.begin(), atoms.end(), tempDih[i].a3) != atoms.end()) {
      dih.push_back(tempDih[i]);
    }
  }
  /*
  for(uint i = 0; i < ang.size(); i++) {
    printf("R:Angle on %d-%d-%d: %d -> %d -> %d \n", a0, a1, a2, ang[i].a0,
  ang[i].a1, ang[i].a2);
  }
  for(uint i = 0; i < dih.size(); i++) {
    printf("R:Angle on %d-%d-%d: %d -> %d -> %d -> %d \n", a0, a1, a2,
  dih[i].a0, dih[i].a1, dih[i].a2, dih[i].a3);
  } */
}

void DCRotateOnAtom::PrepareOld(TrialMol &oldMol, uint molIndex) {
  for (uint a = 0; a < totAtoms; a++) {
    oldMol.ConfirmOldAtom(a);
  }

  XYZ center = oldMol.AtomPosition(a1);
  for (uint i = 0; i < numAtom; i++) {
    // Unwrap the coordinates with respect to a0.
    XYZ temp = oldMol.AtomPosition(atoms[i]);
    XYZ coord = data->axes.UnwrapPBC(temp, oldMol.GetBox(), center);
    // Shift the atoms to origin
    coord -= center;
    multiPosRotions[i].Set(0, coord);
  }
}

void DCRotateOnAtom::PrepareNew(TrialMol &newMol, uint molIndex) {
  for (uint a = 0; a < totAtoms; a++) {
    newMol.ConfirmOldAtom(a);
  }

  XYZ center = newMol.AtomPosition(a1);
  for (uint i = 0; i < numAtom; i++) {
    // Unwrap the coordinates with respect to a0.
    XYZ temp = newMol.AtomPosition(atoms[i]);
    XYZ coord = data->axes.UnwrapPBC(temp, newMol.GetBox(), center);
    // Shift the atoms to origin
    coord -= center;
    multiPosRotions[i].Set(0, coord);
  }
}

void DCRotateOnAtom::BuildOld(TrialMol &oldMol, uint molIndex) {
  PRNG &prng = data->prng;
  uint nLJTrials = data->nLJTrialsNth;
  uint nDihTrials = data->nDihTrials;
  double *torsion = data->angles;
  double *torWeights = data->angleWeights;
  double *torEnergy = data->angleEnergy;
  double *bondedEn = data->bonded;
  double *nonbonded = data->nonbonded;
  double *ljWeights = data->ljWeights;
  double *inter = data->inter;
  double *real = data->real;
  bool *overlap = data->overlap;
  double stepWeight = 0.0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(nonbonded, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  // Set up rotation matrix using.
  // rotaxis = (a0 - a2)
  XYZ center = oldMol.AtomPosition(a1);
  XYZ rotAxis = oldMol.AtomPosition(a0) - oldMol.AtomPosition(a2);
  rotAxis = data->axes.MinImage(rotAxis, oldMol.GetBox());
  rotAxis.Normalize();
  RotationMatrix cross = RotationMatrix::CrossProduct(rotAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotAxis);

  // Spin all nLJTrial except the original coordinates
  for (uint lj = nLJTrials; lj-- > 1;) {
    ChooseTorsion(oldMol, molIndex, cross, tensor);
    ljWeights[lj] = std::accumulate(torWeights, torWeights + nDihTrials, 0.0);
    uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
    bondedEn[lj] = torEnergy[winner];
    // convert chosen torsion to 3D positions
    RotationMatrix spin =
        RotationMatrix::FromAxisAngle(torsion[winner], cross, tensor);

    for (uint a = 0; a < numAtom; a++) {
      // find positions
      multiPosRotions[a].Set(lj, spin.Apply(multiPosRotions[a][0]));
      multiPosRotions[a].Add(lj, center);
    }
  }

  ChooseTorsionOld(oldMol, molIndex, cross, tensor);
  ljWeights[0] = std::accumulate(torWeights, torWeights + nDihTrials, 0.0);
  bondedEn[0] = torEnergy[0];

  for (uint a = 0; a < numAtom; a++) {
    // Shift original coordinate back.
    multiPosRotions[a].Add(0, center);
    // Wrap the atom coordinates
    data->axes.WrapPBC(multiPosRotions[a], oldMol.GetBox());
    // Calculate nonbonded energy
    data->calc.ParticleInter(inter, real, multiPosRotions[a], overlap, atoms[a],
                             molIndex, oldMol.GetBox(), nLJTrials);
    ParticleNonbonded(oldMol, multiPosRotions[a], atoms[a], nLJTrials);
  }

  for (uint trial = 0; trial < nLJTrials; ++trial) {
    // Skip exp() calculation for small values for efficiency and to avoid
    // subnormal values that contribute to differences between processors.
    // Value chosen mathematically: See cppreference.com exp function notes.
    double betaWeight =
        -data->ff.beta * (inter[trial] + real[trial] + nonbonded[trial]);
    if (betaWeight >= num::MIN_EXP_NONZERO_VAL) {
      ljWeights[trial] *= std::exp(betaWeight);
      stepWeight += ljWeights[trial];
    } else {
      ljWeights[trial] = 0.0;
    }
  }
  oldMol.UpdateOverlap(overlap[0]);
  oldMol.MultWeight(stepWeight / nLJTrials);
  oldMol.AddEnergy(
      Energy(bondedEn[0], nonbonded[0], inter[0], real[0], 0.0, 0.0, 0.0));

  for (uint a = 0; a < numAtom; a++) {
    oldMol.AddAtom(atoms[a], multiPosRotions[a][0]);
  }
}

void DCRotateOnAtom::BuildNew(TrialMol &newMol, uint molIndex) {
  PRNG &prng = data->prng;
  uint nLJTrials = data->nLJTrialsNth;
  uint nDihTrials = data->nDihTrials;
  double *torsion = data->angles;
  double *torWeights = data->angleWeights;
  double *torEnergy = data->angleEnergy;
  double *bondedEn = data->bonded;
  double *nonbonded = data->nonbonded;
  double *ljWeights = data->ljWeights;
  double *inter = data->inter;
  double *real = data->real;
  bool *overlap = data->overlap;
  double stepWeight = 0.0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(nonbonded, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  // Set up rotation matrix
  // rotaxis = (a0 - a2)
  XYZ center = newMol.AtomPosition(a1);
  XYZ rotAxis = newMol.AtomPosition(a0) - newMol.AtomPosition(a2);
  rotAxis = data->axes.MinImage(rotAxis, newMol.GetBox());
  rotAxis.Normalize();
  RotationMatrix cross = RotationMatrix::CrossProduct(rotAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotAxis);

  // Go backward to to preserve prototype
  for (uint lj = nLJTrials; lj-- > 0;) {
    ChooseTorsion(newMol, molIndex, cross, tensor);
    ljWeights[lj] = std::accumulate(torWeights, torWeights + nDihTrials, 0.0);
    uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
    bondedEn[lj] = torEnergy[winner];
    // convert chosen torsion to 3D positions
    RotationMatrix spin =
        RotationMatrix::FromAxisAngle(torsion[winner], cross, tensor);

    for (uint a = 0; a < numAtom; a++) {
      // find positions
      multiPosRotions[a].Set(lj, spin.Apply(multiPosRotions[a][0]));
      multiPosRotions[a].Add(lj, center);
    }
  }

  for (uint a = 0; a < numAtom; a++) {
    // Wrap the atom coordinates
    data->axes.WrapPBC(multiPosRotions[a], newMol.GetBox());
    // Calculate nonbonded energy
    data->calc.ParticleInter(inter, real, multiPosRotions[a], overlap, atoms[a],
                             molIndex, newMol.GetBox(), nLJTrials);
    ParticleNonbonded(newMol, multiPosRotions[a], atoms[a], nLJTrials);
  }

  for (uint trial = 0; trial < nLJTrials; trial++) {
    // Skip exp() calculation for small values for efficiency and to avoid
    // subnormal values that contribute to differences between processors.
    // Value chosen mathematically: See cppreference.com exp function notes.
    double betaWeight =
        -data->ff.beta * (inter[trial] + real[trial] + nonbonded[trial]);
    if (betaWeight >= num::MIN_EXP_NONZERO_VAL) {
      ljWeights[trial] *= std::exp(betaWeight);
      stepWeight += ljWeights[trial];
    } else {
      ljWeights[trial] = 0.0;
    }
  }
  uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  newMol.UpdateOverlap(overlap[winner]);
  newMol.MultWeight(stepWeight / nLJTrials);
  newMol.AddEnergy(Energy(bondedEn[winner], nonbonded[winner], inter[winner],
                          real[winner], 0.0, 0.0, 0.0));

  for (uint a = 0; a < numAtom; a++) {
    newMol.AddAtom(atoms[a], multiPosRotions[a][winner]);
  }
}

void DCRotateOnAtom::ChooseTorsion(TrialMol &mol, uint molIndex,
                                   RotationMatrix &cross,
                                   RotationMatrix &tensor) {
  uint nDihTrials = data->nDihTrials;
  double *torsion = data->angles;
  double *torWeights = data->angleWeights;
  double *torEnergy = data->angleEnergy;

  XYZ center = mol.AtomPosition(a1);
  for (uint tor = 0; tor < nDihTrials; ++tor) {
    torsion[tor] = data->prng.rand(2.0 * M_PI);
    // convert chosen torsion to 3D positions
    RotationMatrix spin =
        RotationMatrix::FromAxisAngle(torsion[tor], cross, tensor);
    for (uint a = 0; a < numAtom; a++) {
      XYZ coord = spin.Apply(multiPosRotions[a][0]);
      mol.AddAtom(atoms[a], coord + center);
    }
    double en = CalcIntraBonded(mol, molIndex);
    torEnergy[tor] = en;
    torWeights[tor] = std::exp(-data->ff.beta * en);
  }
}

void DCRotateOnAtom::ChooseTorsionOld(TrialMol &mol, uint molIndex,
                                      RotationMatrix &cross,
                                      RotationMatrix &tensor) {
  uint nDihTrials = data->nDihTrials;
  double *torsion = data->angles;
  double *torWeights = data->angleWeights;
  double *torEnergy = data->angleEnergy;

  XYZ center = mol.AtomPosition(a1);
  for (uint tor = 0; tor < nDihTrials; ++tor) {
    // Use actual coordinate fir first torsion trial
    torsion[tor] = (tor == 0) ? 0.0 : data->prng.rand(2.0 * M_PI);
    // convert chosen torsion to 3D positions
    RotationMatrix spin =
        RotationMatrix::FromAxisAngle(torsion[tor], cross, tensor);
    for (uint a = 0; a < numAtom; a++) {
      XYZ coord = spin.Apply(multiPosRotions[a][0]);
      mol.AddAtom(atoms[a], coord + center);
    }
    double en = CalcIntraBonded(mol, molIndex);
    torEnergy[tor] = en;
    torWeights[tor] = std::exp(-data->ff.beta * en);
  }
}

double DCRotateOnAtom::CalcIntraBonded(TrialMol &mol, uint molIndex) {
  double bondedEn = 0.0;
  uint box = mol.GetBox();
  XYZ b1, b2, b3;
  const XYZArray &coords = mol.GetCoords();
  for (uint i = 0; i < ang.size(); i++) {
    b1 = data->axes.MinImage(coords.Difference(ang[i].a0, ang[i].a1), box);
    b2 = data->axes.MinImage(coords.Difference(ang[i].a2, ang[i].a1), box);
    bondedEn += data->ff.angles->Calc(ang[i].kind, geom::Theta(b1, b2));
  }

  for (uint i = 0; i < dih.size(); i++) {
    // No need to calcminImage since we unwrap it
    b1 = data->axes.MinImage(coords.Difference(dih[i].a1, dih[i].a0), box);
    b2 = data->axes.MinImage(coords.Difference(dih[i].a2, dih[i].a1), box);
    b3 = data->axes.MinImage(coords.Difference(dih[i].a3, dih[i].a2), box);
    bondedEn += data->ff.dihedrals.Calc(dih[i].kind, geom::Phi(b1, b2, b3));
  }
  return bondedEn;
}

void DCRotateOnAtom::ParticleNonbonded(cbmc::TrialMol const &mol,
                                       XYZArray const &trialPos,
                                       const uint partIndex,
                                       const uint trials) {
  ParticleNonbonded1_N(mol, trialPos, partIndex, trials);
  ParticleNonbonded1_4(mol, trialPos, partIndex, trials);
  ParticleNonbonded1_3(mol, trialPos, partIndex, trials);
}

void DCRotateOnAtom::ParticleNonbonded1_N(cbmc::TrialMol const &mol,
                                          XYZArray const &trialPos,
                                          const uint partIndex,
                                          const uint trials) {
  double *nonbonded = data->nonbonded;
  uint box = mol.GetBox();
  const MoleculeKind &kind = mol.GetKind();
  // loop over all partners of the trial particle
  const uint *partner = kind.sortedNB.Begin(partIndex);
  const uint *end = kind.sortedNB.End(partIndex);
  while (partner != end) {
    if (mol.AtomExists(*partner) &&
        (std::find(atoms.begin(), atoms.end(), *partner) == atoms.end())) {
      for (uint t = 0; t < trials; ++t) {
        double distSq;
        if (data->axes.InRcut(distSq, trialPos, t, mol.GetCoords(), *partner,
                              box)) {
          nonbonded[t] += data->ff.particles->CalcEn(
              distSq, kind.AtomKind(partIndex), kind.AtomKind(*partner), 1.0);
          if (data->ff.electrostatic) {
            double qi_qj_Fact = kind.AtomCharge(partIndex) *
                                kind.AtomCharge(*partner) * num::qqFact;
            data->ff.particles->CalcCoulombAdd_1_4(nonbonded[t], distSq,
                                                   qi_qj_Fact, true);
          }
        }
      }
    }
    ++partner;
  }
}

void DCRotateOnAtom::ParticleNonbonded1_4(cbmc::TrialMol const &mol,
                                          XYZArray const &trialPos,
                                          const uint partIndex,
                                          const uint trials) {
  if (!data->ff.OneFour)
    return;

  double *nonbonded = data->nonbonded;
  uint box = mol.GetBox();
  const MoleculeKind &kind = mol.GetKind();
  // loop over all partners of the trial particle
  const uint *partner = kind.sortedNB_1_4.Begin(partIndex);
  const uint *end = kind.sortedNB_1_4.End(partIndex);
  while (partner != end) {
    if (mol.AtomExists(*partner) &&
        (std::find(atoms.begin(), atoms.end(), *partner) == atoms.end())) {
      for (uint t = 0; t < trials; ++t) {
        double distSq;
        if (data->axes.InRcut(distSq, trialPos, t, mol.GetCoords(), *partner,
                              box)) {
          data->ff.particles->CalcAdd_1_4(nonbonded[t], distSq,
                                          kind.AtomKind(partIndex),
                                          kind.AtomKind(*partner));
          if (data->ff.electrostatic) {
            double qi_qj_Fact = kind.AtomCharge(partIndex) *
                                kind.AtomCharge(*partner) * num::qqFact;
            data->ff.particles->CalcCoulombAdd_1_4(nonbonded[t], distSq,
                                                   qi_qj_Fact, false);
          }
        }
      }
    }
    ++partner;
  }
}

void DCRotateOnAtom::ParticleNonbonded1_3(cbmc::TrialMol const &mol,
                                          XYZArray const &trialPos,
                                          const uint partIndex,
                                          const uint trials) {
  if (!data->ff.OneThree)
    return;

  double *nonbonded = data->nonbonded;
  uint box = mol.GetBox();
  const MoleculeKind &kind = mol.GetKind();
  // loop over all partners of the trial particle
  const uint *partner = kind.sortedNB_1_3.Begin(partIndex);
  const uint *end = kind.sortedNB_1_3.End(partIndex);
  while (partner != end) {
    if (mol.AtomExists(*partner) &&
        (std::find(atoms.begin(), atoms.end(), *partner) == atoms.end())) {
      for (uint t = 0; t < trials; ++t) {
        double distSq;
        if (data->axes.InRcut(distSq, trialPos, t, mol.GetCoords(), *partner,
                              box)) {
          data->ff.particles->CalcAdd_1_4(nonbonded[t], distSq,
                                          kind.AtomKind(partIndex),
                                          kind.AtomKind(*partner));
          if (data->ff.electrostatic) {
            double qi_qj_Fact = kind.AtomCharge(partIndex) *
                                kind.AtomCharge(*partner) * num::qqFact;
            data->ff.particles->CalcCoulombAdd_1_4(nonbonded[t], distSq,
                                                   qi_qj_Fact, false);
          }
        }
      }
    }
    ++partner;
  }
}

} // namespace cbmc
