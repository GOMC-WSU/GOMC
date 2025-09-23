/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "DCRotateCOM.h"
#include "DCData.h"
#include "Forcefield.h"
#include "MolSetup.h"
#include "NumLib.h"
#include "PRNG.h"
#include "TrialMol.h"
/*
This file is only called from MoleculeExchange and IntraMoleculeExchange file.
This file depend on ensemble kind, swap the COM of small molecules with one
large molecule.

Molecule is allowed to have trial position in Box or trial position
in cavity (HasCav()), or no trial position (COMFix()).

Molecule is allowed to have rotation around COM or around Z-axis
*/

namespace cbmc {

DCRotateCOM::DCRotateCOM(DCData *data, const mol_setup::MolKind kind)
    : data(data), rotateMatrix(3), invMatrix(3) {
  rotateMatrix.Set(0, 0.0, 0.0, 0.0);
  rotateMatrix.Set(1, 0.0, 0.0, 0.0);
  rotateMatrix.Set(2, 0.0, 0.0, 1.0);
  atomNumber = kind.atoms.size();
  multiPosRotions = new XYZArray[atomNumber];

  if (data->totalTrials < 1) {
    std::cout << "Error: CBMC first or secondary atom trials must be greater "
                 "than 0.\n";
    exit(EXIT_FAILURE);
  }
  for (uint i = 0; i < atomNumber; ++i) {
    multiPosRotions[i] = XYZArray(data->totalTrials);
  }
}

void DCRotateCOM::RandRotateZ() {
  PRNG &prng = data->prng;
  double theta = prng();
  theta *= 2.0 * M_PI;
  theta -= M_PI;
  double cosTheta = cos(theta);
  double sinTheta = sin(theta);
  rotateMatrix.Set(0, cosTheta, -sinTheta, 0.0);
  rotateMatrix.Set(1, sinTheta, cosTheta, 0.0);
}

void DCRotateCOM::PrepareNew(TrialMol &newMol, uint molIndex) {
  newMol.SetWeight(1.0);
  // old center of mass
  oldCOM = newMol.GetCOM();
}

void DCRotateCOM::PickTransferCOMNew(TrialMol &newMol, uint molIndex) {
  PRNG &prng = data->prng;

  if (newMol.COMFix()) {
    COM = newMol.GetCavityCenter();
  } else if (newMol.HasCav()) {
    // new center of mass that need to be transfered
    // Pick molecule in cav dimension
    prng.FillWithRandomInCavity(COM, newMol.GetCavity());
    // rotate using cavity matrix
    COM = newMol.Transform(COM);
    // add center
    COM += newMol.GetCavityCenter();
  } else {
    prng.FillWithRandom(COM, data->axes, newMol.GetBox());
  }

  XYZ diff = COM - oldCOM;

  for (uint p = 0; p < atomNumber; p++) {
    newMol.SetAtomCoords(p, newMol.AtomPosition(p) + diff);
  }

  oldCOM = COM;
}

void DCRotateCOM::PrepareOld(TrialMol &oldMol, uint molIndex) {
  oldMol.SetWeight(1.0);
  // old center of mass
  oldCOM = oldMol.GetCOM();
  COM = oldCOM;
}

void DCRotateCOM::PickTransferCOMOld(TrialMol &oldMol, uint molIndex) {
  PRNG &prng = data->prng;
  // new center of mass that need to be transfered
  if (oldMol.HasCav()) {
    // Pick molecule in cav dimension
    prng.FillWithRandomInCavity(COM, oldMol.GetCavity());
    // rotate using cavity matrix
    COM = oldMol.Transform(COM);
    // add center
    COM += oldMol.GetCavityCenter();
  } else
    prng.FillWithRandom(COM, data->axes, oldMol.GetBox());

  XYZ diff = COM - oldCOM;

  for (uint p = 0; p < atomNumber; p++) {
    oldMol.SetAtomCoords(p, oldMol.AtomPosition(p) + diff);
  }

  oldCOM = COM;
}

void DCRotateCOM::BuildNew(TrialMol &newMol, uint molIndex) {
  PRNG &prng = data->prng;
  const CalculateEnergy &calc = data->calc;
  // const Ewald *calcEwald = data->calcEwald;
  const Forcefield &ff = data->ff;
  uint nLJTrials = data->nLJTrialsNth;
  uint fLJTrials = data->nLJTrialsFirst;
  uint totalTrials = data->totalTrials;
  double *ljWeights = data->ljWeightsT;
  double *inter = data->interT;
  double *real = data->realT;
  bool *overlap = data->overlapT;
  RotationMatrix spin;

  std::fill_n(inter, totalTrials, 0.0);
  std::fill_n(real, totalTrials, 0.0);
  std::fill_n(ljWeights, totalTrials, 0.0);
  std::fill_n(overlap, totalTrials, false);

  if (atomNumber == 1) {
    nLJTrials = 1;
    totalTrials = fLJTrials;
  }

  if (newMol.COMFix()) {
    fLJTrials = 1;
    totalTrials = nLJTrials;
    // if we rotate around backbone we need to calc the rotation matrix
    if (newMol.RotateBB()) {
      // find the inverse matrix of molecule that we insert
      XYZ backBone;
      if (atomNumber != 1) {
        backBone = newMol.GetCoords().Difference(newMol.GetAtomBB(0),
                                                 newMol.GetAtomBB(1));
      } else {
        backBone = prng.RandomUnitVect();
      }
      XYZArray T(3);
      geom::SetBasis(T, backBone);
      geom::TransposeMatrix(invMatrix, T);
    }
  }

  for (uint p = 0; p < fLJTrials; ++p) {
    // Pick a new position for COM and transfer the molecule
    PickTransferCOMNew(newMol, molIndex);
    // get info about existing geometry
    newMol.ShiftBasis(COM);
    const XYZ center = COM;
    uint index = p * nLJTrials;

    for (uint a = 0; a < atomNumber; ++a) {
      multiPosRotions[a].Set(index, newMol.AtomPosition(a));
      multiPosRotions[a].Add(index, -center);
    }

    // Rotational trial the molecule around COM
    for (uint r = nLJTrials; r-- > 0;) {
      if (newMol.RotateBB()) {
        // we only perform rotation around z axis
        RandRotateZ();
      } else {
        // convert chosen torsion to 3D positions
        double u1, u2, u3;
        u1 = prng();
        u2 = prng();
        u3 = prng();
        spin = RotationMatrix::UniformRandom(u1, u2, u3);
      }

      for (uint a = 0; a < atomNumber; ++a) {
        if (newMol.RotateBB()) {
          XYZ coord = multiPosRotions[a][index];
          // transform backbone to z axis
          coord = geom::Transform(invMatrix, coord);
          // rotate around z
          coord = geom::Transform(rotateMatrix, coord);
          // transfer backbone to cavity orientation
          coord = newMol.Transform(coord);
          multiPosRotions[a].Set(index + r, coord);
        } else {
          // find positions
          multiPosRotions[a].Set(index + r,
                                 spin.Apply(multiPosRotions[a][index]));
        }
        multiPosRotions[a].Add(index + r, center);
      }
    }
  }

  for (uint a = 0; a < atomNumber; ++a) {
    data->axes.WrapPBC(multiPosRotions[a], newMol.GetBox());
    calc.ParticleInter(inter, real, multiPosRotions[a], overlap, a, molIndex,
                       newMol.GetBox(), totalTrials);
  }

  double stepWeight = 0.0;
  for (uint lj = 0; lj < totalTrials; ++lj) {
    // Skip exp() calculation for small values for efficiency and to avoid
    // subnormal values that contribute to differences between processors.
    // Value chosen mathematically: See cppreference.com exp function notes.
    // Note: ljWeights prefilled with 0.0, so don't need to initialize it.
    double betaWeight = -ff.beta * (inter[lj] + real[lj]);
    if (betaWeight >= num::MIN_EXP_NONZERO_VAL) {
      ljWeights[lj] = std::exp(betaWeight);
      stepWeight += ljWeights[lj];
    }
  }
  uint winner = prng.PickWeighted(ljWeights, totalTrials, stepWeight);
  for (uint a = 0; a < atomNumber; ++a) {
    newMol.AddAtom(a, multiPosRotions[a][winner]);
  }

  newMol.UpdateOverlap(overlap[winner]);
  newMol.AddEnergy(
      Energy(0.0, 0.0, inter[winner], real[winner], 0.0, 0.0, 0.0));
  newMol.MultWeight(stepWeight / totalTrials);
}

void DCRotateCOM::BuildOld(TrialMol &oldMol, uint molIndex) {
  PRNG &prng = data->prng;
  const CalculateEnergy &calc = data->calc;
  // const Ewald * calcEwald = data->calcEwald;
  const Forcefield &ff = data->ff;
  uint nLJTrials = data->nLJTrialsNth;
  uint fLJTrials = data->nLJTrialsFirst;
  uint totalTrials = data->totalTrials;
  double *ljWeights = data->ljWeightsT;
  double *inter = data->interT;
  double *real = data->realT;
  bool *overlap = data->overlapT;
  RotationMatrix spin;

  std::fill_n(inter, totalTrials, 0.0);
  std::fill_n(real, totalTrials, 0.0);
  std::fill_n(ljWeights, totalTrials, 0.0);
  std::fill_n(overlap, totalTrials, false);

  if (atomNumber == 1) {
    nLJTrials = 1;
    totalTrials = fLJTrials;
  }

  if (oldMol.COMFix()) {
    fLJTrials = 1;
    totalTrials = nLJTrials;
    // if we rotate around backbone of the molecule
    if (oldMol.RotateBB()) {
      // find the inverse matrix of cavity
      oldMol.TransposeMatrix(invMatrix);
    }
  }

  const XYZ orgCenter = COM;

  for (uint p = 0; p < fLJTrials; ++p) {
    // First trial is current configuration
    // get info about existing geometry
    oldMol.ShiftBasis(COM);
    const XYZ center = COM;
    uint index = p * nLJTrials;

    for (uint a = 0; a < atomNumber; ++a) {
      // get position and shift to origin
      multiPosRotions[a].Set(index, oldMol.AtomPosition(a));
      multiPosRotions[a].Add(index, -center);
    }

    // Rotational trial the molecule around COM
    for (uint r = nLJTrials; r-- > 0;) {
      if ((index + r) == 0)
        continue;

      if (oldMol.RotateBB()) {
        // we only perform rotation around z axis
        RandRotateZ();
      } else {
        // convert chosen torsion to 3D positions
        double u1, u2, u3;
        u1 = prng();
        u2 = prng();
        u3 = prng();
        spin = RotationMatrix::UniformRandom(u1, u2, u3);
      }

      for (uint a = 0; a < atomNumber; ++a) {
        if (oldMol.RotateBB()) {
          XYZ coord = multiPosRotions[a][index];
          // transform backbone to z axis
          coord = geom::Transform(invMatrix, coord);
          // rotate around z
          coord = geom::Transform(rotateMatrix, coord);
          // transfer backbone to cavity orientation
          coord = oldMol.Transform(coord);
          multiPosRotions[a].Set(index + r, coord);
        } else {
          // find positions
          multiPosRotions[a].Set(index + r,
                                 spin.Apply(multiPosRotions[a][index]));
        }
        multiPosRotions[a].Add(index + r, center);
      }
    }

    // Pick a new position for COM and transfer the molecule
    PickTransferCOMOld(oldMol, molIndex);
  }

  for (uint a = 0; a < atomNumber; ++a) {
    multiPosRotions[a].Add(0, orgCenter);
    data->axes.WrapPBC(multiPosRotions[a], oldMol.GetBox());
    calc.ParticleInter(inter, real, multiPosRotions[a], overlap, a, molIndex,
                       oldMol.GetBox(), totalTrials);
  }

  double stepWeight = 0.0;
  for (uint lj = 0; lj < totalTrials; ++lj) {
    // Skip exp() calculation for small values for efficiency and to avoid
    // subnormal values that contribute to differences between processors.
    // Value chosen mathematically: See cppreference.com exp function notes.
    double betaWeight = -ff.beta * (inter[lj] + real[lj]);
    if (betaWeight >= num::MIN_EXP_NONZERO_VAL) {
      stepWeight += std::exp(betaWeight);
    }
  }

  for (uint a = 0; a < atomNumber; ++a) {
    oldMol.AddAtom(a, multiPosRotions[a][0]);
  }

  oldMol.UpdateOverlap(overlap[0]);
  oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0], 0.0, 0.0, 0.0));
  oldMol.MultWeight(stepWeight / totalTrials);
}

} // namespace cbmc
