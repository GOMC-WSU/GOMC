/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.11
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCLink.h"
#include "TrialMol.h"
#include "Forcefield.h"
#include "XYZArray.h"
#include "MoleculeKind.h"
#include "MolSetup.h"
#include "NumLib.h"

namespace cbmc
{

DCLink::DCLink(DCData* data, const mol_setup::MolKind kind,
               uint atom, uint focus)

  : data(data), atom(atom), focus(focus), angleFix(false)

{
  //will fail quietly if not a part of a valid linear molecule,
  //but we checked that, right?
  using namespace mol_setup;
  std::vector<Bond> bonds = AtomBonds(kind, atom);
  for (uint i = 0; i < bonds.size(); ++i) {
    if (bonds[i].a0 == focus || bonds[i].a1 == focus) {
      bondLength = data->ff.bonds.Length(bonds[i].kind);
      bond[2] = bondLength;
      bondKind = bonds[i].kind;

      break;
    }
  }

  std::vector<Angle> angles = AtomEndAngles(kind, atom);
  for (uint i = 0; i < angles.size(); ++i) {
    if (angles[i].a1 == focus) {
      angleKind = angles[i].kind;

      if (data->ff.angles->AngleFixed(angleKind)) {
        angleFix = true;
        thetaFix = data->ff.angles->Angle(angleKind);
      }

      break;
    }
  }
  std::vector<Dihedral> dihs = AtomEndDihs(kind, atom);
  for (uint i = 0; i < dihs.size(); ++i) {
    if (dihs[i].a1 == focus) {
      dihKind = dihs[i].kind;
      prev = dihs[i].a2;
      prevprev = dihs[i].a3;
      break;
    }
  }

  std::vector<Bond> bond2 = AtomBonds(kind, focus);
  for (uint i = 0; i < bond2.size(); ++i) {
    if (bond2[i].a0 == prev || bond2[i].a1 == prev) {
      bond[1] = data->ff.bonds.Length(bond2[i].kind);
      break;
    }
  }

  std::vector<Bond> bond3 = AtomBonds(kind, prev);
  for (uint i = 0; i < bond3.size(); ++i) {
    if (bond3[i].a0 == prevprev || bond3[i].a1 == prevprev) {
      bond[0] = data->ff.bonds.Length(bond3[i].kind);
      break;
    }
  }
}

void DCLink::PrepareNew(TrialMol& newMol, uint molIndex)
{
  double* angles = data->angles;
  double* angleEnergy = data->angleEnergy;
  double* angleWeights = data->angleWeights;
  double* nonbonded_1_3 =  data->nonbonded_1_3;
  PRNG& prng = data->prng;
  const Forcefield& ff = data->ff;
  uint count = data->nAngleTrials;
  std::fill_n(nonbonded_1_3, count, 0.0);
  bendWeight = 0;

  for (uint trial = 0; trial < count; trial++) {
    if (angleFix) {
      angles[trial] = thetaFix;
      angleEnergy[trial] = 0.0;
    } else {
      angles[trial] = prng.rand(M_PI);
      angleEnergy[trial] = ff.angles->Calc(angleKind, angles[trial]);
    }


    double distSq = newMol.AngleDist(bond[1], bond[2], angles[trial]);
    nonbonded_1_3[trial] = data->calc.IntraEnergy_1_3(distSq, prev, atom,
                           molIndex);

    if(isnan(nonbonded_1_3[trial]))
      nonbonded_1_3[trial] = num::BIGNUM;

    angleWeights[trial] = exp((angleEnergy[trial] + nonbonded_1_3[trial])
                              * -ff.beta);

    bendWeight += angleWeights[trial];
  }
  uint winner = prng.PickWeighted(angleWeights, count, bendWeight);
  theta = angles[winner];
  bendEnergy = angleEnergy[winner];
  oneThree = nonbonded_1_3[winner];
}

void DCLink::PrepareOld(TrialMol& oldMol, uint molIndex)
{
  PRNG& prng = data->prng;
  const Forcefield& ff = data->ff;
  uint count = data->nAngleTrials - 1;
  bendWeight = 0;

  //set bond distance for old molecule
  double BondDistSq1 = oldMol.OldDistSq(focus, atom);
  double BondDistSq2 = oldMol.OldDistSq(prev, focus);
  double BondDistSq3 = oldMol.OldDistSq(prevprev, prev);
  SetOldMolBond(2, BondDistSq1);
  SetOldMolBond(1, BondDistSq2);
  SetOldMolBond(0, BondDistSq3);

  for (uint trial = 0; trial < count; trial++) {

    double trialAngle;
    double trialEn;

    if (angleFix) {
      trialAngle = thetaFix;
      trialEn = 0.0;
    } else {
      trialAngle = prng.rand(M_PI);
      trialEn = ff.angles->Calc(angleKind, trialAngle);
    }

    double distSq = oldMol.AngleDist(oldBond[1], oldBond[2], trialAngle);


    double tempEn = data->calc.IntraEnergy_1_3(distSq, prev, atom, molIndex);
    if(isnan(tempEn))
      tempEn = num::BIGNUM;

    trialEn += tempEn;


    double trialWeight = exp(-ff.beta * trialEn);
    bendWeight += trialWeight;
  }
}

void DCLink::IncorporateOld(TrialMol& oldMol, uint molIndex)
{
  oldMol.OldThetaAndPhi(atom, focus, theta, phi);
  const Forcefield& ff = data->ff;

  bendEnergy = ff.angles->Calc(angleKind, theta);
  double distSq = oldMol.OldDistSq(prev, atom);
  oneThree = data->calc.IntraEnergy_1_3(distSq, prev, atom, molIndex);
  bendWeight += exp(-ff.beta * (bendEnergy + oneThree));
  //considering bond energy for old molecule. There is no need to calculate
  //for new molecule since we dont sample bond.
  double BondDistSq = oldMol.OldDistSq(focus, atom);
  oldBondEnergy = ff.bonds.Calc(bondKind, sqrt(BondDistSq));
}

void DCLink::AlignBasis(TrialMol& mol)
{
  mol.SetBasis(focus, prev, prevprev);
}

void DCLink::SetOldMolBond(const uint i, const double distSq)
{
  oldBond[i] = sqrt(distSq);
}

void DCLink::BuildOld(TrialMol& oldMol, uint molIndex)
{
  AlignBasis(oldMol);
  IncorporateOld(oldMol, molIndex);
  double* angles = data->angles;
  double* angleEnergy = data->angleEnergy;
  double* angleWeights = data->angleWeights;
  double* ljWeights = data->ljWeights;
  double* torsion = data->bonded;
  double* nonbonded = data->nonbonded;
  double* nonbonded_1_4 = data->nonbonded_1_4;
  double* inter = data->inter;
  double* real = data->real;
  double* oneFour = data->oneFour;
  uint nLJTrials = data->nLJTrialsNth;
  XYZArray& positions = data->positions;
  PRNG& prng = data->prng;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(nonbonded, nLJTrials, 0.0);
  std::fill_n(nonbonded_1_4, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(angles, data->nDihTrials, 0.0);
  std::fill_n(angleWeights, data->nDihTrials, 0.0);
  std::fill_n(angleEnergy, data->nDihTrials, 0.0);
  std::fill_n(torsion, nLJTrials, 0.0);
  std::fill_n(oneFour, nLJTrials, 0.0);

  UseOldDih(oldMol, molIndex, torsion[0], ljWeights[0]);
  oneFour[0] = nonbonded_1_4[0];
  positions.Set(0, oldMol.AtomPosition(atom));

  for (uint trial = 1; trial < nLJTrials; ++trial) {
    ljWeights[trial] = GenerateDihedralsOld(oldMol, molIndex, angles,
                                            angleEnergy, angleWeights);
    uint winner = prng.PickWeighted(angleWeights, data->nDihTrials,
                                    ljWeights[trial]);
    torsion[trial] = angleEnergy[winner];
    oneFour[trial] = nonbonded_1_4[winner];
    positions.Set(trial, oldMol.GetRectCoords(bondLength, theta,
                  angles[winner]));
  }

  data->axes.WrapPBC(positions, oldMol.GetBox());
  data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                           oldMol.GetBox(), nLJTrials);

  data->calc.ParticleNonbonded(nonbonded, oldMol, positions, atom,
                               oldMol.GetBox(), nLJTrials);

  double dihLJWeight = 0;
  for (uint trial = 0; trial < nLJTrials; ++trial) {
    ljWeights[trial] *= exp(-data->ff.beta *
                            (inter[trial] + real[trial] +
                             nonbonded[trial]));
    dihLJWeight += ljWeights[trial];
  }
  oldMol.MultWeight(dihLJWeight * bendWeight);
  oldMol.ConfirmOldAtom(atom);
  oldMol.AddEnergy(Energy(torsion[0] + bendEnergy + oldBondEnergy,
                          nonbonded[0] + oneThree + oneFour[0],
                          inter[0], real[0], 0.0, 0.0, 0.0));
}

void DCLink::BuildNew(TrialMol& newMol, uint molIndex)
{
  AlignBasis(newMol);
  double* angles = data->angles;
  double* angleEnergy = data->angleEnergy;
  double* angleWeights = data->angleWeights;
  double* ljWeights = data->ljWeights;
  double* torsion = data->bonded;
  double* nonbonded = data->nonbonded;
  double* nonbonded_1_4 = data->nonbonded_1_4;
  double* inter = data->inter;
  double* real = data->real;
  double* oneFour = data->oneFour;
  uint nLJTrials = data->nLJTrialsNth;
  XYZArray& positions = data->positions;
  PRNG& prng = data->prng;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(nonbonded, nLJTrials, 0.0);
  std::fill_n(nonbonded_1_4, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(angles, data->nDihTrials, 0.0);
  std::fill_n(angleWeights, data->nDihTrials, 0.0);
  std::fill_n(angleEnergy, data->nDihTrials, 0.0);
  std::fill_n(torsion, nLJTrials, 0.0);
  std::fill_n(oneFour, nLJTrials, 0.0);

  for (uint trial = 0; trial < nLJTrials; ++trial) {
    ljWeights[trial] = GenerateDihedralsNew(newMol, molIndex, angles,
                                            angleEnergy, angleWeights);
    uint winner = prng.PickWeighted(angleWeights, data->nDihTrials,
                                    ljWeights[trial]);
    oneFour[trial] = nonbonded_1_4[winner];
    torsion[trial] = angleEnergy[winner];
    positions.Set(trial, newMol.GetRectCoords(bondLength, theta,
                  angles[winner]));
  }

  data->axes.WrapPBC(positions, newMol.GetBox());
  data->calc.ParticleInter(inter, real, positions, atom, molIndex,
                           newMol.GetBox(), nLJTrials);

  data->calc.ParticleNonbonded(nonbonded, newMol, positions, atom,
                               newMol.GetBox(), nLJTrials);

  double dihLJWeight = 0;
  double beta = data->ff.beta;
  for (uint trial = 0; trial < nLJTrials; ++trial) {
    ljWeights[trial] *= exp(-data->ff.beta *
                            (inter[trial] + real[trial] +
                             nonbonded[trial]));

    dihLJWeight += ljWeights[trial];
  }

  uint winner = prng.PickWeighted(ljWeights, nLJTrials, dihLJWeight);
  newMol.MultWeight(dihLJWeight * bendWeight);
  newMol.AddAtom(atom, positions[winner]);
  newMol.AddEnergy(Energy(torsion[winner] + bendEnergy, nonbonded[winner] +
                          oneThree + oneFour[winner], inter[winner],
                          real[winner], 0.0, 0.0, 0.0));
}

double DCLink::GenerateDihedralsNew(TrialMol& newMol, uint molIndex,
                                    double* angles, double* angleEnergy,
                                    double* angleWeights)
{
  double* nonbonded_1_4 = data->nonbonded_1_4;
  double stepWeight = 0.0;
  PRNG& prng = data->prng;
  const Forcefield& ff = data->ff;
  std::fill_n(nonbonded_1_4, data->nDihTrials, 0.0);

  double theta1 = newMol.GetTheta(prevprev, prev, focus);

  for (uint trial = 0, count = data->nDihTrials; trial < count; ++trial) {
    angles[trial] = prng.rand(2 * M_PI);
    angleEnergy[trial] = ff.dihedrals.Calc(dihKind, angles[trial]);
    double distSq = newMol.DihedDist(bond[0], bond[1], bond[2], theta1,
                                     theta, angles[trial]);
    nonbonded_1_4[trial] = data->calc.IntraEnergy_1_4(distSq, prevprev,
                           atom, molIndex);
    if(isnan(nonbonded_1_4[trial]))
      nonbonded_1_4[trial] = num::BIGNUM;

    angleWeights[trial] = exp(-ff.beta * (angleEnergy[trial] +
                                          nonbonded_1_4[trial]));
    stepWeight += angleWeights[trial];
  }
  return stepWeight;
}

double DCLink::GenerateDihedralsOld(TrialMol& oldMol, uint molIndex,
                                    double* angles, double* angleEnergy,
                                    double* angleWeights)
{
  double* nonbonded_1_4 = data->nonbonded_1_4;
  double stepWeight = 0.0;
  PRNG& prng = data->prng;
  const Forcefield& ff = data->ff;
  std::fill_n(nonbonded_1_4, data->nDihTrials, 0.0);

  double theta1 = oldMol.GetTheta(prevprev, prev, focus);

  for (uint trial = 0, count = data->nDihTrials; trial < count; ++trial) {
    angles[trial] = prng.rand(2 * M_PI);
    angleEnergy[trial] = ff.dihedrals.Calc(dihKind, angles[trial]);
    double distSq = oldMol.DihedDist(oldBond[0], oldBond[1], oldBond[2],
                                     theta1, theta, angles[trial]);


    nonbonded_1_4[trial] = data->calc.IntraEnergy_1_4(distSq, prevprev,
                           atom, molIndex);

    if(isnan(nonbonded_1_4[trial]))
      nonbonded_1_4[trial] = num::BIGNUM;

    angleWeights[trial] = exp(-ff.beta * (angleEnergy[trial] +
                                          nonbonded_1_4[trial]));
    stepWeight += angleWeights[trial];
  }
  return stepWeight;
}

void DCLink::UseOldDih(TrialMol& oldMol, uint molIndex, double& energy,
                       double& weight)
{
  double* nonbonded_1_4 = data->nonbonded_1_4;
  PRNG& prng = data->prng;
  const Forcefield& ff = data->ff;
  double beta = data->ff.beta;
  std::fill_n(nonbonded_1_4, data->nDihTrials, 0.0);

  double theta0 = oldMol.GetTheta(prevprev, prev, focus);

  energy = ff.dihedrals.Calc(dihKind, phi);
  double distSq = oldMol.OldDistSq(prevprev, atom);

  nonbonded_1_4[0] = data->calc.IntraEnergy_1_4(distSq, prevprev, atom,
                     molIndex);
  weight = exp(-beta * (energy + nonbonded_1_4[0]));

  for (uint trial = data->nDihTrials - 1; trial-- > 0;) {
    double trialPhi = prng.rand(2 * M_PI);
    double distSq = oldMol.DihedDist(oldBond[0], oldBond[1], oldBond[2],
                                     theta0, theta, trialPhi);

    double tempEn = data->calc.IntraEnergy_1_4(distSq, prevprev, atom, molIndex);
    if(isnan(tempEn))
      tempEn = num::BIGNUM;

    double trialEnergy = ff.dihedrals.Calc(dihKind, trialPhi) + tempEn;
    weight += exp(-beta * trialEnergy);
  }
}

}
