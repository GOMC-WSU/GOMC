/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCCrankShaftAng.h"
#include "TrialMol.h"
#include "PRNG.h"
#include "MolSetup.h"
#include "Forcefield.h"
#include "NumLib.h"
#include "CalculateEnergy.h"
#include "XYZArray.h"
#include <numeric>
#include <cassert>

namespace cbmc
{
  struct FindA1 {
  FindA1(uint x) : x(x) {};
  bool operator()(const mol_setup::Bond& b) {
    return (b.a1 == x);
  }
  uint x;
};

DCCrankShaftAng::DCCrankShaftAng(DCData* data, const mol_setup::MolKind& kind,
                                uint a0, uint a1, uint a2) : 
                                data(data), a0(a0), a1(a1), a2(a2)
{
  using namespace mol_setup;
  using namespace std;
  vector<bool> visited(kind.atoms.size(), false);
  totAtoms = kind.atoms.size();
  //Find all the atoms that bonds with atoms a1
  vector<Bond> bonds = AtomBonds(kind, a1);
  //Remove the a0-a1 and a1-a2 bond
  bonds.erase(remove_if(bonds.begin(), bonds.end(), FindA1(a0)), bonds.end());
  bonds.erase(remove_if(bonds.begin(), bonds.end(), FindA1(a2)), bonds.end());
  //Store the a1 index
  atoms.push_back(a1);
  visited[a1] = true;

  //Loop through other atoms that are bonded to a1
  for(uint b = 0; b < bonds.size(); b++) {
    //Store the atom index if it doesnot exist
    if(!visited[bonds[b].a0]) {
      atoms.push_back(bonds[b].a0);
      visited[bonds[b].a0] = true;
    }

    vector<Bond> temp = AtomBonds(kind, bonds[b].a1);
    for(uint i = 0; i < temp.size(); i++) {
      if(!visited[temp[i].a0]) {
        bonds.push_back(temp[i]);
      }
    }
  }

  numAtom = atoms.size();

  multiPosRotions = new XYZArray[numAtom];
  for(uint i = 0; i < numAtom; ++i) {
    multiPosRotions[i] = XYZArray(data->nLJTrialsNth);
  }

  if(data->nLJTrialsNth < 1) {
      std::cout << "Error: CBMC secondary atom trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }

  if(data->nDihTrials < 1) {
    std::cout << "Error: CBMC dihedral trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

void DCCrankShaftAng::PrepareOld(TrialMol& oldMol, uint molIndex)
{
  for (uint a = 0; a < totAtoms; a++) {                        
    oldMol.ConfirmOldAtom(a);
  }

  XYZ center = oldMol.AtomPosition(a0);
  for(uint i = 0; i < numAtom; i++) {
    oldMol.UnConfirmOldAtom(atoms[i]);
    //Unwrap the coordinates with respect to a0.
    XYZ temp = oldMol.AtomPosition(atoms[i]);
    XYZ coord = data->axes.UnwrapPBC(temp, oldMol.GetBox(), center);
    //Shift the atoms to origin
    coord -= center;
    multiPosRotions[i].Set(0, coord); 
  }
  //Calculate ol bonded and intrabonded energy wihtout atoms
  oldEnergy = data->calc.MoleculeIntra(oldMol, molIndex);
}

void DCCrankShaftAng::PrepareNew(TrialMol& newMol, uint molIndex)
{
  for (uint a = 0; a < totAtoms; a++) {                        
    newMol.ConfirmOldAtom(a);
  }

  XYZ center = newMol.AtomPosition(a0);
  for(uint i = 0; i < numAtom; i++) {
    newMol.UnConfirmOldAtom(atoms[i]);
    //Unwrap the coordinates with respect to a0.
    XYZ temp = newMol.AtomPosition(atoms[i]);
    XYZ coord = data->axes.UnwrapPBC(temp, newMol.GetBox(), center);
    //Shift the atoms to origin
    coord -= center;
    multiPosRotions[i].Set(0, coord);  
  }
  //Calculate ol bonded and intrabonded energy wihtout atoms
  oldEnergy = data->calc.MoleculeIntra(newMol, molIndex);
}

    
void DCCrankShaftAng::BuildOld(TrialMol& oldMol, uint molIndex)
{
  PRNG& prng = data->prng;
  uint nLJTrials = data->nLJTrialsNth;
  uint nDihTrials = data->nDihTrials;
  double* torsion = data->angles;
  double* torWeights = data->angleWeights;
  double* torEnergy = data->angleEnergy;
  double* bondedEn = data->bonded;
  double* nonbonded = data->nonbonded;
  double* intraNonbonded = data->nonbonded_1_4;
  double* ljWeights = data->ljWeights;
  double* inter = data->inter;
  double* real = data->real;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);

  //Set up rotation matrix using a0-a3 axis.
  XYZ center = oldMol.AtomPosition(a0);
  XYZ rotAxis = oldMol.AtomPosition(a0) - oldMol.AtomPosition(a2);
  rotAxis = data->axes.MinImage(rotAxis, oldMol.GetBox());
  rotAxis.Normalize();
  RotationMatrix cross = RotationMatrix::CrossProduct(rotAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotAxis);

  //Spin all nLJTrial except the original coordinates
  for (uint lj = nLJTrials; lj-- > 1;) {
    ChooseTorsion(oldMol, molIndex, cross, tensor);
    ljWeights[lj] = std::accumulate(torWeights,
                                    torWeights + nDihTrials, 0.0);
    uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
    bondedEn[lj] = torEnergy[winner];
    nonbonded[lj] = intraNonbonded[winner];
    //convert chosen torsion to 3D positions
    RotationMatrix spin = RotationMatrix::FromAxisAngle(torsion[winner],
                          cross, tensor);

    for (uint a = 0; a < numAtom; a++) {
      //find positions
      multiPosRotions[a].Set(lj, spin.Apply(multiPosRotions[a][0]));
      multiPosRotions[a].Add(lj, center);
    }
  }

  ChooseTorsionOld(oldMol, molIndex, cross, tensor);
  ljWeights[0] = std::accumulate(torWeights, torWeights + nDihTrials, 0.0);
  bondedEn[0] = torEnergy[0];
  nonbonded[0] = intraNonbonded[0];

  for (uint a = 0; a < numAtom; a++) {
    //Shift original coordinate back.
    multiPosRotions[a].Add(0, center);
    //Wrap the atom coordinates
    data->axes.WrapPBC(multiPosRotions[a], oldMol.GetBox());
    //Calculate nonbonded energy
    data->calc.ParticleInter(inter, real, multiPosRotions[a], atoms[a], molIndex,
                             oldMol.GetBox(), nLJTrials);
  }

  for (uint trial = 0; trial < nLJTrials; ++trial) {
    ljWeights[trial] *= exp(-1 * data->ff.beta *
                           (inter[trial] + real[trial]));
    stepWeight += ljWeights[trial];
  }
  oldMol.MultWeight(stepWeight / nLJTrials);
  oldMol.AddEnergy(Energy(oldEnergy.intraBond + bondedEn[0], 
                          oldEnergy.intraNonbond + nonbonded[0],
                          inter[0], real[0],
                          0.0, 0.0, 0.0));
}

void DCCrankShaftAng::BuildNew(TrialMol& newMol, uint molIndex)
{
  PRNG& prng = data->prng;
  uint nLJTrials = data->nLJTrialsNth;
  uint nDihTrials = data->nDihTrials;
  double* torsion = data->angles;
  double* torWeights = data->angleWeights;
  double* torEnergy = data->angleEnergy;
  double* bondedEn = data->bonded;
  double* nonbonded = data->nonbonded;
  double* intraNonbonded = data->nonbonded_1_4;
  double* ljWeights = data->ljWeights;
  double* inter = data->inter;
  double* real = data->real;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);

  //Set up rotation matrix using a0-a3 axis.
  XYZ center = newMol.AtomPosition(a0);
  XYZ rotAxis = newMol.AtomPosition(a0) - newMol.AtomPosition(a2);
  rotAxis = data->axes.MinImage(rotAxis, newMol.GetBox());
  rotAxis.Normalize();
  RotationMatrix cross = RotationMatrix::CrossProduct(rotAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotAxis);

  //Go backward to to preserve prototype
  for (uint lj = nLJTrials; lj-- > 0;) {
    ChooseTorsion(newMol, molIndex, cross, tensor);
    ljWeights[lj] = std::accumulate(torWeights,
                                    torWeights + nDihTrials, 0.0);
    uint winner = prng.PickWeighted(torWeights, nDihTrials, ljWeights[lj]);
    bondedEn[lj] = torEnergy[winner];
    nonbonded[lj] = intraNonbonded[winner];
    //convert chosen torsion to 3D positions
    RotationMatrix spin = RotationMatrix::FromAxisAngle(torsion[winner],
                          cross, tensor);

    for (uint a = 0; a < numAtom; a++) {
      //find positions
      multiPosRotions[a].Set(lj, spin.Apply(multiPosRotions[a][0]));
      multiPosRotions[a].Add(lj, center);
    }
  }

  for (uint a = 0; a < numAtom; a++) {
    //Wrap the atom coordinates
    data->axes.WrapPBC(multiPosRotions[a], newMol.GetBox());
    //Calculate nonbonded energy
    data->calc.ParticleInter(inter, real, multiPosRotions[a], atoms[a], molIndex,
                             newMol.GetBox(), nLJTrials);
  }
  
  for (uint trial = 0; trial < nLJTrials; trial++) {
    ljWeights[trial] *= exp(-1 * data->ff.beta *
                           (inter[trial] + real[trial]));
    stepWeight += ljWeights[trial];
  }
  uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  newMol.MultWeight(stepWeight / nLJTrials);
  newMol.AddEnergy(Energy(oldEnergy.intraBond + bondedEn[winner],
                          oldEnergy.intraNonbond + nonbonded[winner],
                          inter[winner], real[winner],
                          0.0, 0.0, 0.0));

  for (uint a = 0; a < numAtom; a++) { 
    newMol.AddAtom(atoms[a], multiPosRotions[a][winner]);
  }

}

void DCCrankShaftAng::ChooseTorsion(TrialMol& mol, uint molIndex,
                                    RotationMatrix& cross,
                                    RotationMatrix& tensor)
{ 
  uint nDihTrials = data->nDihTrials;
  double* torsion = data->angles;
  double* torWeights = data->angleWeights;
  double* torEnergy = data->angleEnergy;
  double* intraNonbonded = data->nonbonded_1_4;

  XYZ center = mol.AtomPosition(a0);
  for (uint tor = 0; tor < nDihTrials; ++tor) {
    torsion[tor] = data->prng.rand(M_PI * 2);
    //convert chosen torsion to 3D positions
    RotationMatrix spin = RotationMatrix::FromAxisAngle(torsion[tor],
                          cross, tensor);
    for (uint a = 0; a < numAtom; a++) {
      XYZ coord = spin.Apply(multiPosRotions[a][0]);
      mol.AddAtom(atoms[a], coord + center);
    }
    Energy en = data->calc.MoleculeIntra(mol, molIndex);
    //Get bonded energy difference due to rotating of atoms
    en -= oldEnergy;
    torEnergy[tor] = en.intraBond;
    intraNonbonded[tor] = en.intraNonbond;
    torWeights[tor] = exp(-1 * data->ff.beta * (torEnergy[tor] + intraNonbonded[tor]));
  }
}

void DCCrankShaftAng::ChooseTorsionOld(TrialMol& mol, uint molIndex,
                                      RotationMatrix& cross,
                                      RotationMatrix& tensor)
{ 
  uint nDihTrials = data->nDihTrials;
  double* torsion = data->angles;
  double* torWeights = data->angleWeights;
  double* torEnergy = data->angleEnergy;
  double* intraNonbonded = data->nonbonded_1_4;

  XYZ center = mol.AtomPosition(a0);
  for (uint tor = 0; tor < nDihTrials; ++tor) {
    //Use actual coordinate fir first torsion trial
    torsion[tor] = (tor == 0) ? 0.0 : data->prng.rand(M_PI * 2);
    //convert chosen torsion to 3D positions
    RotationMatrix spin = RotationMatrix::FromAxisAngle(torsion[tor],
                          cross, tensor);
    for (uint a = 0; a < numAtom; a++) {
      XYZ coord = spin.Apply(multiPosRotions[a][0]);
      mol.AddAtom(atoms[a], coord + center);
    }
    Energy en = data->calc.MoleculeIntra(mol, molIndex);
    //Get bonded energy difference due to rotating of atoms
    en -= oldEnergy;
    torEnergy[tor] = en.intraBond;
    intraNonbonded[tor] = en.intraNonbond;
    torWeights[tor] = exp(-1 * data->ff.beta * (torEnergy[tor] + intraNonbonded[tor]));
  }
}


}
