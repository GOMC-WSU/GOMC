/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCCrankShaftDih.h"
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

DCCrankShaftDih::DCCrankShaftDih(DCData* data, const mol_setup::MolKind& kind,
                                uint a0, uint a1, uint a2, uint a3) : 
                                data(data), a0(a0), a1(a1), a2(a2), a3(a3)
{
  using namespace mol_setup;
  using namespace std;
  vector<bool> visited(kind.atoms.size(), false);
  totAtoms = kind.atoms.size();
  //Find all the atoms that bonds with atoms a1
  vector<Bond> bonds = AtomBonds(kind, a1);
  //Remove the a0-a1 bond
  bonds.erase(remove_if(bonds.begin(), bonds.end(), FindA1(a0)), bonds.end());

  //Loop through atoms that are bonded to a1
  for(uint b = 0; b < bonds.size(); b++) {
    //Store the atom index if it doesnot exist
    if(!visited[bonds[b].a0]) {
      atoms.push_back(bonds[b].a0);
      visited[bonds[b].a0] = true;
    }

    vector<Bond> temp = AtomBonds(kind, bonds[b].a1);
    //Remove a2-a3 bonds
    temp.erase(remove_if(temp.begin(), temp.end(), FindA1(a3)), temp.end());
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
}

void DCCrankShaftDih::PrepareOld(TrialMol& oldMol, uint molIndex)
{
  XYZ center = oldMol.AtomPosition(a0);
  for(uint i = 0; i < numAtom; i++) {
    //Unwrap the coordinates with respect to a0.
    XYZ temp = oldMol.AtomPosition(atoms[i]);
    XYZ coord = data->axes.UnwrapPBC(temp, oldMol.GetBox(), center);
    //Shift the atoms to origin
    coord -= center;
    multiPosRotions[i].Set(0, coord); 
  }
}

void DCCrankShaftDih::PrepareNew(TrialMol& newMol, uint molIndex)
{
  XYZ center = newMol.AtomPosition(a0);
  for(uint i = 0; i < numAtom; i++) {
    //Unwrap the coordinates with respect to a0.
    XYZ temp = newMol.AtomPosition(atoms[i]);
    XYZ coord = data->axes.UnwrapPBC(temp, newMol.GetBox(), center);
    //Shift the atoms to origin
    coord -= center;
    multiPosRotions[i].Set(0, coord);  
  }
}

    
void DCCrankShaftDih::BuildOld(TrialMol& oldMol, uint molIndex)
{
  PRNG& prng = data->prng;
  uint nLJTrials = data->nLJTrialsNth;
  double* inter = data->inter;
  double* real = data->real;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);

  //Set up rotation matrix using a0-a3 axis.
  XYZ center = oldMol.AtomPosition(a0);
  XYZ rotAxis = oldMol.AtomPosition(a0) - oldMol.AtomPosition(a3);
  rotAxis = data->axes.MinImage(rotAxis, oldMol.GetBox());
  rotAxis.Normalize();
  RotationMatrix cross = RotationMatrix::CrossProduct(rotAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotAxis);

  //Spin all nLJTrial except the original coordinates
  for (uint lj = nLJTrials; lj-- > 1;) {
    double theta = data->prng.rand(M_PI * 2.0);
    RotationMatrix spin = RotationMatrix::FromAxisAngle(theta, cross, tensor);

    for (uint a = 0; a < numAtom; a++) {
      //find positions
      multiPosRotions[a].Set(lj, spin.Apply(multiPosRotions[a][0]));
      multiPosRotions[a].Add(lj, center);
    }
  }

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
    stepWeight += exp(-1 * data->ff.beta *
                      (inter[trial] + real[trial]));
  }
  oldMol.MultWeight(stepWeight / nLJTrials);
  oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0],
                          0.0, 0.0, 0.0));

  for (uint a = 0; a < totAtoms; a++) {                        
    oldMol.ConfirmOldAtom(a);
  }
}

void DCCrankShaftDih::BuildNew(TrialMol& newMol, uint molIndex)
{
  PRNG& prng = data->prng;
  uint nLJTrials = data->nLJTrialsNth;
  double* inter = data->inter;
  double* real = data->real;
  double* ljWeights = data->ljWeights;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);

  //Set up rotation matrix using a0-a3 axis.
  XYZ center = newMol.AtomPosition(a0);
  XYZ rotAxis = newMol.AtomPosition(a0) - newMol.AtomPosition(a3);
  rotAxis = data->axes.MinImage(rotAxis, newMol.GetBox());
  rotAxis.Normalize();
  RotationMatrix cross = RotationMatrix::CrossProduct(rotAxis);
  RotationMatrix tensor = RotationMatrix::TensorProduct(rotAxis);

  //Go backward to to preserve prototype
  for (uint lj = nLJTrials; lj-- > 0;) {
    double theta = data->prng.rand(M_PI * 2.0);
    RotationMatrix spin = RotationMatrix::FromAxisAngle(theta, cross, tensor);

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
    ljWeights[trial] = exp(-1 * data->ff.beta *
                           (inter[trial] + real[trial]));
    stepWeight += ljWeights[trial];
  }
  uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  newMol.MultWeight(stepWeight / nLJTrials);
  newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner],
                          0.0, 0.0, 0.0));

  for (uint a = 0; a < numAtom; a++) { 
    newMol.AddAtom(atoms[a], multiPosRotions[a][winner]);
  }

  for (uint a = 0; a < totAtoms; a++) {                                                      
     newMol.ConfirmOldAtom(a);
  }
}
}
