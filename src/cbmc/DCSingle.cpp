/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "DCSingle.h"
#include "TrialMol.h"
#include "DCData.h"
#include "PRNG.h"
#include "CalculateEnergy.h"
#include "XYZArray.h"
#include "Forcefield.h"

namespace cbmc
{
DCSingle::DCSingle(DCData* data, uint atom) : data(data), atom(atom)
{
  if(data->nLJTrialsFirst < 1) {
    std::cout << "Error: CBMC first atom trials must be greater than 0.\n";
    exit(EXIT_FAILURE);
  }
}

void DCSingle::BuildOld(TrialMol& oldMol, uint molIndex)
{
  PRNG& prng = data->prng;
  XYZArray& positions = data->positions;
  uint nLJTrials = data->nLJTrialsFirst;
  double* inter = data->inter;
  double* real = data->real;
  bool* overlap = data->overlap;
  double stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  if(oldMol.COMFix()) {
    nLJTrials = 1;
    std::cerr << "Error: In DCSingle::entered comfix branch" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    prng.FillWithRandom(positions, nLJTrials, data->axes, oldMol.GetBox());
  }
  positions.Set(0, data->axes.WrapPBC(oldMol.AtomPosition(atom), oldMol.GetBox()));
  data->calc.ParticleInter(inter, real, positions, overlap, atom, molIndex,
                           oldMol.GetBox(), nLJTrials);

  for (uint trial = 0; trial < nLJTrials; ++trial) {
    stepWeight += exp(-1 * data->ff.beta *
                      (inter[trial] + real[trial]));
  }
  oldMol.UpdateOverlap(overlap[0]);
  oldMol.MultWeight(stepWeight / nLJTrials);
  oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real[0],
                          0.0, 0.0, 0.0));
  oldMol.ConfirmOldAtom(atom);
}

void DCSingle::BuildNew(TrialMol& newMol, uint molIndex)
{
  PRNG& prng = data->prng;
  XYZArray& positions = data->positions;
  uint nLJTrials = data->nLJTrialsFirst;
  double* inter = data->inter;
  double* real = data->real;
  double* ljWeights = data->ljWeights;
  bool* overlap = data->overlap;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  if(newMol.COMFix()) {
    nLJTrials = 1;
    positions.Set(0, data->axes.WrapPBC(newMol.GetCavityCenter(), newMol.GetBox()));
    std::cerr << "Error: In DCSingle::entered comfix branch" << std::endl;
    exit(EXIT_FAILURE);
  } else {
    prng.FillWithRandom(positions, nLJTrials, data->axes, newMol.GetBox());
  }
  data->calc.ParticleInter(inter, real, positions, overlap, atom, molIndex,
                           newMol.GetBox(), nLJTrials);

  double stepWeight = 0;
  for (uint trial = 0; trial < nLJTrials; ++trial) {
    ljWeights[trial] = exp(-1 * data->ff.beta *
                           (inter[trial] + real[trial]));
    stepWeight += ljWeights[trial];
  }
  uint winner = prng.PickWeightedSpecial(ljWeights, nLJTrials, stepWeight, inter, real, positions, overlap, atom, molIndex,
                           newMol.GetBox(), nLJTrials, newMol.COMFix());
  newMol.UpdateOverlap(overlap[winner]);
  newMol.MultWeight(stepWeight / nLJTrials);
  newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real[winner],
                          0.0, 0.0, 0.0));
  newMol.AddAtom(atom, positions[winner]);
}
}
