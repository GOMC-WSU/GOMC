/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
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
  real* inter = data->inter;
  real* real_en = data->real_en;
  bool* overlap = data->overlap;
  real stepWeight = 0;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real_en, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  if(oldMol.COMFix()) {
    nLJTrials = 1;
  } else {
    prng.FillWithRandom(positions, nLJTrials, data->axes, oldMol.GetBox());
  }
  positions.Set(0, data->axes.WrapPBC(oldMol.AtomPosition(atom), oldMol.GetBox()));
  data->calc.ParticleInter(inter, real_en, positions, overlap, atom, molIndex,
                           oldMol.GetBox(), nLJTrials);

  for (uint trial = 0; trial < nLJTrials; ++trial) {
    stepWeight += exp(-1 * data->ff.beta *
                      (inter[trial] + real_en[trial]));
  }
  oldMol.UpdateOverlap(overlap[0]);
  oldMol.MultWeight(stepWeight / nLJTrials);
  oldMol.AddEnergy(Energy(0.0, 0.0, inter[0], real_en[0],
                          0.0, 0.0, 0.0));
  oldMol.ConfirmOldAtom(atom);
}

void DCSingle::BuildNew(TrialMol& newMol, uint molIndex)
{
  PRNG& prng = data->prng;
  XYZArray& positions = data->positions;
  uint nLJTrials = data->nLJTrialsFirst;
  real* inter = data->inter;
  real* real_en = data->real_en;
  real* ljWeights = data->ljWeights;
  bool* overlap = data->overlap;

  std::fill_n(inter, nLJTrials, 0.0);
  std::fill_n(real_en, nLJTrials, 0.0);
  std::fill_n(ljWeights, nLJTrials, 0.0);
  std::fill_n(overlap, nLJTrials, false);

  if(newMol.COMFix()) {
    nLJTrials = 1;
    positions.Set(0, data->axes.WrapPBC(newMol.GetCavityCenter(), newMol.GetBox()));
  } else {
    prng.FillWithRandom(positions, nLJTrials, data->axes, newMol.GetBox());
  }
  data->calc.ParticleInter(inter, real_en, positions, overlap, atom, molIndex,
                           newMol.GetBox(), nLJTrials);

  real stepWeight = 0;
  for (uint trial = 0; trial < nLJTrials; ++trial) {
    ljWeights[trial] = exp(-1 * data->ff.beta *
                           (inter[trial] + real_en[trial]));
    stepWeight += ljWeights[trial];
  }
  uint winner = prng.PickWeighted(ljWeights, nLJTrials, stepWeight);
  newMol.UpdateOverlap(overlap[winner]);
  newMol.MultWeight(stepWeight / nLJTrials);
  newMol.AddEnergy(Energy(0.0, 0.0, inter[winner], real_en[winner],
                          0.0, 0.0, 0.0));
  newMol.AddAtom(atom, positions[winner]);
}
}
