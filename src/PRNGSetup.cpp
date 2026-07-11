/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "PRNGSetup.h" //header spec

#include <stdio.h> //for exit

#include "ConfigSetup.h" //For restart, gen kind settings
#include "Reader.h"      //For file I/O
#include "StrLib.h"      //for Compare

namespace prng_setup {
const std::string STEP_STR = "STEP";

void PRNGInitData::HandleError(std::string const &mode) {
  std::cerr << "ERROR: Attempt to initialize generator in mode \"" << mode
            << "\" failed!" << std::endl
            << std::endl;
  exit(1);
}

// Init from file.
void PRNGInitData::Init(Reader &prngSeed, const ulong resStep) {
  Goto(prngSeed, resStep);
  Read(prngSeed);
  prng = new MTRand();
  prng->load(loadArray);
}

void PRNGInitData::Goto(Reader &prngSeed, const ulong resStep) {
  std::string varName;
  ulong step;
  while (prngSeed.Read(varName)) {
    if (str::compare(varName, STEP_STR))
      prngSeed.file >> step;
    if (step == resStep)
      break;
  }
}

void PRNGInitData::Read(Reader &prngSeed) {
  uint temp, i = 0;
  // Init seed array.
  loadArray = new MTRand::uint32[MTRand::N + 1];
  // Read seed;
  for (; i < MTRand::N; i++) {
    prngSeed.file >> temp;
    loadArray[i] = (MTRand::uint32)(temp);
  }
  // Read left
  prngSeed.file >> temp;
  loadArray[i] = (MTRand::uint32)(temp);
}
} // namespace prng_setup

const std::string PRNGSetup::seedFileAlias =
    "Seed file for Mersenne Twister pseudo-random number generator";

void PRNGSetup::Init(config_setup::RestartSettings const &restart,
                     config_setup::PRNGKind const &genConf,
                     std::string const &name) {
  if (genConf.IsRand())
    prngMaker.Init();
  else if (genConf.IsSeed())
    prngMaker.Init(genConf.seed);

  if (prngMaker.prng == NULL)
    prngMaker.HandleError(genConf.kind);
}
