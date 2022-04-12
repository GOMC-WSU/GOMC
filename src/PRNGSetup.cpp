/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in the COPYRIGHT.txt
along with this program, also can be found at <https://opensource.org/licenses/MIT>.
********************************************************************************/
#include <stdio.h> //for exit

#include "Reader.h" //For file I/O
#include "PRNGSetup.h" //header spec
#include "ConfigSetup.h" //For restart, gen kind settings
#include "StrLib.h" //for Compare

namespace prng_setup
{
const std::string STEP_STR = "STEP";

void PRNGInitData::HandleError(std::string const& mode)
{
  std::cerr << "ERROR: Attempt to initialize generator in mode \""
            << mode << "\" failed!" << std::endl << std::endl;
  exit(1);
}

//Init from file.
void PRNGInitData::Init(Reader & prngSeed, const ulong resStep)
{
  Goto(prngSeed, resStep);
  Read(prngSeed);
  prng = new MTRand();
  prng->load(loadArray);
}

void PRNGInitData::Goto(Reader & prngSeed, const ulong resStep)
{
  std::string varName;
  ulong step;
  while (prngSeed.Read(varName)) {
    if (str::compare(varName, STEP_STR))
      prngSeed.file >> step;
    if (step == resStep) break;
  }
}

void PRNGInitData::Read(Reader & prngSeed)
{
  uint temp, i = 0;
  //Init seed array.
  loadArray = new MTRand::uint32[MTRand::N + 1];
  //Read seed;
  for (; i < MTRand::N; i++) {
    prngSeed.file >> temp;
    loadArray[i] = (MTRand::uint32)(temp);
  }
  //Read left
  prngSeed.file >> temp;
  loadArray[i] = (MTRand::uint32)(temp);
}
}

const std::string PRNGSetup::seedFileAlias =
  "Seed file for Mersenne Twister pseudo-random number generator";

void PRNGSetup::Init(config_setup::RestartSettings const& restart,
                     config_setup::PRNGKind const& genConf,
                     std::string const& name)
{
  if (genConf.IsRand())
    prngMaker.Init();
  else if (genConf.IsSeed())
    prngMaker.Init(genConf.seed);
  else if (genConf.IsRestart()) {
    Reader prngSeedFile(name, seedFileAlias);
    prngMaker.Init(prngSeedFile, restart.step);
  }

  if (prngMaker.prng == NULL)
    prngMaker.HandleError(genConf.kind);
}
