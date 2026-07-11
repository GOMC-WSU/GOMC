/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef PRNG_SETUP_H
#define PRNG_SETUP_H

#include "BasicTypes.h" //For uint/ulong
#include "MersenneTwister.h"

class Reader;
namespace config_setup {
struct RestartSettings;
struct PRNGKind;
} // namespace config_setup

namespace prng_setup {
struct PRNGInitData {
  PRNGInitData() : loadArray(NULL) {}
  // WARNING/NOTE: Object is passed off, NOT deleted.
  ~PRNGInitData(void) { delete[] loadArray; }

  // Init from file.
  void Init(Reader &prngSeed, const ulong resStep);

  // Init from unsigned integer seed.
  void Init(const uint seed) { prng = new MTRand((MTRand::uint32)(seed)); }

  // Init pure random.
  void Init(void) { prng = new MTRand(); }

  void HandleError(std::string const &mode);

  MTRand *prng;

private:
  // Goto correct step in seed dump file.
  void Goto(Reader &prngSeed, const ulong resStep);

  // Read in saved seed;
  void Read(Reader &prngSeed);

  MTRand::uint32 *loadArray;
};
} // namespace prng_setup

class PRNGSetup {
public:
  prng_setup::PRNGInitData prngMaker;
  void Init(config_setup::RestartSettings const &restart,
            config_setup::PRNGKind const &genConf, std::string const &name);

  static const std::string seedFileAlias;
};

#endif /*PRNG_SETUP_H*/
