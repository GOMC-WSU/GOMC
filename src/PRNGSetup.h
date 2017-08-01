/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.0
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PRNG_SETUP_H
#define PRNG_SETUP_H

#include "MersenneTwister.h"
#include "BasicTypes.h" //For uint/ulong

class Reader;
namespace config_setup {
struct RestartSettings;
struct PRNGKind;
}

namespace prng_setup
{
   struct PRNGInitData
   {
      PRNGInitData() : loadArray(NULL) {}
      //WARNING/NOTE: Object is passed off, NOT deleted.
      ~PRNGInitData(void) { delete[] loadArray; }

      //Init from file.
      void Init(Reader & prngSeed, const ulong resStep);

      //Init from unsigned integer seed.
      void Init(const uint seed) { prng = new MTRand((MTRand::uint32)(seed)); }

      //Init pure random.
      void Init(void) { prng = new MTRand(); }

      void HandleError(std::string const& mode);


      MTRand * prng;
    private:

      //Goto correct step in seed dump file.
      void Goto(Reader & prngSeed, const ulong resStep);

      //Read in saved seed;
      void Read(Reader & prngSeed);

      MTRand::uint32 * loadArray;
   };
}

class PRNGSetup
{
public:
   prng_setup::PRNGInitData prngMaker;
   void Init(config_setup::RestartSettings const& restart,
	     config_setup::PRNGKind const& genConf,
	     std::string const& name);

   static const std::string seedFileAlias;
};

#endif /*PRNG_SETUP_H*/
