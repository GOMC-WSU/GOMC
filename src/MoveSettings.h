/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOVE_SETTINGS_H
#define MOVE_SETTINGS_H

#include "EnsemblePreprocessor.h" //For BOX_TOTAL
#include "BasicTypes.h"           //for uint
#include "OutputVars.h"
#include "PDBSetup.h" //Primary source of volume.
#include "MoveConst.h"           //For sizes of arrays.
#include <vector>

class StaticVals;                 //For various initialization constants.
class BoxDimensions;              //For axis sizes

using namespace std;

class MoveSettings
{
public:
  friend class OutputVars;
  MoveSettings(BoxDimensions & dim) : boxDimRef(dim)
  {
    acceptPercent.resize(BOX_TOTAL);
    scale.resize(BOX_TOTAL);
    accepted.resize(BOX_TOTAL);
    tries.resize(BOX_TOTAL);
    tempAccepted.resize(BOX_TOTAL);
    tempTries.resize(BOX_TOTAL);
    for(uint b = 0; b < BOX_TOTAL; b++) {
      acceptPercent[b].resize(mv::MOVE_KINDS_TOTAL);
      scale[b].resize(mv::MOVE_KINDS_TOTAL);
      accepted[b].resize(mv::MOVE_KINDS_TOTAL);
      tries[b].resize(mv::MOVE_KINDS_TOTAL);
      tempAccepted[b].resize(mv::MOVE_KINDS_TOTAL);
      tempTries[b].resize(mv::MOVE_KINDS_TOTAL);
    }
  }

  MoveSettings& operator=(MoveSettings const& rhs)
  {
    return *this;
  }

  void Init(StaticVals const& statV, pdb_setup::Remarks const& remarks,
            const uint tkind);

  void Update(const uint move, const bool isAccepted, const uint step,
              const uint box, const uint kind = 0);

  void AdjustMoves(const uint step);

  void Adjust(const uint box, const uint move, const uint kind);

  real Scale(const uint box, const uint move, const uint kind = 0) const
  {
    return scale[box][move][kind];
  }

  real GetAccept(const uint box, const uint move, const uint kind = 0) const
  {
    return acceptPercent[box][move][kind];
  }

  real GetTrial(const uint box, const uint move, const uint kind = 0) const
  {
    return tries[box][move][kind];
  }

  uint GetAcceptTot(const uint box, const uint move) const;
  uint GetTrialTot(const uint box, const uint move) const;
  real GetScaleTot(const uint box, const uint move) const;

private:

  vector< vector< vector<real> > > scale, acceptPercent;
  vector< vector< vector<uint> > > accepted, tries, tempAccepted, tempTries;

  uint perAdjust;
  uint totKind;

#if ENSEMBLE == GEMC
  uint GEMC_KIND;
#endif

  BoxDimensions & boxDimRef;

  static const real TARGET_ACCEPT_FRACT;
  static const real TINY_AMOUNT;

  // make checkopintoutput and checkpointsetup a friend class to have access to
  // private data
  friend class CheckpointOutput;
  friend class CheckpointSetup;
};


#endif /*MOVE_SETTINGS_H*/
