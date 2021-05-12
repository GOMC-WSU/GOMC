/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef MOVE_SETTINGS_H
#define MOVE_SETTINGS_H

#include "EnsemblePreprocessor.h" //For BOX_TOTAL
#include "BasicTypes.h"           //For uint
#include "OutputVars.h"
#include "PDBSetup.h"             //Primary source of volume.
#include "MoveConst.h"            //For array sizes
#include <vector>

namespace mp
{
const int MPDISPLACE = 0;
const int MPROTATE = 1;
const int MPTOTALTYPES = 2;
const double TARGET_ACCEPT_FRACT = 0.3;
}

class StaticVals;                 //For various initialization constants.
class BoxDimensions;              //For axis sizes

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
    mp_r_max.resize(BOX_TOTAL);
    mp_t_max.resize(BOX_TOTAL);
    mp_accepted.resize(BOX_TOTAL);
    mp_tries.resize(BOX_TOTAL);
    mp_interval_accepted.resize(BOX_TOTAL);
    mp_interval_tries.resize(BOX_TOTAL);
    for(uint b = 0; b < BOX_TOTAL; b++) {
      acceptPercent[b].resize(mv::MOVE_KINDS_TOTAL);
      scale[b].resize(mv::MOVE_KINDS_TOTAL);
      accepted[b].resize(mv::MOVE_KINDS_TOTAL);
      tries[b].resize(mv::MOVE_KINDS_TOTAL);
      tempAccepted[b].resize(mv::MOVE_KINDS_TOTAL);
      tempTries[b].resize(mv::MOVE_KINDS_TOTAL);
      mp_accepted[b].resize(mp::MPTOTALTYPES);
      mp_tries[b].resize(mp::MPTOTALTYPES);
      mp_interval_accepted[b].resize(mp::MPTOTALTYPES);
      mp_interval_tries[b].resize(mp::MPTOTALTYPES);
    }
  }

  void Init(StaticVals const& statV, pdb_setup::Remarks const& remarks,
            const uint tkind);

  void Update(const uint move, const bool isAccepted,
              const uint box, const uint kind = 0);

  void AdjustMoves(const ulong step);

  void Adjust(const uint box, const uint move, const uint kind);

  void AdjustMultiParticle(const uint box, const uint typePick);

  void UpdateMoveSettingMultiParticle(uint box, bool isAccept, uint typePick);

  inline double Scale(const uint box, const uint move, const uint kind = 0) const
  {
    return scale[box][move][kind];
  }

  inline double GetAccept(const uint box, const uint move, const uint kind = 0) const
  {
    return acceptPercent[box][move][kind];
  }

  inline double GetTrial(const uint box, const uint move, const uint kind = 0) const
  {
    return tries[box][move][kind];
  }

  inline double GetRMAX(const uint box) const
  {
    return mp_r_max[box];
  }

  inline double GetTMAX(const uint box) const
  {
    return mp_t_max[box];
  }

  uint GetAcceptTot(const uint box, const uint move) const;
  uint GetTrialTot(const uint box, const uint move) const;
  double GetScaleTot(const uint box, const uint move) const;
  
  inline bool GetSingleMoveAccepted(uint box) const
  {
    return isSingleMoveAccepted[box];
  }
  
  inline void SetSingleMoveAccepted(uint box)
  {
    isSingleMoveAccepted[box] = true;
  }

  inline void UnsetSingleMoveAccepted(uint box)
  {
    isSingleMoveAccepted[box] = false;
  }

private:

  std::vector< std::vector< std::vector<double> > > scale, acceptPercent;
  std::vector< std::vector< std::vector<uint> > > accepted, tries, tempAccepted, tempTries;
  std::vector< std::vector< uint > > mp_accepted, mp_tries, mp_interval_accepted, mp_interval_tries;
  std::vector< double > mp_r_max;
  std::vector< double > mp_t_max;

  uint perAdjust;
  uint totKind;
  bool isSingleMoveAccepted[BOXES_WITH_U_NB];

#if ENSEMBLE == GEMC
  uint GEMC_KIND;
#endif

  BoxDimensions & boxDimRef;

  static const double TARGET_ACCEPT_FRACT;
  static const double TINY_AMOUNT;
  //The alpha values are used to control how much weight is placed on the number of accepted
  //moves from the most recent adjustment period versus on the total number of accepted moves
  static const double r_alpha, t_alpha;
  //If the MultiParticle acceptance percentage is within mp_accept_tol, we don't adjust the max
  static const double mp_accept_tol;

  // make CheckpointOutput and CheckpointSetup friend classes to have access to
  // private data
  friend class CheckpointOutput;
  friend class CheckpointSetup;
};

#endif /*MOVE_SETTINGS_H*/
