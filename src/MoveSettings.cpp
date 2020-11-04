/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoveSettings.h"    //header spec
#include "BoxDimensions.h"   //For axis sizes
#include "BoxDimensionsNonOrth.h"
#include "StaticVals.h"      //For init info
#include <cmath>
#include "NumLib.h"          //For bounding functions
#include "GeomLib.h"         //For M_PI

const double MoveSettings::TARGET_ACCEPT_FRACT = 0.5;
const double MoveSettings::TINY_AMOUNT = 1.0e-7;
const double MoveSettings::r_alpha = 0.2;
const double MoveSettings::t_alpha = 0.2;

void MoveSettings::Init(StaticVals const& statV,
                        pdb_setup::Remarks const& remarks,
                        const uint tkind)
{
  isSingleMoveAccepted = true;
  totKind = tkind;
  perAdjust = statV.simEventFreq.perAdjust;
  for(uint b = 0; b < BOX_TOTAL; b++) {
    for(uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++) {
      acceptPercent[b][m].resize(totKind, 0);
      scale[b][m].resize(totKind, 0);
      accepted[b][m].resize(totKind, 0);
      tries[b][m].resize(totKind, 0);
      tempAccepted[b][m].resize(totKind, 0);
      tempTries[b][m].resize(totKind, 0);
    }
  }


  for(uint b = 0; b < BOX_TOTAL; b++) {
    for(uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++) {
      for(uint k = 0; k < totKind; k++) {
        if(m == mv::DISPLACE) {
          if(remarks.restart && remarks.disp[b] > 0.0) {
            scale[b][m][k] = remarks.disp[b];
          } else {
            scale[b][m][k] = boxDimRef.axis.Min(b) * 0.25;
          }
        } else if (m == mv::ROTATE) {
          if(remarks.restart && remarks.rotate[b] > 0.0) {
            scale[b][m][k] = remarks.rotate[b];
          } else {
            scale[b][m][k] = M_PI_4;
          }
        }
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
        else if (m == mv::VOL_TRANSFER) {
          if(remarks.restart && remarks.vol[b] > 0.0) {
            scale[b][m][k] = remarks.vol[b];
          } else {
            scale[b][m][k] = 500;
          }
        }
#endif
      }
    }
  }

  // Initialize MultiParticle settings
  for(int b = 0; b < BOX_TOTAL; b++) {
    mp_r_max[b] = 0.01 * M_PI;
    mp_t_max[b] = 0.02;
    for(int m = 0; m < mp::MPTOTALTYPES; m++) {
      mp_tries[b][m] = 0;
      mp_accepted[b][m] = 0;
      mp_interval_tries[b][m] = 0;
      mp_interval_accepted[b][m] = 0;
    }
  }
}

//Process results of move we just did in terms of acceptance counters
void MoveSettings::Update(const uint move, const bool isAccepted,
                          const uint step, const uint box, const uint kind)
{
  tries[box][move][kind]++;
  tempTries[box][move][kind]++;
  if(isAccepted) {
    tempAccepted[box][move][kind]++;
    accepted[box][move][kind]++;

    if(move != mv::MULTIPARTICLE || move != mv::MULTIPARTICLE_BM)
      isSingleMoveAccepted = true;
  }

  if(move == mv::MULTIPARTICLE || move != mv::MULTIPARTICLE_BM)
    isSingleMoveAccepted = false;

  acceptPercent[box][move][kind] = (double)(accepted[box][move][kind]) /
                                   (double)(tries[box][move][kind]);

  //for any move that we don't care about kind of molecule, it should be included
  //in the if condition
  if (move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE || move == mv::MULTIPARTICLE_BM
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
      || move == mv::MEMC
#endif
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
      || move == mv::VOL_TRANSFER
#endif
     ) {
    for (uint k = 1; k < totKind; k++) {
      tries[box][move][k]++;
      tempTries[box][move][k]++;
      if(isAccepted) {
        tempAccepted[box][move][k]++;
        accepted[box][move][k]++;
      }
      acceptPercent[box][move][k] = (double)(accepted[box][move][k]) /
                                    (double)(tries[box][move][k]);
    }
  }
}

void MoveSettings::AdjustMoves(const uint step)
{
  //Check whether we need to adjust this move's scaling.
  if ((step + 1) % perAdjust == 0) {
    for(uint b = 0; b < BOX_TOTAL; b++) {
      for(uint m = 0; m < mv::MOVE_KINDS_TOTAL; ++m) {
        for(uint k = 0; k < totKind; k++) {
          Adjust(b, m, k);
        }
      }
      for(int m = 0; m < mp::MPTOTALTYPES; m++) {
        AdjustMultiParticle(b, m);
      }
    }
  }
}

void MoveSettings::AdjustMultiParticle(const uint box, const uint typePick)
{
  //Update totals for the entire simulation
  for(int m = 0; m < mp::MPTOTALTYPES; ++m) {
    mp_tries[box][m] += mp_interval_tries[box][m];
    mp_accepted[box][m] += mp_interval_accepted[box][m];
  }
  //Make sure we tried some displacements, otherwise mp_t_max will be nan
  if((double)mp_interval_tries[box][mp::MPDISPLACE] != 0) {
    double fractOfTotalAccept = ((double)mp_accepted[box][mp::MPDISPLACE] /
                                 (double)mp_tries[box][mp::MPDISPLACE]) / mp::TARGET_ACCEPT_FRACT;
    double fractOfIntervalAccept = ((double)mp_interval_accepted[box][mp::MPDISPLACE] /
                                  (double)mp_interval_tries[box][mp::MPDISPLACE]) / mp::TARGET_ACCEPT_FRACT;
    if (fractOfIntervalAccept > 0.0) {
      mp_t_max[box] *= ((1.0-t_alpha) * fractOfTotalAccept + t_alpha * fractOfIntervalAccept);
    } else {
      mp_t_max[box] *= 0.5;
    }
    num::Bound<double>(mp_t_max[box], 0.001, (boxDimRef.axis.Min(box) * 0.5) - 0.001);
    //Reset totals for next adjustment period
    mp_interval_accepted[box][mp::MPDISPLACE] = 0;
    mp_interval_tries[box][mp::MPDISPLACE] = 0;
  }

  //Make sure we tried some rotations, otherwise mp_r_max will be nan
  if((double)mp_interval_tries[box][mp::MPROTATE] != 0) {
    double fractOfTotalAccept = ((double)mp_accepted[box][mp::MPROTATE] /
                                 (double)mp_tries[box][mp::MPROTATE]) / mp::TARGET_ACCEPT_FRACT;
    double fractOfIntervalAccept = ((double)mp_interval_accepted[box][mp::MPROTATE] /
                                  (double)mp_interval_tries[box][mp::MPROTATE]) / mp::TARGET_ACCEPT_FRACT;
    if (fractOfIntervalAccept > 0.0) {
      mp_r_max[box] *= ((1.0-r_alpha) * fractOfTotalAccept + r_alpha * fractOfIntervalAccept);
    } else {
      mp_r_max[box] *= 0.5;
    }
    num::Bound<double>(mp_r_max[box], 0.001, M_PI - 0.001);
    //Reset totals for next adjustment period
    mp_interval_accepted[box][mp::MPROTATE] = 0;
    mp_interval_tries[box][mp::MPROTATE] = 0;
  }
}

void MoveSettings::UpdateMoveSettingMultiParticle(const uint box, bool isAccept,
    const uint typePick)
{
  mp_interval_tries[box][typePick]++;
  if(isAccept) {
    mp_interval_accepted[box][typePick]++;
  }
}

//Adjust responsibly
void MoveSettings::Adjust(const uint box, const uint move, const uint kind)
{
  if(move == mv::DISPLACE) {
    if(tempTries[box][move][kind] > 0) {
      double currentAccept = (double)(tempAccepted[box][move][kind]) /
                             (double)(tempTries[box][move][kind]);
      double fractOfTargetAccept = currentAccept / TARGET_ACCEPT_FRACT;
      if (fractOfTargetAccept > 0.0) {
        scale[box][move][kind] *= fractOfTargetAccept;
      } else {
        scale[box][move][kind] *= 0.5;
      }
    }
    num::Bound<double>(scale[box][move][kind], 1.0e-10,
                       (boxDimRef.axis.Min(box) * 0.5) - TINY_AMOUNT);
  } else if(move == mv::ROTATE) {
    if(tempTries[box][move][kind] > 0) {
      double currentAccept = (double)(tempAccepted[box][move][kind]) /
                             (double)(tempTries[box][move][kind]);
      double fractOfTargetAccept = currentAccept / TARGET_ACCEPT_FRACT;
      if (fractOfTargetAccept > 0.0) {
        scale[box][move][kind] *= fractOfTargetAccept;
      } else {
        scale[box][move][kind] *= 0.5;
      }
    }
    num::Bound<double>(scale[box][move][kind], 1.0e-6, M_PI - TINY_AMOUNT);
  }
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
  else if(move == mv::VOL_TRANSFER) {
    if(tempTries[box][move][kind] > 0) {
      double currentAccept = (double)(tempAccepted[box][move][kind]) /
                             (double)(tempTries[box][move][kind]);
      double fractOfTargetAccept = currentAccept / TARGET_ACCEPT_FRACT;
      if (fractOfTargetAccept > 0.0) {
        scale[box][move][kind] *= fractOfTargetAccept;
      } else {
        scale[box][move][kind] *= 0.5;
      }
    }
    //Warning: This will lead to have acceptance > 50%
    double maxVolExchange = boxDimRef.volume[box] - boxDimRef.minVol[box];
    num::Bound<double>(scale[box][move][kind], 0.001,  maxVolExchange - 0.001);
  }
#endif
  tempAccepted[box][move][kind] = 0;
  tempTries[box][move][kind] = 0;
}

uint MoveSettings::GetAcceptTot(const uint box, const uint move) const
{
  uint sum = 0;
  for(uint k = 0; k < totKind; k++) {
    sum += accepted[box][move][k];
  }

  if(move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE || move == mv::MULTIPARTICLE_BM
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
      || move == mv::MEMC
#endif
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
      || move == mv::VOL_TRANSFER
#endif
    ) {
    sum /= totKind;
  }
  return sum;
}

double MoveSettings::GetScaleTot(const uint box, const uint move) const
{
  double sum = 0.0;
  for(uint k = 0; k < totKind; k++) {
    sum += scale[box][move][k];
  }
  return sum / (double)(totKind);
}

uint MoveSettings::GetTrialTot(const uint box, const uint move) const
{
  uint sum = 0;
  for(uint k = 0; k < totKind; k++) {
    sum += tries[box][move][k];
  }

  if(move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE || move == mv::MULTIPARTICLE_BM
#if ENSEMBLE == GEMC || ENSEMBLE == GCMC
      || move == mv::MEMC
#endif
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
      || move == mv::VOL_TRANSFER
#endif
    ) {
    sum /= totKind;
  }

  return sum;
}
