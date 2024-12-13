/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "MoveSettings.h"  //header spec
#include "BoxDimensions.h" //For axis sizes
#include "BoxDimensionsNonOrth.h"
#include "GeomLib.h"    //For M_PI
#include "NumLib.h"     //For bounding functions
#include "StaticVals.h" //For init info

const double MoveSettings::TARGET_ACCEPT_FRACT = 0.5;
const double MoveSettings::TINY_AMOUNT = 1.0e-7;
const double MoveSettings::r_alpha = 0.2;
const double MoveSettings::t_alpha = 0.2;
const double MoveSettings::mp_accept_tol = 0.1;

void MoveSettings::Init(StaticVals const &statV,
                        pdb_setup::Remarks const &remarks, const uint tkind,
                        bool restartFromCheckpoint) {
  // Set to true so that we calculate the forces for the current system, even if
  // a MultiParticle move is called before any other moves are accepted.

  totKind = tkind;
  perAdjust = statV.simEventFreq.perAdjust;
  if (!restartFromCheckpoint) {
    for (uint b = 0; b < BOXES_WITH_U_NB; b++) {
      SetSingleMoveAccepted(b);
    }
    for (uint b = 0; b < BOX_TOTAL; b++) {
      for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++) {
        acceptPercent[b][m].resize(totKind, 0);
        scale[b][m].resize(totKind, 0);
        accepted[b][m].resize(totKind, 0);
        tries[b][m].resize(totKind, 0);
        tempAccepted[b][m].resize(totKind, 0);
        tempTries[b][m].resize(totKind, 0);
      }
    }

    for (uint b = 0; b < BOX_TOTAL; b++) {
      for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; m++) {
        for (uint k = 0; k < totKind; k++) {
          if (m == mv::DISPLACE) {
            if (remarks.restart && remarks.disp[b] > 0.0) {
              scale[b][m][k] = remarks.disp[b];
            } else {
              scale[b][m][k] = boxDimRef.axis.Min(b) * 0.25;
            }
          } else if (m == mv::ROTATE) {
            if (remarks.restart && remarks.rotate[b] > 0.0) {
              scale[b][m][k] = remarks.rotate[b];
            } else {
              scale[b][m][k] = M_PI_4;
            }
          }
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
          else if (m == mv::VOL_TRANSFER) {
            if (remarks.restart && remarks.vol[b] > 0.0) {
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
    for (int b = 0; b < BOX_TOTAL; b++) {
      // Start with smaller maximum adjustments for Brownian Motion MP
      // Also need to check for NeMTMC move using Brownian Motion Multiparticle.
      if (statV.movePerc[mv::MULTIPARTICLE_BM] > TINY_AMOUNT ||
          (statV.neMTMCVal.enable && statV.neMTMCVal.MPBEnable)) {
        mp_r_max[b] = 0.0005 * M_PI;
        mp_t_max[b] = 0.002;
      } else {
        mp_r_max[b] = 0.01 * M_PI;
        mp_t_max[b] = 0.02;
      }
      for (int m = 0; m < mp::MPTOTALTYPES; m++) {
        mp_tries[b][m] = 0;
        mp_accepted[b][m] = 0;
        mp_interval_tries[b][m] = 0;
        mp_interval_accepted[b][m] = 0;
      }
    }
  }
}

// Process results of move we just did in terms of acceptance counters
void MoveSettings::Update(const uint move, const bool isAccepted,
                          const uint box, const uint kind) {
  tries[box][move][kind]++;
  tempTries[box][move][kind]++;
  if (isAccepted) {
    tempAccepted[box][move][kind]++;
    accepted[box][move][kind]++;

    if (move != mv::MULTIPARTICLE || move != mv::MULTIPARTICLE_BM) {
      SetSingleMoveAccepted(box);
#if ENSEMBLE == GEMC
      // GEMC has multiple boxes and this move changed both boxes, so we
      // need to also mark the other box as having a move accepted.
      if (move == mv::MEMC || move == mv::MOL_TRANSFER || move == mv::NE_MTMC ||
          move == mv::VOL_TRANSFER || move == mv::TARGETED_SWAP) {
        // Simple way to figure out which box is the other one. 0 -->1 and 1-->0
        // Assumes just two boxes.
        uint otherBox = box == 0;
        SetSingleMoveAccepted(otherBox);
      }
#endif
    }
  }

  if (move == mv::MULTIPARTICLE || move == mv::MULTIPARTICLE_BM)
    UnsetSingleMoveAccepted(box);

  acceptPercent[box][move][kind] =
      (double)(accepted[box][move][kind]) / (double)(tries[box][move][kind]);

  // for any move that we don't care about kind of molecule, it should be
  // included in the if condition
  if (move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE ||
      move == mv::MULTIPARTICLE_BM
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
      if (isAccepted) {
        tempAccepted[box][move][k]++;
        accepted[box][move][k]++;
      }
      acceptPercent[box][move][k] =
          (double)(accepted[box][move][k]) / (double)(tries[box][move][k]);
    }
  }
}

void MoveSettings::AdjustMoves(const ulong step) {
  // Check whether we need to adjust this move's scaling.
  if ((step + 1) % perAdjust == 0) {
    for (uint b = 0; b < BOX_TOTAL; b++) {
      for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; ++m) {
        for (uint k = 0; k < totKind; k++) {
          Adjust(b, m, k);
        }
      }
      for (int m = 0; m < mp::MPTOTALTYPES; m++) {
        AdjustMultiParticle(b, m);
      }
    }
  }
}

void MoveSettings::AdjustMultiParticle(const uint box, const uint typePick) {
  // Make sure we tried some moves of this move type, otherwise move max will be
  // NaN
  if (mp_interval_tries[box][typePick] > 0) {
    // Update totals for the entire simulation
    mp_tries[box][typePick] += mp_interval_tries[box][typePick];
    mp_accepted[box][typePick] += mp_interval_accepted[box][typePick];
    double fractOfTotalAccept =
        ((double)mp_accepted[box][typePick] / (double)mp_tries[box][typePick]) /
        mp::TARGET_ACCEPT_FRACT;
    double fractOfIntervalAccept =
        ((double)mp_interval_accepted[box][typePick] /
         (double)mp_interval_tries[box][typePick]) /
        mp::TARGET_ACCEPT_FRACT;
    if (typePick == mp::MPDISPLACE) {
      if (fractOfIntervalAccept == 0.0) {
        mp_t_max[box] *= 0.5;
      } else if (std::fabs(fractOfIntervalAccept - mp::TARGET_ACCEPT_FRACT) >
                 mp_accept_tol) {
        mp_t_max[box] *= ((1.0 - t_alpha) * fractOfTotalAccept +
                          t_alpha * fractOfIntervalAccept);
      }
      num::Bound<double>(mp_t_max[box], 0.001,
                         (boxDimRef.axis.Min(box) * 0.5) - 0.001);
    } else {
      if (fractOfIntervalAccept == 0.0) {
        mp_r_max[box] *= 0.5;
      } else if (std::fabs(fractOfIntervalAccept - mp::TARGET_ACCEPT_FRACT) >
                 mp_accept_tol) {
        mp_r_max[box] *= ((1.0 - r_alpha) * fractOfTotalAccept +
                          r_alpha * fractOfIntervalAccept);
      }
      num::Bound<double>(mp_r_max[box], 0.001, M_PI - 0.001);
    }
    // Reset totals for next adjustment period
    mp_interval_accepted[box][typePick] = 0;
    mp_interval_tries[box][typePick] = 0;
  }
}

void MoveSettings::UpdateMoveSettingMultiParticle(const uint box, bool isAccept,
                                                  const uint typePick) {
  mp_interval_tries[box][typePick]++;
  if (isAccept) {
    mp_interval_accepted[box][typePick]++;
  }
}

// Adjust responsibly
void MoveSettings::Adjust(const uint box, const uint move, const uint kind) {
  if (move == mv::DISPLACE) {
    if (tempTries[box][move][kind] > 0) {
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
  } else if (move == mv::ROTATE) {
    if (tempTries[box][move][kind] > 0) {
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
  else if (move == mv::VOL_TRANSFER) {
    if (tempTries[box][move][kind] > 0) {
      double currentAccept = (double)(tempAccepted[box][move][kind]) /
                             (double)(tempTries[box][move][kind]);
      double fractOfTargetAccept = currentAccept / TARGET_ACCEPT_FRACT;
      if (fractOfTargetAccept > 0.0) {
        scale[box][move][kind] *= fractOfTargetAccept;
      } else {
        scale[box][move][kind] *= 0.5;
      }
    }
    // Warning: This will lead to have acceptance > 50%
    double maxVolExchange = boxDimRef.volume[box] - boxDimRef.minVol[box];
    num::Bound<double>(scale[box][move][kind], 0.001, maxVolExchange - 0.001);
  }
#endif
  tempAccepted[box][move][kind] = 0;
  tempTries[box][move][kind] = 0;
}

uint MoveSettings::GetAcceptTot(const uint box, const uint move) const {
  uint sum = 0;
  for (uint k = 0; k < totKind; k++) {
    sum += accepted[box][move][k];
  }

  if (move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE ||
      move == mv::MULTIPARTICLE_BM
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

double MoveSettings::GetScaleTot(const uint box, const uint move) const {
  double sum = 0.0;
  for (uint k = 0; k < totKind; k++) {
    sum += scale[box][move][k];
  }
  return sum / (double)(totKind);
}

uint MoveSettings::GetTrialTot(const uint box, const uint move) const {
  uint sum = 0;
  for (uint k = 0; k < totKind; k++) {
    sum += tries[box][move][k];
  }

  if (move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE ||
      move == mv::MULTIPARTICLE_BM
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

#if GOMC_GTEST

bool MoveSettings::operator==(const MoveSettings &rhs) {
  bool result = true;

  result &= (scale == rhs.scale);
  result &= (acceptPercent == rhs.acceptPercent);
  // index [BOX_TOTAL * kind + box] is the first element of that kind/box in
  // molLookup
  // index [BOX_TOTAL * kind + box + 1] is the element after the end
  // of that kind/box
  result &= (accepted == rhs.accepted);
  result &= (tries == rhs.tries);
  result &= (tempAccepted == rhs.tempAccepted);
  result &= (tempTries == rhs.tempTries);
  result &= (mp_accepted == rhs.mp_accepted);
  result &=
      (mp_tries == rhs.mp_tries); // Kinds that can move intra and inter box
  result &= (mp_interval_accepted ==
             rhs.mp_interval_accepted); // Kinds that can move intra box only
  result &=
      (mp_interval_tries == rhs.mp_interval_tries); // stores the molecule index
                                                    // for global atom index
  result &= (mp_r_max ==
             rhs.mp_r_max); // stores the local atom index for global atom index
  result &= (mp_t_max ==
             rhs.mp_t_max); // stores the molecule kind for global atom index
  result &= (isSingleMoveAccepted == rhs.isSingleMoveAccepted);

  return result;
}
#endif