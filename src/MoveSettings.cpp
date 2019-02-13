/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.31
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoveSettings.h" //header spec.
#include "BoxDimensions.h" //For axis sizes
#include "BoxDimensionsNonOrth.h"
#include "StaticVals.h" //For init info.

#include "NumLib.h" //For bounding functions.
#include "GeomLib.h"    //For M_PI

const double MoveSettings::TARGET_ACCEPT_FRACT = 0.50;
const double MoveSettings::TINY_AMOUNT = 0.0000001;

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
      for(uint k = 0; k < totKind; k++){
        if(m == mv::DISPLACE) {
          if(remarks.restart && remarks.disp[b] > 0.0) {
            scale[b][m][k] = remarks.disp[b];
          } else {
            scale[b][m][k] = boxDimRef.axis.Min(b) / 4;
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
    
    if(move != mv::MULTIPARTICLE)
      isSingleMoveAccepted = true;
  }

  if(move == mv::MULTIPARTICLE)
    isSingleMoveAccepted = false;

  acceptPercent[box][move][kind] = (double)(accepted[box][move][kind]) /
                                  (double)(tries[box][move][kind]);

  //for any move that we dont care about kind of molecule, it should be included
  //in the if condition  
  if (move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE
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
      acceptPercent[box][move][kind] = (double)(accepted[box][move][kind]) /
                                      (double)(tries[box][move][kind]);
    }
  }
}

void MoveSettings::AdjustMoves(const uint step)
{
  //Check whether we need to adjust this move's scaling.
  if ((step + 1) % perAdjust == 0) {
    for(uint b = 0; b < BOX_TOTAL; b++) {
      for (uint m = 0; m < mv::MOVE_KINDS_TOTAL; ++m) {
        for(uint k = 0; k < totKind; k++) {
          Adjust(b, m, k);
        }
      }
    }
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
    num::Bound<double>(scale[box][move][kind], 0.0000000001,
                       (boxDimRef.axis.Min(box) / 2) - TINY_AMOUNT);
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
    num::Bound<double>(scale[box][move][kind], 0.000001, M_PI - TINY_AMOUNT);
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
    //Warning: This will lead to have acceptance > %50
    double maxVolExchange = boxDimRef.volume[box] - boxDimRef.minVol[box];
    num::Bound<double>(scale[box][move][kind], 0.001,  maxVolExchange - 0.001);
  }
#endif
  tempAccepted[box][move][kind] = 0;
  tempTries[box][move][kind] = 0;
}

uint MoveSettings::GetAcceptTot(const uint box, const uint move) const
{
  uint sum= 0;
  for(uint k = 0; k < totKind; k++) {
    sum += accepted[box][move][k];
  }

  if(move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE
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
  uint sum = 0.0;
  for(uint k = 0; k < totKind; k++) {
    sum += tries[box][move][k];
  }

  if(move == mv::INTRA_MEMC || move == mv::MULTIPARTICLE
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
