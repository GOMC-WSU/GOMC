/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.1
Copyright (C) 2016  GOMC Group
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
                        pdb_setup::Remarks const& remarks)
{
  perAdjust = statV.simEventFreq.perAdjust;
  for (uint m = 0; m < mv::COUNT; m++) {
    if (m < mv::SCALEABLE) {
      tempTries[m] = 0;
      tempAccepted[m] = 0;
    }
    acceptPercent[m] = 0.0;
    accepted[m] = tries[m] = 0;
  }

  if(remarks.restart) {
    for(uint b = 0; b < BOX_TOTAL; b++) {
      uint disp = mv::GetMoveSubIndex(mv::DISPLACE, b);
      scale[disp] = remarks.disp[b];
      uint rotate = mv::GetMoveSubIndex(mv::ROTATE, b);
      scale[rotate] = remarks.rotate[b];
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
      uint volume = mv::GetMoveSubIndex(mv::VOL_TRANSFER, b);
      scale[volume] = remarks.vol[b];
#endif
    }
  } else {
    for(uint b = 0; b < BOX_TOTAL; b++) {
      uint disp = mv::GetMoveSubIndex(mv::DISPLACE, b);
      scale[disp] = boxDimRef.axis.Min(b) / 4;;
      uint rotate = mv::GetMoveSubIndex(mv::ROTATE, b);
      scale[rotate] = M_PI_4;
#if ENSEMBLE == NPT || ENSEMBLE == GEMC
      uint volume = mv::GetMoveSubIndex(mv::VOL_TRANSFER, b);
      scale[volume] = 500;
#endif
    }
  }

}

//Process results of move we just did in terms of acceptance counters
void MoveSettings::Update(const bool isAccepted, const uint moveIndex,
                          const uint step)
{
  if (moveIndex < mv::SCALEABLE) {
    tempTries[moveIndex]++;
  } else {
    tries[moveIndex]++;
  }

  //Store to appropriate acceptance counter, if accepted.
  if (isAccepted && moveIndex < mv::SCALEABLE) {
    tempAccepted[moveIndex]++;
  } else if (isAccepted) {
    accepted[moveIndex]++;
  }

  //Refresh acceptance percentage appropriately
  if (moveIndex < mv::SCALEABLE)
    acceptPercent[moveIndex] =
      (double)(tempAccepted[moveIndex] + accepted[moveIndex]) /
      (double)(tempTries[moveIndex] + tries[moveIndex]);
  else
    acceptPercent[moveIndex] = (double)(accepted[moveIndex]) /
                               (double)(tries[moveIndex]);

}

void MoveSettings::AdjustMoves(const uint step)
{
  //Check whether we need to adjust this move's scaling.
  if ((step + 1) % perAdjust == 0) {
    for (uint m = 0; m < mv::SCALEABLE; ++m) {
      uint majMoveKind = m / BOX_TOTAL;

      uint b = m - (majMoveKind * BOXES_WITH_U_NB);
      Adjust(majMoveKind, m, b);
    }
  }
}

//Adjust responsibly
void MoveSettings::Adjust(const uint majMoveKind,
                          const uint moveIndex, const uint b)
{
  if (tempTries[moveIndex] > 0) {
    double currentAccept = (double)(tempAccepted[moveIndex]) /
                           (double)(tempTries[moveIndex]);

    double fractOfTargetAccept = currentAccept / TARGET_ACCEPT_FRACT;


    if (fractOfTargetAccept > 0.0) {
      scale[moveIndex] *= fractOfTargetAccept;
    }

    accepted[moveIndex] += tempAccepted[moveIndex];
    tries[moveIndex] += tempTries[moveIndex];
#if 0
    if (majMoveKind == mv::VOL_TRANSFER)
      std::cout << "scale[moveIndex]: " << scale[moveIndex] << std::endl;
#endif
  }
  tempAccepted[moveIndex] = tempTries[moveIndex] = 0;
  //Bound our values to prevent to big or too small a move.
  switch (majMoveKind) {
  case mv::DISPLACE :
    num::Bound<double>(scale[moveIndex], 0.0000000001,
                       (boxDimRef.axis.Min(b) / 2) - TINY_AMOUNT);
    break;
  case mv::ROTATE :
    num::Bound<double>(scale[moveIndex], 0.000001, M_PI - TINY_AMOUNT);
    break;
  default:
    break;
  }
}
