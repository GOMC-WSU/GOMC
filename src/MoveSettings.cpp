/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoveSettings.h" //header spec.
#include "BoxDimensions.h" //For axis sizes
#include "StaticVals.h" //For init info.

#include "../lib/NumLib.h" //For bounding functions.
#include "../lib/GeomLib.h"    //For M_PI

const double MoveSettings::TARGET_ACCEPT_FRACT = 0.50;
const double MoveSettings::TINY_AMOUNT = 0.0000001;

void MoveSettings::Init(StaticVals const& statV)
{
   uint baseAdjust = statV.simEventFreq.perAdjust;
   uint maj = 0, subDiv = 0;
   perAdjust = statV.simEventFreq.perAdjust;
   for (uint m = 0; m < mv::COUNT; m++)
   {
      if (m < mv::SCALEABLE)
      {
	 tempTries[m] = 0;
	 tempAccepted[m] = 0;
	 mv::GetMoveMajIndex(maj, subDiv, m);	 
	 //baseAdjust * 
	 /*   ( statV.movePerc[maj] / subDiv / statV.totalPerc );
	 if (perAdjust[m]==0)
	 {
	    perAdjust[m] = statV.simEventFreq.total;
	    }*/
      }
      acceptPercent[m] = 0.0;
      accepted[m] = tries[m] = 0;
   }

#if ENSEMBLE == NVT || ENSEMBLE == GCMC 
   scale[mv::DISPLACE] = boxDimRef.axis.Min(0)/4;
   scale[mv::ROTATE] = M_PI_4;
#elif ENSEMBLE == GEMC
   scale[mv::DISPLACE] = boxDimRef.axis.Min(0)/4;
   scale[mv::DISPLACE+1] = boxDimRef.axis.Min(1)/4;
   scale[mv::ROTATE*BOX_TOTAL] = M_PI_4;
   scale[mv::ROTATE*BOX_TOTAL+1] = M_PI_4;
   scale[mv::VOL_TRANSFER*BOX_TOTAL] = 500;
   scale[mv::VOL_TRANSFER*BOX_TOTAL+1] = 500;
   GEMC_KIND = statV.kindOfGEMC;
#endif

}

//Process results of move we just did in terms of acceptance counters
void MoveSettings::Update(const bool isAccepted, const uint moveIndex,
			  const uint step)
{ 
   bool adjust = ((step + 1) % perAdjust == 0);

   if (moveIndex < mv::SCALEABLE)
   {
      tempTries[moveIndex]++;
   }
   else
   {
      tries[moveIndex]++;
   }
   //Store to appropriate acceptance counter, if accepted.
   if (isAccepted && moveIndex < mv::SCALEABLE)
   {
      tempAccepted[moveIndex]++;
   }
   else if (isAccepted)
   {
      accepted[moveIndex]++;
   }

   //Refresh acceptance percentage appropriately
   if (moveIndex < mv::SCALEABLE)
      acceptPercent[moveIndex] = 
	 (double)(tempAccepted[moveIndex]+accepted[moveIndex]) / 
	 (double)(tempTries[moveIndex]+tries[moveIndex]);
   else
      acceptPercent[moveIndex] = (double)(accepted[moveIndex]) / 
	 (double)(tries[moveIndex]);

   //Check whether we need to adjust this move's scaling.
   if (adjust)
   {
      for (uint m = 0; m < mv::SCALEABLE; ++m)
      {
#if ENSEMBLE == GCMC
         uint majMoveKind = m;
         if (m >= mv::MOL_TRANSFER)
            --majMoveKind;
#else
	 uint majMoveKind = m/BOX_TOTAL;
#endif
	 uint b = m - (majMoveKind*BOXES_WITH_U_NB);
	 Adjust(majMoveKind, m, b);
      }
   }
}

//Adjust responsibly
void MoveSettings::Adjust(const uint majMoveKind,
			  const uint moveIndex, const uint b)
{   
   if (tempTries[moveIndex] > 0)
   {
      double currentAccept =
	 (double)(tempAccepted[moveIndex])/(double)(tempTries[moveIndex]),
	 fractOfTargetAccept = currentAccept/TARGET_ACCEPT_FRACT;

#if 0
      if (majMoveKind == mv::VOL_TRANSFER)
	 std::cout << std::endl
		   << "currentAccept: " << currentAccept << std::endl
		   << "fractOfTargetAccept: " << fractOfTargetAccept << std::endl
		   << "tempAccepted[moveIndex]: " << tempAccepted[moveIndex]
		   << std::endl
		   << "tempTries[moveIndex]: " << tempTries[moveIndex]
		   << std::endl
		   << "scale[moveIndex]: " << scale[moveIndex] << std::endl;
#endif

      if (fractOfTargetAccept > 0.0)
      {
#if 0
	 //ENSEMBLE == GEMC
	 if (majMoveKind == mv::VOL_TRANSFER)
	    fractOfTargetAccept = pow(fractOfTargetAccept, (1.0/3.0));
#endif
         scale[moveIndex] *= fractOfTargetAccept;
      }
#if 0
      else if (fractOfTargetAccept > 0.0)
      {
         if (fractOfTargetAccept > 1)
	    scale[moveIndex] *= 1.05;
	 else
	    scale[moveIndex] *= 0.95;
      }
#endif
      accepted[moveIndex] += tempAccepted[moveIndex];
      tries[moveIndex] += tempTries[moveIndex];
#if 0
      if (majMoveKind == mv::VOL_TRANSFER)
	 std::cout << "scale[moveIndex]: " << scale[moveIndex] << std::endl;
#endif
   }
   tempAccepted[moveIndex] = tempTries[moveIndex] = 0;
   //Bound our values to prevent to big or too small a move.
   switch (majMoveKind)
   {
   case mv::DISPLACE :
      num::Bound<double>(scale[moveIndex], 0.0000000001, 
			 (boxDimRef.axis.Min(b)/2) - TINY_AMOUNT);
      break; 
   case mv::ROTATE :  
      num::Bound<double>(scale[moveIndex], 0.000001, M_PI-TINY_AMOUNT);      
      break;
   default:
      break;
   }
}

