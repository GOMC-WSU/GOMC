#include "MoveSettings.h" //header spec.
#include "BoxDimensions.h" //For axis sizes
#include "StaticVals.h" //For init info.

#include "NumLib.h" //For bounding functions.
#include "GeomLib.h"    //For M_PI

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
      }
      acceptPercent[m] = 0.0;
      accepted[m] = tries[m] = 0;
   }

#if ENSEMBLE == NVT || ENSEMBLE == GCMC
   scale[mv::DISPLACE] = boxDimRef.axis.Min(0)/4;
   scale[mv::ROTATE] = M_PI_4;
#elif ENSEMBLE == NPT
   scale[mv::DISPLACE] = boxDimRef.axis.Min(0)/4;
   scale[mv::ROTATE] = M_PI_4;
   scale[mv::VOL_TRANSFER] = 500;
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
 
}

void MoveSettings::AdjustMoves(const uint step)
{
   //Check whether we need to adjust this move's scaling.
   if ((step + 1) % perAdjust == 0)
   {
      for (uint m = 0; m < mv::SCALEABLE; ++m)
      {
	 uint majMoveKind = m/BOX_TOTAL;

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
      double currentAccept = (double)(tempAccepted[moveIndex])/
	(double)(tempTries[moveIndex]);

      double fractOfTargetAccept = currentAccept/TARGET_ACCEPT_FRACT;


      if (fractOfTargetAccept > 0.0)
      {
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
