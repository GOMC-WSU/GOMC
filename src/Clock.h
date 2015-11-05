/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CLOCK_H
#define CLOCK_H

//clock() function; CLOCKS_PER_SEC constant
#include <time.h>
#include "../lib/BasicTypes.h"             //uint, ulong
#include <iostream> //for cout

struct Clock 
{
   Clock(): stepsPerOut(0), prevStep(0), lastStep(0), lastTime(0.0) {}
   void Init(const ulong steps, const ulong totSt)
   { stepsPerOut = steps; strt = clock(); lastStep = totSt - 1; }
   void CheckTime(const uint step);
private:
   double TimeInSec(const double strt, const double stp)
   { return (double(stp)-double(strt))/CLOCKS_PER_SEC; }
   clock_t strt, stop;
   double lastTime;
   ulong stepsPerOut, prevStep, lastStep;
};

inline void Clock::CheckTime(const uint step)
{
   uint stepDelta = step - prevStep;
   if (stepDelta == stepsPerOut && step != lastStep)
   {
      double currTime = clock();
      std::cout << "Steps/sec. : " 
		<< stepDelta/TimeInSec(lastTime, currTime) << std::endl;
      prevStep = step;
      lastTime = currTime;
   }
   else if (step == lastStep)
   {
      stop = clock();
      std::cout << "Simulation Time (total): " << TimeInSec(strt, stop)
		<< "sec." << std::endl;
   }
}

#endif /*CLOCK_H*/

