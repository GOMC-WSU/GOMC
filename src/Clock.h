#ifndef CLOCK_H
#define CLOCK_H

//clock() function; CLOCKS_PER_SEC constant
#include <time.h>
#include "BasicTypes.h"             //uint, ulong
#include <iostream> //for cout
#include <sys/time.h> //for timing

struct Clock
{
  Clock(): stepsPerOut(0), prevStep(0), lastStep(0), lastTime(0.0) {}
  void Init(const ulong steps, const ulong totSt)
  {
    stepsPerOut = steps;
    gettimeofday(&tv, &tz);
    strt = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
    lastTime = strt;
    lastStep = totSt - 1;
  }
  void CheckTime(const uint step);
private:
  double TimeInSec(const double strt, const double stp)
  {
    return (double(stp)-double(strt))/CLOCKS_PER_SEC;
  }
 
  struct timeval tv;
  struct timezone tz;
  double strt, stop, lastTime;
  ulong stepsPerOut, prevStep, lastStep;
};

inline void Clock::CheckTime(const uint step)
{
  uint stepDelta = step - prevStep;
  if (stepDelta == stepsPerOut && step != lastStep)
  {
    gettimeofday(&tv, &tz);
    double currTime = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
    std::cout << "Steps/sec. : "
              << stepDelta/(currTime - lastTime) << std::endl;
    prevStep = step;
    lastTime = currTime;
  }
  else if (step == lastStep)
  {
    gettimeofday(&tv, &tz);
    stop = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
    std::cout << "Simulation Time (total): " << (stop - strt)
              << "sec." << std::endl;
  }
}

#endif /*CLOCK_H*/
