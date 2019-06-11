/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CLOCK_H
#define CLOCK_H

//clock() function; CLOCKS_PER_SEC constant
#include <time.h>
#include "BasicTypes.h"             //uint, ulong
#include <iostream> //for cout
#include <cstdio>
#if defined(__linux__) || defined(__APPLE__)
#include <sys/time.h> //for timing
#elif (_WIN32) || (__CYGWIN__)
#include <time.h>
#endif

struct Clock {
  Clock(): stepsPerOut(0), prevStep(0), lastStep(0), lastTime(0.0) {}
  void Init(const ulong steps, const ulong totSt, const ulong startStep)
  {
    stepsPerOut = steps;
    prevStep = startStep;
#if defined(__linux__) || defined(__APPLE__)
    gettimeofday(&tv, &tz);
    strt = (real)tv.tv_sec + (real)tv.tv_usec / 1000000;
#elif (_WIN32) || (__CYGWIN__)
    strt = clock();
#endif
    lastTime = strt;
    lastStep = totSt - 1;
  }
  void CheckTime(const ulong step);
  void SetStart();
  void SetStop();
  real GetTimDiff();
  void CompletionTime(uint &day, uint &hr, uint &min);

private:

#if defined(__linux__) || defined(__APPLE__)
  struct timeval tv;
  struct timezone tz;
  real strt, stop, lastTime;
#elif (_WIN32) || (__CYGWIN__)
  clock_t strt, stop, lastTime;
#endif
  ulong stepsPerOut, prevStep, lastStep;
};

inline void Clock::CheckTime(const ulong step)
{
  ulong stepDelta = step - prevStep;
  real speed = 0.0;
  if (stepDelta == stepsPerOut && step != lastStep) {
#if defined(__linux__) || defined(__APPLE__)
    gettimeofday(&tv, &tz);
    real currTime = (real)tv.tv_sec + (real)tv.tv_usec / 1000000;
    speed = stepDelta / (currTime - lastTime);
#elif (_WIN32) || (__CYGWIN__)
    clock_t currTime = clock();
    speed = stepDelta / (((real)currTime - lastTime) / CLOCKS_PER_SEC);
#endif
    uint day, hr, min;
    prevStep = step;
    lastTime = currTime;
    CompletionTime(day, hr, min);
    printf("Steps/sec: %7.3f, Simulation ends in: %3d d: %3d h: %3d m \n\n",
           speed, day, hr, min);

  } else if (step == lastStep) {
#if defined(__linux__) || defined(__APPLE__)
    gettimeofday(&tv, &tz);
    stop = (real)tv.tv_sec + (real)tv.tv_usec / 1000000;
    std::cout << "Simulation Time (total): " << (stop - strt)
              << " sec." << std::endl;
#elif (_WIN32) || (__CYGWIN__)
    stop = clock();
    std::cout << "Simulation Time (total): "
              << (((real)stop - strt) / CLOCKS_PER_SEC)
              << " sec." << std::endl;
#endif

  }
}

inline void Clock::SetStart()
{
#if defined(__linux__) || defined(__APPLE__)
  gettimeofday(&tv, &tz);
  strt = (real)tv.tv_sec + (real)tv.tv_usec / 1000000;
#elif (_WIN32) || (__CYGWIN__)
  strt = clock();
#endif
}

inline void Clock::SetStop()
{
#if defined(__linux__) || defined(__APPLE__)
  gettimeofday(&tv, &tz);
  stop = (real)tv.tv_sec + (real)tv.tv_usec / 1000000;
#elif (_WIN32) || (__CYGWIN__)
  stop = clock();
#endif
}

inline real Clock::GetTimDiff()
{
#if defined(__linux__) || defined(__APPLE__)
  return (stop - strt);
#elif (_WIN32) || (__CYGWIN__)
  return (real)(stop - strt) / CLOCKS_PER_SEC;
#endif
}

inline void Clock::CompletionTime(uint &day, uint &hr, uint &min)
{
  real speed = 0.0;
#if defined(__linux__) || defined(__APPLE__)
  speed = (real)(prevStep) / (lastTime - strt);
#elif (_WIN32) || (__CYGWIN__)
  speed = (real)(prevStep) / ((real)(lastTime - strt) / CLOCKS_PER_SEC);
#endif
  ulong rem = lastStep - prevStep;
  ulong sec = rem / speed;
  day = sec / 86400;
  sec = sec % 86400;
  hr = sec / 3600;
  sec = sec % 3600;
  min = ceil(sec / 60.0);
}

#endif /*CLOCK_H*/
