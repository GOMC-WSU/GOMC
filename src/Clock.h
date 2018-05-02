/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.20
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
#if defined(__linux__) || defined(__APPLE__)
#include <sys/time.h> //for timing
#elif _WIN32
#include <time.h>
#endif

struct Clock {
  Clock(): stepsPerOut(0), prevStep(0), lastStep(0), lastTime(0.0) {}
  void Init(const ulong steps, const ulong totSt)
  {
    stepsPerOut = steps;
#if defined(__linux__) || defined(__APPLE__)
    gettimeofday(&tv, &tz);
    strt = (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
#elif _WIN32
    strt = clock();
#endif
    lastTime = strt;
    lastStep = totSt - 1;
  }
  void CheckTime(const uint step);
  void SetStart();
  void SetStop();
  double GetTimDiff();

private:
  double TimeInSec(const double strt, const double stp)
  {
    return (double(stp) - double(strt)) / CLOCKS_PER_SEC;
  }

#if defined(__linux__) || defined(__APPLE__)
  struct timeval tv;
  struct timezone tz;
  double strt, stop, lastTime;
#elif _WIN32
  clock_t strt, stop, lastTime;
#endif
  ulong stepsPerOut, prevStep, lastStep;
};

inline void Clock::CheckTime(const uint step)
{
  uint stepDelta = step - prevStep;
  if (stepDelta == stepsPerOut && step != lastStep) {
#if defined(__linux__) || defined(__APPLE__)
    gettimeofday(&tv, &tz);
    double currTime = (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
    std::cout << "Steps/sec. : "
              << stepDelta / (currTime - lastTime) << std::endl;
#elif _WIN32
    clock_t currTime = clock();
    std::cout << "Steps/sec. : "
              << stepDelta / ((double)currTime - lastTime / CLOCKS_PER_SEC)
              << std::endl;
#endif
    prevStep = step;
    lastTime = currTime;
  } else if (step == lastStep) {
#if defined(__linux__) || defined(__APPLE__)
    gettimeofday(&tv, &tz);
    stop = (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
    std::cout << "Simulation Time (total): " << (stop - strt)
              << " sec." << std::endl;
#elif _WIN32
    stop = clock();
    std::cout << "Simulation Time (total): "
              << ((double)stop - strt / CLOCKS_PER_SEC)
              << " sec." << std::endl;
#endif

  }
}

inline void Clock::SetStart()
{
#if defined(__linux__) || defined(__APPLE__)
  gettimeofday(&tv, &tz);
  strt = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
#elif _WIN32
  strt = clock();
#endif
}

inline void Clock::SetStop()
{
#if defined(__linux__) || defined(__APPLE__)
  gettimeofday(&tv, &tz);
  stop = (double)tv.tv_sec + (double)tv.tv_usec/1000000;
#elif _WIN32
  stop = clock();
#endif
}

inline double Clock::GetTimDiff()
{
  return (stop - strt);
}

#endif /*CLOCK_H*/
