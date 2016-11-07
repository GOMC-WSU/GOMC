#ifndef MOVE_SETTINGS_H
#define MOVE_SETTINGS_H

#include "EnsemblePreprocessor.h" //For BOX_TOTAL
#include "BasicTypes.h"           //for uint
#include "OutputVars.h"

#include "MoveConst.h"           //For sizes of arrays.

class StaticVals;                 //For various initialization constants.
class BoxDimensions;              //For axis sizes

class MoveSettings
{
public:
  friend class OutputVars;
  MoveSettings(BoxDimensions & dim) : boxDimRef(dim) {}

  MoveSettings& operator=(MoveSettings const& rhs)
  {
    return *this;
  }

  void Init(StaticVals const& statV);

  void Update(const bool isAccepted, const uint moveIndex, const uint step);

  void Adjust(const uint majMoveKind, const uint moveIndex, const uint b);

  double Scale(const uint move) const
  {
    return scale[move];
  }

private:
  double scale[mv::SCALEABLE];
  double acceptPercent[mv::COUNT];
  uint accepted[mv::COUNT];
  uint tries[mv::COUNT];
  uint perAdjust;
  uint tempAccepted[mv::SCALEABLE], tempTries[mv::SCALEABLE];

#if ENSEMBLE == GEMC
  uint GEMC_KIND;
#endif

  BoxDimensions & boxDimRef;

  static const double TARGET_ACCEPT_FRACT;
  static const double TINY_AMOUNT;
};


#endif /*MOVE_SETTINGS_H*/
