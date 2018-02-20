#pragma once

#include "MoveBase.h"
#include "System.h"
#include "StaticVals.h"

class MultiParticle : public MoveBase
{
public:
  MultiParticle(System &sys, StaticVals const& statV);

  virtual uint Prep(const double subDraw, const double movPerc);
  virtual void CalcEn();
  virtual uint Transform();
  virtual void Accept(const uint rejectState, const uint step);
private:
  uint bPick;
  SystemPotential sysPotNew;
};

inline MultiParticle::MultiParticle(System &sys, StaticVals const &statV) :
  MoveBase(sys, statV)
{
}

inline uint MultiParticle::Prep(const double subDraw, const double movPerc)
{
  prng.PickBox(bPick, subDraw, movPerc);
  sysPotRef = calcEnRef.BoxInter(sysPotRef, coordCurrRef, comCurrRef, 
                                 boxDimRef, bPick);
}

inline uint Transform()
{
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.

}

inline void MultiParticle::CalcEn() 
{
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function 
}

inline void Accept(const uint rejectState, const uint step)
{
  // Here we compare the values of reference and trial and decide whether to 
  // accept or reject the move
}