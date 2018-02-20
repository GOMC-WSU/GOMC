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
};

inline MultiParticle::MultiParticle(System &sys, StaticVals const &statV) :
  MoveBase(sys, statV)
{
}

inline uint MultiParticle::Prep(const double subDraw, const double movPerc)
{
  
}