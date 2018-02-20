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
  virtual double GetCoeff();
private:
  uint bPick;
  SystemPotential sysPotNew;
  XYZArray& atomForceRef;
  XYZArray& atomForceNew;
  XYZArray& atomTorqueRef;
  XYZArray& atomTorqueNew;
  Coordinates newMolsPos;
  COM newCOMs;
};

inline MultiParticle::MultiParticle(System &sys, StaticVals const &statV) :
  MoveBase(sys, statV), atomForceRef(sys.atomForcesOld),
  atomForceNew(sys.atomForcesNew), atomTorqueRef(sys.atomTorqueOld),
  atomTorqueNew(sys.atomTorqueOld)
{
}

inline uint MultiParticle::Prep(const double subDraw, const double movPerc)
{
  prng.PickBox(bPick, subDraw, movPerc);
  sysPotRef = calcEnRef.BoxInter(sysPotRef, coordCurrRef, comCurrRef, 
                                 boxDimRef, bPick);
}

inline uint MultiParticle::Transform()
{
  // Based on the reference force decided whether to displace or rotate each
  // individual particle.

  // move particles according to force and torque and store them in the new pos
  return 0;
}

inline void MultiParticle::CalcEn() 
{
  // Calculate the new force and energy and we will compare that to the
  // reference values in Accept() function 

  // copy forces to the old one
  atomForceNew.CopyRange(atomForceRef, 0, 0, atomForceNew.Count());
  atomTorqueNew.CopyRange(atomTorqueRef, 0, 0, atomTorqueNew.Count());
  sysPotNew = calcEnRef.BoxInter(sysPotNew, newMolPos, newCOMs, boxDimRef,
                                 bPick);

  return;
}

inline double MultiParticle::GetCoeff() const
{
  // calculate (w_new->old/w_old->new) and return it.
}

inline void MultiParticle::Accept(const uint rejectState, const uint step)
{
  // Here we compare the values of reference and trial and decide whether to 
  // accept or reject the move
  return;
}