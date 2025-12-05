/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "TrialMol.h"

#include <algorithm>
#include <utility> //swap

#include "BoxDimensions.h"
#include "GeomLib.h" //for Theta
#ifndef NDEBUG
#include <iostream>
#endif

namespace cbmc {

TrialMol::TrialMol(const MoleculeKind &k, const BoxDimensions &ax, uint box)
    : kind(&k), axes(&ax), box(box), tCoords(k.NumAtoms()), cavMatrix(3),
      bCoords(k.NumAtoms()), totalWeight(1.0), bonds(k.bondList) {
  cavMatrix.Set(0, 1.0, 0.0, 0.0);
  cavMatrix.Set(1, 0.0, 1.0, 0.0);
  cavMatrix.Set(2, 0.0, 0.0, 1.0);
  growthToWorld.LoadIdentity();
  backbone[0] = backbone[1] = 0;
  growingAtomIndex = 0;
  comInCav = false;
  comFix = false;
  rotateBB = false;
  overlap = false;
  atomBuilt = new bool[k.NumAtoms()];
  std::fill_n(atomBuilt, k.NumAtoms(), false);
}

TrialMol::TrialMol()
    : kind(NULL), axes(NULL), box(0), tCoords(0), cavMatrix(3), bCoords(0),
      comInCav(false), comFix(false), rotateBB(false), overlap(false),
      atomBuilt(NULL), bonds() {
  backbone[0] = backbone[1] = 0;
  growingAtomIndex = 0;
  cavMatrix.Set(0, 1.0, 0.0, 0.0);
  cavMatrix.Set(1, 0.0, 1.0, 0.0);
  cavMatrix.Set(2, 0.0, 0.0, 1.0);
}

TrialMol::TrialMol(const TrialMol &other)
    : kind(other.kind), axes(other.axes), box(other.box),
      tCoords(other.tCoords), cavMatrix(other.cavMatrix),
      bCoords(other.bCoords), en(other.en), totalWeight(other.totalWeight),
      basisPoint(other.basisPoint), bonds(other.bonds) {
  atomBuilt = new bool[kind->NumAtoms()];
  std::copy(other.atomBuilt, other.atomBuilt + kind->NumAtoms(), atomBuilt);
  bonds.Unset();
  backbone[0] = backbone[1] = 0;
  growingAtomIndex = 0;
  comInCav = false;
  comFix = false;
  rotateBB = false;
  overlap = false;
  cavMatrix.Set(0, 1.0, 0.0, 0.0);
  cavMatrix.Set(1, 0.0, 1.0, 0.0);
  cavMatrix.Set(2, 0.0, 0.0, 1.0);
}

TrialMol &TrialMol::operator=(TrialMol other) {
  swap(*this, other);
  return *this;
}

void swap(TrialMol &a, TrialMol &b) {
  using std::swap;
  swap(a.kind, b.kind);
  swap(a.axes, b.axes);
  swap(a.box, b.box);
  swap(a.tCoords, b.tCoords);
  swap(a.bCoords, b.bCoords);
  swap(a.en, b.en);
  swap(a.totalWeight, b.totalWeight);
  swap(a.atomBuilt, b.atomBuilt);
  swap(a.growthToWorld, b.growthToWorld);
  swap(a.worldToGrowth, b.worldToGrowth);
  swap(a.bonds, b.bonds);
  a.bonds.Unset();
  b.bonds.Unset();
  a.comInCav = false;
  b.comInCav = false;
  a.comFix = false;
  b.comFix = false;
  a.rotateBB = false;
  b.rotateBB = false;
  a.overlap = false;
  b.overlap = false;
  a.cavMatrix.Set(0, 1.0, 0.0, 0.0);
  a.cavMatrix.Set(1, 0.0, 1.0, 0.0);
  a.cavMatrix.Set(2, 0.0, 0.0, 1.0);
  b.cavMatrix.Set(0, 1.0, 0.0, 0.0);
  b.cavMatrix.Set(1, 0.0, 1.0, 0.0);
  b.cavMatrix.Set(2, 0.0, 0.0, 1.0);
}

TrialMol::~TrialMol() { delete[] atomBuilt; }

void TrialMol::AddAtom(const uint index, const XYZ &loc) {
  tCoords.Set(index, loc);
  atomBuilt[index] = true;
}

void TrialMol::SetAtomCoords(uint index, const XYZ &loc) {
  tCoords.Set(index, loc);
}

void TrialMol::ConfirmOldAtom(uint i) { atomBuilt[i] = true; }

//! Returns rectangular coordinates of an addition
//! Determines coordinates with respect to current basis
XYZ TrialMol::GetRectCoords(double bond, double theta, double phi) const {
  // rotation/translation from growth to sim coordinates
  return axes->WrapPBC(RawRectCoords(bond, theta, phi) + basisPoint, box);
}

// conversion from polar coordinates to rectangular growth coordinates
XYZ TrialMol::RawRectCoords(double bond, double theta, double phi) const {
  XYZ result(bond, bond, bond);
  double sinTh = sin(theta);
  result.x *= sinTh * cos(phi);
  result.y *= sinTh * sin(phi);
  result.z *= cos(theta);
  return growthToWorld.Apply(result);
}

void TrialMol::OldThetaAndPhi(const uint atom, const uint lastAtom,
                              double &theta, double &phi) const {
  XYZ diff = tCoords.Difference(atom, lastAtom);
  diff = axes->MinImage(diff, box);
  XYZ growthCoords = worldToGrowth.Apply(diff);
  theta = acos(growthCoords.z / growthCoords.Length());
  phi = atan2(growthCoords.y, growthCoords.x);
  return;
}

double TrialMol::OldDistSq(const uint lastAtom, const uint atom) {
  XYZ diff = tCoords.Difference(atom, lastAtom);
  diff = axes->MinImage(diff, box);
  double distSq = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
  return distSq;
}

double TrialMol::DistSq(const XYZ &a, const XYZ &b) {
  XYZ diff = a - b;
  diff = axes->MinImage(diff, box);
  double distSq = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
  return distSq;
}

//! Return angle in radians between confirmed atoms a-b-c
double TrialMol::GetTheta(uint a, uint b, uint c) const {
  return geom::Theta(axes->MinImage(tCoords.Difference(a, b), box),
                     axes->MinImage(tCoords.Difference(c, b), box));
}

//! Return dihedral in radians between confirmed atoms a-b-c-d
double TrialMol::GetPhi(uint a, uint b, uint c, uint d) const {
  return geom::Phi(axes->MinImage(tCoords.Difference(b, a), box),
                   axes->MinImage(tCoords.Difference(c, b), box),
                   axes->MinImage(tCoords.Difference(d, c), box));
}

void TrialMol::SetBasis(const uint p1, const uint p2, const uint p3) {
  using namespace geom;
  // W is unit vec of p1->p2
  XYZ wVec = axes->MinImage(tCoords.Difference(p2, p1), box);
  wVec.Normalize();
  // U will be unit projection of p2->p3 onto plane normal to W
  XYZ uVec = axes->MinImage(tCoords.Difference(p3, p2), box);
  // V is unit vec perpendicular to both W and U
  XYZ vVec = Cross(wVec, uVec);
  vVec.Normalize();
  // Finish X'
  uVec = Cross(vVec, wVec);

  growthToWorld.BasisRotation(uVec, vVec, wVec);
  worldToGrowth = growthToWorld.Inverse();
  basisPoint = tCoords.Get(p1);
}

void TrialMol::SetBasis(const uint p1, const uint p2) {
  using namespace geom;
  // W is unit vec of p1->p2
  XYZ wVec = axes->MinImage(tCoords.Difference(p2, p1), box);
  wVec.Normalize();
  XYZ uVec;
  // check to make sure our W isn't in line with the standard X Axis
  if (std::abs(wVec.x) < 0.8) {
    // V will be W x the standard X unit vec
    uVec = XYZ(1.0, 0.0, 0.0);
  } else {
    // V will be W x the standard Y unit vec
    uVec = XYZ(0.0, 1.0, 0.0);
  }
  XYZ vVec = Cross(wVec, uVec);
  vVec.Normalize();
  // U is unit vec perpendicular to both V and W
  uVec = Cross(vVec, wVec);
  growthToWorld.BasisRotation(uVec, vVec, wVec);
  worldToGrowth = growthToWorld.Inverse();
  basisPoint = tCoords.Get(p1);
}

void TrialMol::ShiftBasis(const uint p1) { basisPoint = tCoords.Get(p1); }

void TrialMol::ShiftBasis(const XYZ cent) { basisPoint = cent; }

void TrialMol::ResetBasis() {
  growthToWorld.LoadIdentity();
  worldToGrowth.LoadIdentity();
  basisPoint = XYZ(0, 0, 0);
}

double TrialMol::PhiBetweenAngles(double theta1, double theta2,
                                  double interior) {
  double y =
      (cos(interior) - cos(theta1) * cos(theta2)) / (sin(theta1) * sin(theta2));
  return M_PI_2 - atan2(y, sqrt(1 - y * y));
}

void TrialMol::SetCoords(const XYZArray &coords, uint start) {
  coords.CopyRange(tCoords, start, 0, tCoords.Count());
}

void TrialMol::SetBCoords(const XYZArray &coords, uint start) {
  coords.CopyRange(bCoords, start, 0, bCoords.Count());
}

double TrialMol::AngleDist(const double b1, const double b2,
                           const double theta) {
  if (!kind->oneThree)
    return 0.0;
  else {
    double v = b1 * b1 - 2 * b1 * b2 * cos(theta) + b2 * b2;
    return v;
  }
}

double TrialMol::DihedDist(const double b1, const double b2, const double b3,
                           const double theta1, const double theta2,
                           const double phi) {
  if (!kind->oneFour)
    return 0.0;
  else {
    double i = b1 * cos(theta1) - b2 + b3 * cos(theta2);
    double j = b3 * sin(theta2) * sin(phi);
    double k = -b1 * sin(theta1) + b3 * sin(theta2) * cos(phi);
    return (i * i + j * j + k * k);
  }
}

void TrialMol::SetCavMatrix(const XYZArray &matrix) {
  matrix.CopyRange(cavMatrix, 0, 0, 3);
}

void TrialMol::SetSeed(const XYZ &coords, const XYZ &cav, const bool inCav,
                       const bool fixCOM, const bool rotBB) {
  cavityCenter = coords;
  cavity = cav;
  comInCav = inCav;
  comFix = fixCOM;
  rotateBB = rotBB;
}

void TrialMol::SetSeed(const bool inCav, const bool fixCOM, const bool rotBB) {
  comInCav = inCav;
  comFix = fixCOM;
  rotateBB = rotBB;
}

void TrialMol::SetBackBone(const int bb[2]) {
  backbone[0] = bb[0];
  backbone[1] = bb[1];
}

XYZ TrialMol::GetCOM() {
  XYZ tcom;
  uint atomNumber = tCoords.Count();
  XYZArray temp(tCoords);
  axes->UnwrapPBC(temp, box, tCoords.Get(0));
  tCoords = temp;

  for (uint p = 0; p < atomNumber; p++) {
    tcom += temp.Get(p);
  }
  tcom *= (1.0 / (double)(atomNumber));
  // Unwrap with respect to COM
  axes->UnwrapPBC(tCoords, box, tcom);

  return tcom;
}

} // namespace cbmc
