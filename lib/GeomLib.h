/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef GEOM_LIB_H
#define GEOM_LIB_H

// Standard way to get pi constant on most platforms
// Needs to be defined _before_ including cmath
// so that the PI constants come from cmath
//#define _USE_MATH_DEFINES
//#include <cmath>  //For sqrt, fabs, M_PI

#include "BasicTypes.h" //For uint, XYZ
#include "XYZArray.h"
#include <limits> //for double limits

/////////////////////////////////////////////////////////////
//  DEFINES  //
///////////////

// Just in case any of these weren't included from cmath
#ifndef M_PI
// From Mathematica:
// N[Pi, 75]
#define M_PI                                                                   \
  3.14159265358979323846264338327950288419716939937510582097494459230781640629
#endif
#ifndef M_1_PI
// Reciprocal of PI:
#define M_1_PI 1.0 / M_PI
#endif
#ifndef M_PI_2
// From Mathematica:
// N[Pi/2, 75]
#define M_PI_2                                                                 \
  1.57079632679489661923132169163975144209858469968755291048747229615390820314
#endif
#ifndef M_PI_4
#define M_PI_4                                                                 \
  0.785398163397448309615660845819875721049292349843776455243736148076954101572
#endif
#ifndef M_2_SQRTPI
#define M_2_SQRTPI 2.0 / std::sqrt(M_PI)
#endif

#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 * M_1_PI) // Same as 180/PI

/////////////////////////////////////////////////////////////
//  FUNCTIONS  //
/////////////////

namespace geom {
inline double RadToDeg(const double ang) { return ang * RAD_TO_DEG; }
inline double DegToRad(const double ang) { return ang * DEG_TO_RAD; }

// NOTE:
// Length functions are now directly in XYZ type.

// Returns cross product of two vectors A and B that share a vertex.
// NOTE:
// To use this on three topologically connected points, we must
// shift all three points such that the shared vertex is at the origin.
inline XYZ Cross(const double x1, const double y1, const double z1,
                 const double x2, const double y2, const double z2) {
  return XYZ(y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2);
}

// Wrapper for cross product using vector types
inline XYZ Cross(XYZ const &v1, XYZ const &v2) {
  return Cross(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z);
}

// Geometric dot product
inline double Dot(XYZ const &v1, XYZ const &v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// Geometric dot product for atomcoordinates.K
inline double Dot(const uint atom, const double kx, const double ky,
                  const double kz, const XYZArray &Coords) {
  return (Coords.x[atom] * kx + Coords.y[atom] * ky + Coords.z[atom] * kz);
}

// Find two vectors that are perpendicular to V
inline void SetBasis(XYZArray &basis, XYZ const &V) {
  XYZ v1 = V;
  v1.Normalize();
  XYZ v2;
  if (fabs(v1.x) < 0.8) {
    // v3 will be v1 x the standard X unit vec
    v2 = XYZ(1.0, 0.0, 0.0);
  } else {
    // v3 will be v1 x the standard Y unit vec
    v2 = XYZ(0.0, 1.0, 0.0);
  }
  XYZ v3 = Cross(v1, v2);
  v3.Normalize();
  // v2 is unit vec perpendicular to both v3 and v2
  v2 = Cross(v3, v1);
  // set v1 az z axis of cell basis
  basis.Set(0, v3);
  basis.Set(1, v2);
  basis.Set(2, v1);
}

inline void TransposeMatrix(XYZArray &Inv, XYZArray const &matrix) {
  Inv.x[0] = matrix.x[0];
  Inv.x[1] = matrix.y[0];
  Inv.x[2] = matrix.z[0];

  Inv.y[0] = matrix.x[1];
  Inv.y[1] = matrix.y[1];
  Inv.y[2] = matrix.z[1];

  Inv.z[0] = matrix.x[2];
  Inv.z[1] = matrix.y[2];
  Inv.z[2] = matrix.z[2];
}

// Calculate transform of A using dot product
inline XYZ Transform(const XYZArray &basis, const XYZ &A) {
  XYZ temp;
  temp.x = A.x * basis.x[0] + A.y * basis.x[1] + A.z * basis.x[2];
  temp.y = A.x * basis.y[0] + A.y * basis.y[1] + A.z * basis.y[2];
  temp.z = A.x * basis.z[0] + A.y * basis.z[1] + A.z * basis.z[2];
  return temp;
}

// Generates angle between two vectors sharing a common vertex.
// NOTE:
// Vectors must be pre-shifted to the origin.
inline double Theta(XYZ const &v1, XYZ const &v2) {
  double value = Dot(v1, v2) / sqrt(v1.LengthSq() * v2.LengthSq());
  if (value < -1.00)
    return M_PI;
  else
    return acos(value);
}

// Returns the dihedral angle between v1 and v3 along v2
// Derived from http://math.stackexchange.com/questions/47059
// See: answer from Rahul Narain
//
//             ^            ^ b3
//            / b3         /
//    ---b2->/        b2  *)--- this angle ^ is phi
//   ^                    |     (NOTE: this is important for chirality!)
//  /                     | b1
// / b1                   v
//
// 1.)  Normal vectors to component planes:
// n1: < b1 x b2 >  n2: < b2 x b3 > ;
// 2.) Orthonormal frame basis:
// (n1, <b2>, m1) where m1: n1 x <b2>
// 3.) phi = atan2(y,x)
//
// NOTE:
// cos-1 route avoided as it's inaccurate near 0 or pi.
inline double Phi(XYZ const &b1, XYZ const &b2, XYZ const &b3) {
  XYZ n1 = Cross(b1, b2).Normalize(), n2 = Cross(b2, b3).Normalize(),
      m1 = Cross(n1, XYZ(b2).Normalize()).Normalize();
  return atan2(Dot(m1, n2), Dot(n1, n2));
}

} // namespace geom

#endif /*GEOM_LIB_H*/
