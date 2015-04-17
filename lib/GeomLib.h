/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef GEOM_LIB_H
#define GEOM_LIB_H

#include <math.h> //For sqrt, fabs, M_PI
#include <limits> //for double limits

#include "BasicTypes.h" //For uint, XYZ

/////////////////////////////////////////////////////////////
//  DEFINES  //
///////////////

//Standard way to get pi constant on most platforms
#define _USE_MATH_DEFINES

//In case that didn't work
#ifndef M_PI
//From Mathematica: 
//N[Pi, 75]
#define M_PI \
   3.14159265358979323846264338327950288419716939937510582097494459230781640629
#endif
#ifndef M_PI_2
//From Mathematica:
//N[Pi/2, 75]
#define M_PI_2 \
   1.57079632679489661923132169163975144209858469968755291048747229615390820314
#endif
#ifndef M_PI_4
#define M_PI_4 \
   0.785398163397448309615660845819875721049292349843776455243736148076954101572
#endif

#define DEG_TO_RAD (M_PI/180.0)
#define RAD_TO_DEG (180.0/M_PI)

/////////////////////////////////////////////////////////////
//  FUNCTIONS  //
/////////////////

namespace geom
{
   inline double RadToDeg(const double ang) { return ang * RAD_TO_DEG; }
   inline double DegToRad(const double ang) { return ang * DEG_TO_RAD; }

   //NOTE:
   //Length functions are now directly in XYZ type.

   //Returns cross product of two vectors A and B that share a vertex.
   //NOTE: 
   //To use this on three topologically connected points, we must
   //shift all three points such that the shared vertex is at the origin.
   inline XYZ Cross(const double x1, const double y1, const double z1, 
                    const double x2, const double y2, const double z2)
   { return XYZ(y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2); }

   //Wrapper for cross product using vector types
   inline XYZ Cross(XYZ const& v1, XYZ const& v2)
   { return Cross(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z); }

   //Geometric dot product
   inline double Dot(XYZ const& v1, XYZ const& v2)
   { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }

   //Generates angle between two vectors sharing a common vertex.
   //NOTE:
   //Vectors must be pre-shifted to the origin.
   inline double Theta(XYZ const& v1, XYZ const& v2)
   { return acos(Dot(v1, v2) / sqrt(v1.LengthSq() * v2.LengthSq())); }

   //Returns the dihedral angle between v1 and v3 along v2
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
   inline double Phi(XYZ const& b1, XYZ const& b2, XYZ const& b3)
   {
      XYZ n1 = Cross(b1, b2).Normalize(), n2 = Cross(b2, b3).Normalize(),
         m1 = Cross(n1, XYZ(b2).Normalize()).Normalize();
      return atan2(Dot(m1, n2), Dot(n1, n2));
   }

}

#endif /*GEOM_LIB_H*/

