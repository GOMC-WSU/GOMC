/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.8
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <cstddef>
#include <math.h>

typedef unsigned int uint;
typedef unsigned long int ulong;

#define UNUSED(x) (void)(x)

//single XYZ for use as a temporary and return type
struct XYZ
{
   double x, y, z;

   XYZ() : x(0.0), y(0.0), z(0.0) {}
   XYZ(double xVal, double yVal, double zVal) : x(xVal), y(yVal), z(zVal) {}

   XYZ& operator+=(XYZ const& rhs)
   { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
   XYZ& operator-=(XYZ const& rhs)
   { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
   XYZ& operator*=(XYZ const& rhs)
   { x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this; }

   XYZ& operator*=(const double a)
   { x *= a; y *= a; z *= a; return *this; }

   XYZ operator+(XYZ const& rhs) const
   { return XYZ(*this) += rhs; }
   XYZ operator-(XYZ const& rhs) const
   { return XYZ(*this) -= rhs; }
   XYZ operator*(XYZ const& rhs) const
   { return XYZ(*this) *= rhs; }


   XYZ operator*(const double a) const
   { return XYZ(*this) *= a; }

   XYZ operator-() const { return XYZ(*this) * -1.0; }

   double Length() const { return sqrt(LengthSq()); }
   double LengthSq() const { return x * x + y * y + z * z; }
   XYZ& Normalize()
   {
      *this *= (1 / Length());
      return *this;
   }

};

#endif /*BASIC_TYPES_H*/
