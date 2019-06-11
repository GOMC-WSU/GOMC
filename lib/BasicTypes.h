/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <cstddef>
#include <math.h>
#include <ostream>

typedef unsigned int uint;
typedef unsigned long int ulong;

#define UNUSED(x) (void)(x)

//single XYZ for use as a temporary and return type
struct XYZ {
  real x, y, z;

  XYZ() : x(0.0), y(0.0), z(0.0) {}
  XYZ(real xVal, real yVal, real zVal) : x(xVal), y(yVal), z(zVal) {}

  XYZ& operator=(XYZ const& rhs)
  {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
    return *this;
  }
  XYZ& operator+=(XYZ const& rhs)
  {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }
  XYZ& operator-=(XYZ const& rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }
  XYZ& operator*=(XYZ const& rhs)
  {
    x *= rhs.x;
    y *= rhs.y;
    z *= rhs.z;
    return *this;
  }
  XYZ& operator/=(XYZ const& rhs)
  {
    x /= rhs.x;
    y /= rhs.y;
    z /= rhs.z;
    return *this;
  }

  XYZ& operator*=(const real a)
  {
    x *= a;
    y *= a;
    z *= a;
    return *this;
  }

  XYZ operator+(XYZ const& rhs) const
  {
    return XYZ(*this) += rhs;
  }
  XYZ operator-(XYZ const& rhs) const
  {
    return XYZ(*this) -= rhs;
  }
  XYZ operator*(XYZ const& rhs) const
  {
    return XYZ(*this) *= rhs;
  }
  XYZ operator/(XYZ const& rhs) const
  {
    return XYZ(*this) /= rhs;
  }


  XYZ operator*(const real a) const
  {
    return XYZ(*this) *= a;
  }

  XYZ operator-() const
  {
    return XYZ(*this) * -1.0;
  }

  void Inverse()
  {
    x = 1.0 / x;
    y = 1.0 / y;
    z = 1.0 / z;
  }

  real Length() const
  {
    return sqrt(LengthSq());
  }
  real LengthSq() const
  {
    return x * x + y * y + z * z;
  }

  XYZ& Normalize()
  {
    *this *= (1 / Length());
    return *this;
  }

  real Max() const
  {
    real m = x;
    if(y > m)
      m = y;
    if(z > m)
      m = z;
    return m;
  }

  real Min() const
  {
    real m = x;
    if(y < m)
      m = y;
    if(z < m)
      m = z;
    return m;
  }
};

inline std::ostream& operator << (std::ostream & stream, const XYZ& p)
{
  stream << "[" << p.x << ", " << p.y << ", " << p.z << "]";
  return stream;
}

#endif /*BASIC_TYPES_H*/
