/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef NUMERIC_LIB_H
#define NUMERIC_LIB_H

#include <limits> //for double limits
#include <vector> //for vector average
#include "BasicTypes.h" //For uint, XYZ

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623158e+308
#endif
#ifndef FLT_MAX
#define FLT_MAX 3.402823e+38
#endif

namespace num
{
static const real dbl_margin = 0.00001;
static const real qqFact = 167000.00;
#if defined(GOMC_DOUBLE)
static const double BIGNUM = DBL_MAX;
#else
static const float BIGNUM = FLT_MAX;
#endif
static const uint VDW_STD_KIND = 0, VDW_SHIFT_KIND = 1, VDW_SWITCH_KIND = 2;

template <typename T>
inline void BoundGt(real & val, const T bound)
{
  if (val > bound) val = bound;
}

template <typename T>
inline void BoundLt(real & val, const T bound)
{
  if (val < bound) val = bound;
}

template <typename T>
inline void BoundNZDecimal(T & val, const int mult)
{
  BoundLt<T>(val, std::numeric_limits<T>::min() * mult);
}

template <typename T>
inline void Bound(T & val, const T lower, const T upper)
{
  BoundLt<T>(val, lower);
  BoundGt<T>(val, upper);
}

//Arithmetic mean.
inline real MeanA(const uint v1, const uint v2)
{
  return (v1 + v2) / 2.0;
}
//Geometric mean.
inline real MeanG(const real v1, const real v2)
{
  return sqrt(v1 * v2);
}
//Arithmetic mean.
inline real MeanA(std::vector<real> const& v1,
                    std::vector<real> const& v2,
                    const uint ix1, const uint ix2)
{
  return (v1[ix1] + v2[ix2]) * 0.5;
}
//Arithmetic mean.
inline real MeanA(std::vector<uint> const& v1,
                    std::vector<uint> const& v2,
                    const uint ix1, const uint ix2)
{
#ifdef MIE_INT_ONLY
  return (v1[ix1] + v2[ix2]) / 2;
#else
  return ((real)(v1[ix1] + v2[ix2])) / 2.0;
#endif
}
//Geometric mean.
inline real MeanG(std::vector<real> const& v1,
                    std::vector<real> const& v2,
                    const uint ix1, const uint ix2)
{
  return sqrt(v1[ix1] * v2[ix2]);
}

//return n!
inline real Factorial(const uint n)
{
  real result = 1.0;
  for(uint i = 2; i <= n; i++) {
    result *= i;
  }
  return result;
}

//return (n+count)!/n!
inline real Factorial(const uint n, const uint count)
{
  real result = 1.0;
  for(uint i = 1; i <= count; i++) {
    result *= n + i;
  }
  return result;
}

template <class Type>
inline Type Sq(const Type v)
{
  return v * v;
}
template <class Type>
inline Type Cb(const Type v)
{
  return v * v * v;
}
template <class Type>
inline void Cb(Type & s, Type & c, const Type v)
{
  s = v * v;
  c = s * v;
}

inline real POW(const real d2, const real d4, const real d6,
                  uint e)
{
  real result = (e & 0x1 ? sqrt(d2) : 1.0);
  e >>= 1;
  switch (e) {
  case 0:
    break;
  case 1:
    result *= d2;
    break;
  case 2:
    result *= d4;
    break;
  case 3:
    result *= d6;
    break;
  case 4:
    result *= Sq(d4);
    break;
  case 5:
    result *= d4 * d6;
    break;
  case 6:
    result *= Sq(d6);
    break;
  case 7:
    result *= Sq(d6) * d2;
    break;
  case 8:
    result *= Sq(d6) * d4;
    break;
  case 9:
    result *= Sq(d6) * d6;
    break;
  case 10:
    result *= Sq(d6) * Sq(d4);
    break;
  case 11:
    result *= Sq(Sq(d4)) * d6;
    break;
  case 12:
    result *= Sq(Sq(d6));
    break;
  case 13:
    result *= Sq(Sq(d6)) * d2;
    break;
  case 14:
    result *= Sq(Sq(d6)) * d4;
    break;
  case 15:
    result *= Sq(Sq(d6)) * d6;
    break;
  case 16:
    result *= Sq(Sq(Sq(d4)));
    break;
  case 17:
    result *= Sq(Sq(Sq(d4))) * d2;
    break;
  case 18:
    result *= Sq(Sq(d6) * d6);
    break;
  case 19:
    result *= Sq(Sq(d6) * d6) * d2;
    break;
  case 20:
    result *= Sq(Sq(d6 * d4));
    break;
  case 21:
    result *= Sq(Sq(d6 * d4)) * d2;
    break;
  case 22:
    result *= Sq(Sq(Sq(d4)) * d6);
    break;
  case 23:
    result *= Sq(Sq(Sq(d4)) * d6) * d2;
    break;
  case 24:
    result *= Sq(Sq(Sq(d6)));
    break;
  case 25:
    result *= Sq(Sq(Sq(d6))) * d2;
    break;
  }
  return result;
}
}

#endif /*NUMERIC_LIB_H*/
