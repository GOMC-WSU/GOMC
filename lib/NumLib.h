/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef NUMERIC_LIB_H
#define NUMERIC_LIB_H

#include "BasicTypes.h" //For uint, XYZ
#include <iostream>
#include <limits> //for double limits
#include <vector> //for vector average

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623158e+308
#endif

#ifndef SMALL_WEIGHT
#define SMALL_WEIGHT 1.0e-38
#endif
namespace num {
static const double qqFact = 167103.208067979;
static const double MIN_EXP_NONZERO_VAL = -708.4;
static const double BIGNUM = DBL_MAX;
static const uint VDW_STD_KIND = 0, VDW_SHIFT_KIND = 1, VDW_SWITCH_KIND = 2;

inline bool approximatelyEqual(double a, double b, double epsilon) {
  if (abs(a) < 1.0 || abs(b) < 1.0) {
    return abs(a - b) <= epsilon;
  } else {
    return std::abs(a - b) <=
           ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * epsilon);
  }
}

template <typename T> inline void BoundGt(double &val, const double bound) {
  if (val > bound)
    val = bound;
}

template <typename T> inline void BoundLt(double &val, const double bound) {
  if (val < bound)
    val = bound;
}

template <typename T> inline void BoundNZDecimal(T &val, const int mult) {
  BoundLt<T>(val, std::numeric_limits<T>::min() * mult);
}

template <typename T>
inline void Bound(T &val, const double lower, const double upper) {
  BoundLt<T>(val, lower);
  BoundGt<T>(val, upper);
}

// Arithmetic mean.
inline double MeanA(const uint v1, const uint v2) { return (v1 + v2) / 2.0; }
// Geometric mean.
inline double MeanG(const double v1, const double v2) { return sqrt(v1 * v2); }
// Arithmetic mean.
inline double MeanA(std::vector<double> const &v1,
                    std::vector<double> const &v2, const uint ix1,
                    const uint ix2) {
  return (v1[ix1] + v2[ix2]) * 0.5;
}
// Arithmetic mean.
inline double MeanA(std::vector<uint> const &v1, std::vector<uint> const &v2,
                    const uint ix1, const uint ix2) {
  return ((double)(v1[ix1] + v2[ix2])) * 0.5;
}
// Geometric mean.
inline double MeanG(std::vector<double> const &v1,
                    std::vector<double> const &v2, const uint ix1,
                    const uint ix2) {
  return sqrt(v1[ix1] * v2[ix2]);
}

inline double MeanG(std::vector<uint> const &v1, std::vector<uint> const &v2,
                    const uint ix1, const uint ix2) {
  return sqrt(v1[ix1] * v2[ix2]);
}

// return n!
inline double Factorial(const uint n) {
  double result = 1.0;
  for (uint i = 2; i <= n; i++) {
    result *= i;
  }
  return result;
}

// return (n+count)!/n!
inline double Factorial(const uint n, const uint count) {
  double result = 1.0;
  for (uint i = 1; i <= count; i++) {
    result *= n + i;
  }
  return result;
}

template <class Type> inline Type Sq(const Type v) { return v * v; }
template <class Type> inline Type Cb(const Type v) { return v * v * v; }
template <class Type> inline void Cb(Type &s, Type &c, const Type v) {
  s = v * v;
  c = s * v;
}

inline double POW(const double d2, const double d4, const double d6, uint e) {
  double result = (e & 0x1 ? sqrt(d2) : 1.0);
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

// Class to define the function used in Zbrent
class Exp6Fun {
public:
  Exp6Fun(const float a, const float s, const float r = 0.0)
      : sigma(s), alpha(a), rmin(r) {}
  virtual ~Exp6Fun(){};
  virtual float operator()(float x) = 0;

  float sigma, alpha, rmin;
};

class RminFun : public Exp6Fun {
public:
  RminFun(double a, double s) : Exp6Fun(a, s) {}
  virtual ~RminFun(){};
  virtual float operator()(float x) {
    double rep = (6.0 / alpha) * exp(alpha * (1.0 - sigma / x));
    double at = pow(x / sigma, 6.0);
    return (float)(rep - at);
  }
};

class RmaxFun : public Exp6Fun {
public:
  RmaxFun(double a, double s, double r) : Exp6Fun(a, s, r) {}
  virtual ~RmaxFun(){};
  virtual float operator()(float x) {
    double rep = (-1.0 / rmin) * exp(alpha * (1.0 - x / rmin));
    double at = pow(rmin / x, 6.0) / x;
    return (float)(rep + at);
  }
};

// Using Brent’s method, find the root of a function func known to lie between
// x1 and x2. The root, returned as zbrent, will be refined until its accuracy
// is tol. Brent, R.P. 1973, Algorithms for Minimization without Derivatives
// (Englewood Cliffs, NJ: Prentice-Hall) Forsythe, G.E., Malcolm, M.A., and
// Moler, C.B. 1977, Computer Methods for Mathematical Computations (Englewood
// Cliffs, NJ: Prentice-Hall)
inline double Zbrent(Exp6Fun *func, float x1, float x2, float tol) {
  const int ITMAX = 100;    // Maximum allowed number of iterations.
  const float EPS = 3.0e-8; // Machine floating-point precision.
  int iter;
  float a = x1, b = x2, c = x2;
  float d, e, min1, min2;
  float fa = (*func)(a), fb = (*func)(b);
  float fc, p, q, r, s, tol1, xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    std::cout << "Root must be bracketed in zbrent!\n";
    exit(EXIT_FAILURE);
  }

  fc = fb;

  for (iter = 1; iter <= ITMAX; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a; // Rename a, b, c and adjust bounding interval
      fc = fa;
      e = d = b - a;
    }

    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol; // Convergence check.
    xm = 0.5 * (c - b);
    if (fabs(xm) <= tol1 || fb == 0.0) {
      return (double)b;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb / fa; // Attempt inverse quadratic interpolation.
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }

      if (p > 0.0) {
        q = -q; // Check whether in bounds.
      }

      p = fabs(p);
      min1 = 3.0 * xm * q - fabs(tol1 * q);
      min2 = fabs(e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        e = d; // Accept interpolation.
        d = p / q;
      } else {
        d = xm; // Interpolation failed, use bisection.
        e = d;
      }
    } else {
      // Bounds decreasing too slowly, use bisection.
      d = xm;
      e = d;
    }

    a = b; // Move last best guess to a.
    fa = fb;
    if (fabs(d) > tol1) { // Evaluate new trial root.
      b += d;
    } else {
      // b += SIGN(tol1, xm);
      b += ((xm) >= 0.0 ? fabs(tol1) : -fabs(tol1));
    }

    fb = (*func)(b);
  }
  std::cout << "Maximum number of iterations exceeded in zbrent!\n";
  exit(EXIT_FAILURE);

  return 0.0; // Never get here.
}

} // namespace num

#endif /*NUMERIC_LIB_H*/
