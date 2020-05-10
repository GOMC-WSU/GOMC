#include <cmath>

#include <complex>

#include "spline.h"

using std::floor;
using std::complex;

const complex<double> I(0.0, 1.0);

spline4::spline4() {}

double spline4::operator()(double x) const {
  if (0.0 < x && x < 1.0)
    return 1.0 / 6.0 * x * x * x;
  else if (1.0 <= x && x < 2.0)
    return 2.0 / 3.0 + (-2.0 + (2.0 - 1.0 / 2.0 * x) * x) * x;
  else if (2.0 <= x && x < 3.0)
    return -22.0 / 3.0 + (10.0 + (-4.0 + 1.0 / 2.0 * x) * x) * x;
  else if (3.0 <= x && x < 4.0)
    return 32.0 / 3.0 + (-8.0 + (2.0 - 1.0 / 6.0 * x) * x) * x;
  else
    return 0.0;
}

double spline4::diff(double x) const {
  if (0.0 < x && x < 1.0)
    return 1.0 / 2.0 * x * x;
  else if (1.0 <= x && x < 2.0)
    return -2.0 + (4.0 - 3.0 / 2.0 * x) * x;
  else if (2.0 <= x && x < 3.0)
    return 10.0 + (-8.0 + 3.0 / 2.0 * x) * x;
  else if (3.0 <= x && x < 4.0)
    return -8.0 + (4.0 - 1.0 / 2.0 * x) * x;
  else
    return 0.0;
}

complex<double> spline4::factor(double m, double K) const {
  const double m_over_K = m / K;

  complex<double> denom = 0.0;
  for (double k = 0.0; k <= order - 2.0; k++)
    denom += operator()(k + 1.0) * exp(2.0 * M_PI * I * m_over_K * k);

  return exp(2.0 * M_PI * I * (order - 1.0) * m_over_K) / denom;
}
