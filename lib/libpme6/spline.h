#ifndef SPLINE_H
#define SPLINE_H

#include <complex>

struct spline4 {
  spline4();

  double operator()(double x) const;

  double diff(double x) const;

  std::complex<double> factor(double m, double K) const;
  
  static const unsigned int order = 4;
};

#endif
