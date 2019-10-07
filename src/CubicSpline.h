#pragma once
#include <vector>
#include "BasicTypes.h"

#define HALF_DOUBLE_MAX 8.988466e+307

class CubicSpline {
public:
  CubicSpline() { }
  ~CubicSpline();
  void Reconstruct(int tableLength, double tableSpacing, double start);
  void InitializeSpecificPoint(double functionValue0, double functionValue1, double derivativeValue0, double derivativeValue1, int index);
  double operator()(double value);
  double GetDerivativeValue(double value);

private:
  double * Y;
  double * F;
  double * G;
  double * H;
  int tableLength;
  double tableSpacing;
  double tableStart;
};

inline double CubicSpline::operator()(double value)
{
  int index = static_cast<int>((value - tableStart) / tableSpacing);
  double x = tableSpacing * index + tableStart;
  double eps = value - x;

  double temp1 = H[index] * eps + G[index];
  double temp2 = temp1 * eps + F[index];
  double temp3 = temp2 * eps + Y[index];
  return (isnan(temp3)) ? HALF_DOUBLE_MAX : temp3;
}

inline double CubicSpline::GetDerivativeValue(double value)
{
  int index = static_cast<int>((value - tableStart) / tableSpacing);
  double eps = value - (tableSpacing * index + tableStart);

  double temp1 = 3.0 * H[index];
  double temp2 = 2.0 * G[index];

  double temp3 = temp1 * eps + temp2;
  double temp4 = temp3 * eps + F[index];
  return (temp4 / tableSpacing);
}