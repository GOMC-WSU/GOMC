#pragma once
#include <vector>

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