#pragma once
#include <vector>

class CubicSpline {
public:
  CubicSpline(int tableLength, double tableSpacing, double start);
  CubicSpline() { }
  void Reconstruct(int tableLength, double tableSpacing, double start);
  void InitializeSpecificPoint(double functionValue0, double functionValue1, double derivativeValue0, double derivativeValue1, int index);
  double operator()(double value);
  double GetDerivativeValue(double value);

private:
  std::vector<double> Y;
  std::vector<double> F;
  std::vector<double> G;
  std::vector<double> H;
  int tableLength;
  double tableSpacing;
  double tableStart;
};