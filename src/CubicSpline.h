#pragma once
#include <vector>
#include <BasicTypes.h>

class CubicSpline {
public:
  CubicSpline() { }
  ~CubicSpline();
  void Reconstruct(int tableLength, uint kindTotalSq, double tableSpacing, double start);
  void InitializeSpecificPoint(double functionValue0, double functionValue1, double derivativeValue0, double derivativeValue1, int index, uint kind);
  double operator()(double value, uint kind);
  double ReturnArray(std::vector<double> &values, int size, std::vector<int> &kinds, double threshold);
  double GetDerivativeValue(double value, uint kind);

private:
  double * Y;
  double * F;
  double * G;
  double * H;
  int tableLength;
  double tableSpacing;
  double tableStart;
};