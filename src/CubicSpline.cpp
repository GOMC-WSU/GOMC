#include "CubicSpline.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cassert>

#define HALF_DOUBLE_MAX 8.988466e+307

CubicSpline::~CubicSpline()
{
  delete [] Y;
  delete [] F;
  delete [] G;
  delete [] H;
}

void CubicSpline::Reconstruct(int tableLength, uint kindTotalSq, double tableSpacing, double start)
{
  this->tableLength = tableLength;
  this->tableSpacing = tableSpacing;
  this->tableStart = start;

  Y = new double[tableLength * kindTotalSq];
  F = new double[tableLength * kindTotalSq];
  G = new double[tableLength * kindTotalSq];
  H = new double[tableLength * kindTotalSq];
}

void CubicSpline::InitializeSpecificPoint(double functionValue0, double functionValue1, double derivativeValue0, double derivativeValue1, int index, uint kind)
{
  Y[kind * tableLength + index] =  functionValue0;
  F[kind * tableLength + index] =  tableSpacing * derivativeValue0;
  G[kind * tableLength + index] =  3.0*( functionValue1 - functionValue0) - tableSpacing * (derivativeValue1 + 2.0 * derivativeValue0);
  H[kind * tableLength + index] = -2.0*( functionValue1 - functionValue0) + tableSpacing * (derivativeValue1 + derivativeValue0);
}

double CubicSpline::operator()(double value, uint kind)
{
  int index = static_cast<int>((value - tableStart) / tableSpacing);
  index = kind * tableLength + index;
  double x = tableSpacing * index + tableStart;
  double eps = value - x;

  double temp1 = H[index] * eps + G[index];
  double temp2 = temp1 * eps + F[index];
  double temp3 = temp2 * eps + Y[index];
  return (isnan(temp3)) ? HALF_DOUBLE_MAX : temp3;
}

double CubicSpline::ReturnArray(std::vector<double> &values, int size, std::vector<int> &kinds, double threshold)
{
  double sum = 0.0;
  for (int i = 0; i < size; i++) {
    if (values[i] > threshold)
      continue;

    int kind = kinds[i];
    int index = static_cast<int>((values[i] - tableStart) / tableSpacing);
    index = kind * tableLength + index;
    double x = tableSpacing * index + tableStart;
    double eps = values[i] - x;

    double temp1 = H[index] * eps + G[index];
    double temp2 = temp1 * eps + F[index];
    double temp3 = temp2 * eps + Y[index];
    if (isnan(temp3)) {
      return HALF_DOUBLE_MAX;
    }
    sum = sum + temp3;
  }
  return sum;
}

double CubicSpline::GetDerivativeValue(double value, uint kind)
{
  int index = static_cast<int>((value - tableStart) / tableSpacing);
  index = kind * tableLength + index;
  double eps = value - (tableSpacing * index + tableStart);

  double temp1 = 3.0 * H[index];
  double temp2 = 2.0 * G[index];

  double temp3 = temp1 * eps + temp2;
  double temp4 = temp3 * eps + F[index];
  return (temp4 / tableSpacing);
}