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

void CubicSpline::Reconstruct(int tableLength, double tableSpacing, double start)
{
  this->tableLength = tableLength;
  this->tableSpacing = tableSpacing;
  this->tableStart = start;

  Y = new double[tableLength];
  F = new double[tableLength];
  G = new double[tableLength];
  H = new double[tableLength];
}

void CubicSpline::InitializeSpecificPoint(double functionValue0, double functionValue1, double derivativeValue0, double derivativeValue1, int index)
{
  Y[index] =  functionValue0;
  F[index] =  tableSpacing * derivativeValue0;
  G[index] =  3.0*( functionValue1 - functionValue0) - tableSpacing * (derivativeValue1 + 2.0 * derivativeValue0);
  H[index] = -2.0*( functionValue1 - functionValue0) + tableSpacing * (derivativeValue1 + derivativeValue0);
}

double CubicSpline::operator()(double value)
{
  int index = static_cast<int>((value - tableStart) / tableSpacing);
  double x = tableSpacing * index + tableStart;
  double eps = value - x;

  double temp1 = H[index] * eps + G[index];
  double temp2 = temp1 * eps + F[index];
  double temp3 = temp2 * eps + Y[index];
  return (isnan(temp3)) ? HALF_DOUBLE_MAX : temp3;
}

double CubicSpline::GetDerivativeValue(double value)
{
  int index = static_cast<int>((value - tableStart) / tableSpacing);
  double eps = value - (tableSpacing * index + tableStart);

  double temp1 = 3.0 * H[index];
  double temp2 = 2.0 * G[index];

  double temp3 = temp1 * eps + temp2;
  double temp4 = temp3 * eps + F[index];
  return (temp4 / tableSpacing);
}