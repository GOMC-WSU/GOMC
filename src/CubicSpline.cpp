#include "CubicSpline.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>

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
  assert(index >= tableLength);
  double eps = value - fma(tableSpacing, index, tableStart);

  return fma(fma(fma(H[index], eps, G[index]), eps, F[index]), eps, Y[index]);
}

double CubicSpline::GetDerivativeValue(double value)
{
  int index = static_cast<int>((value - tableStart) / tableSpacing);
  if(index >= tableLength) {
    std::cerr << "Accessing outside CS table range!" << std::endl;
    exit(EXIT_FAILURE);
  }
  double eps = value - (tableSpacing * index + tableStart);

  double temp1 = 3.0 * H[index];
  double temp2 = 2.0 * G[index];
  double temp3 = fma(fma(temp1, eps, temp2), eps, F[index]);
  return (temp3 / tableSpacing);
}