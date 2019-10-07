#include "CubicSpline.h"
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cassert>

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