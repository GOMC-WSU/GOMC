#include <cstdlib>
#include <cstring>

#include <sys/types.h>

#include "linear-algebra.h"

namespace ewald {

void matrix_vector_product(double const* A, double* X) {
  double Y[] = { 0.0, 0.0, 0.0 };

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      Y[i] += A[3*i+j] * X[j];

  memcpy(X, Y, 3 * sizeof(double));
}

void matrix_vector_product_transpose(double const* A, double* X) {
  double Y[] = { 0.0, 0.0, 0.0 };

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      Y[i] += A[3*j+i] * X[j];

  memcpy(X, Y, 3 * sizeof(double));
}

void matrix_matrix_product(double const* A, double const* B, double* C) {
  memset(C, 0, 3 * 3 * sizeof(double));

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        C[3*i+j] += A[3*i+k] * B[3*k+j];
}

double determinant(double const* A) {
  return A[3*0+0] * (A[3*1+1] * A[3*2+2] - A[3*1+2] * A[3*2+1])
       - A[3*0+1] * (A[3*1+0] * A[3*2+2] - A[3*1+2] * A[3*2+0])
       + A[3*0+2] * (A[3*1+0] * A[3*2+1] - A[3*1+1] * A[3*2+0]);
}

void matrix_inverse(double const* A, double* Ainv) {
  const double invdet = 1.0 / determinant(A);
  memset(Ainv, 0, 3 * 3 * sizeof(double));
  Ainv[3*0+0] = ( A[3*1+1]*A[3*2+2] - A[3*1+2]*A[3*2+1]) * invdet;
  Ainv[3*0+1] = (-A[3*0+1]*A[3*2+2] + A[3*0+2]*A[3*2+1]) * invdet;
  Ainv[3*0+2] = ( A[3*0+1]*A[3*1+2] - A[3*0+2]*A[3*1+1]) * invdet;
  Ainv[3*1+0] = (-A[3*1+0]*A[3*2+2] + A[3*1+2]*A[3*2+0]) * invdet;
  Ainv[3*1+1] = ( A[3*0+0]*A[3*2+2] - A[3*0+2]*A[3*2+0]) * invdet;
  Ainv[3*1+2] = (-A[3*0+0]*A[3*1+2] + A[3*0+2]*A[3*1+0]) * invdet;
  Ainv[3*2+0] = ( A[3*1+0]*A[3*2+1] - A[3*1+1]*A[3*2+0]) * invdet;
  Ainv[3*2+1] = (-A[3*0+0]*A[3*2+1] + A[3*0+1]*A[3*2+0]) * invdet;
  Ainv[3*2+2] = ( A[3*0+0]*A[3*1+1] - A[3*0+1]*A[3*1+0]) * invdet;
}

double dot_product(double const* u, double const* v) {
  double result = 0.0;

  for (size_t k = 0; k < 3; ++k)
    result += u[k] * v[k];

  return result;
}

}
