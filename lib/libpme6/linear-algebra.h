#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

namespace ewald {

void matrix_vector_product(double const* A, double* X);

void matrix_vector_product_transpose(double const* A, double* X);

void matrix_matrix_product(double const* A, double const* B, double* C);

double determinant(double const* A);

void matrix_inverse(double const* A, double* Ainv);

double dot_product(double const* u, double const* v);

}

#endif
