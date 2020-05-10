#include <cassert>

#include <algorithm>

#include "cube.h"

namespace ewald {

using std::complex;

template <typename T>
cube<T>::cube(int d0, int d1, int d2)
    : dim0(d0), dim1(d1), dim2(d2) {
  memory = new T[dim0 * dim1 * dim2];
}

template <typename T>
cube<T>::~cube() {
  delete [] memory;
}

template <typename T>
void cube<T>::zeros() {
  std::fill_n(memory, dim0 * dim1 * dim2, T(0));
}

template <typename T>
T cube<T>::operator()(int k0, int k1, int k2) const {
  assert(0 <= k0 && k0 < dim0);
  assert(0 <= k1 && k1 < dim1);
  assert(0 <= k2 && k2 < dim2);

  return memory[k0 * dim2 * dim1 + k1 * dim2 + k2];
}

template <typename T>
T& cube<T>::operator()(int k0, int k1, int k2) {
  assert(0 <= k0 && k0 < dim0);
  assert(0 <= k1 && k1 < dim1);
  assert(0 <= k2 && k2 < dim2);

  return memory[k0 * dim2 * dim1 + k1 * dim2 + k2];
}

template <typename T>
cube<T>& cube<T>::operator*=(cube const& rhs) {
  assert(rhs.dim0 == dim0);
  assert(rhs.dim1 == dim1);
  assert(rhs.dim2 == dim2);

  for (int k0 = 0; k0 < dim0; ++k0)
    for (int k1 = 0; k1 < dim1; ++k1)
      for (int k2 = 0; k2 < dim2; ++k2)
        operator()(k0, k1, k2) *= rhs(k0, k1, k2);

  return *this;
}

template <typename T>
T* cube<T>::memptr() {
  return memory;
}

}
