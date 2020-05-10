#ifndef CUBE_H
#define CUBE_H

#include <complex>

namespace ewald {

template <typename T>
class cube {
 public:
  cube(int d0, int d1, int d2);
  ~cube();

  void zeros();

  T operator()(int k0, int k1, int k2) const;
  T& operator()(int k0, int k1, int k2);

  cube<T>& operator*=(cube const& rhs);

  T* memptr();

 private:
  const int dim0, dim1, dim2;
  T* memory;
};

typedef cube<double> real_cube;
typedef cube<std::complex<double> > complex_cube;

}

#endif
