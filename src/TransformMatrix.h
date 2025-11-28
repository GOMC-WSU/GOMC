/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef TRANSFORM_MATRIX_H
#define TRANSFORM_MATRIX_H

#include "BasicTypes.h"

class TransformMatrix;
typedef TransformMatrix RotationMatrix;

class TransformMatrix {
public:
  TransformMatrix() { LoadIdentity(); }
  explicit TransformMatrix(double d) { SetDiagonal(d); }
  //! Set the rotation matrix to the identity matrix.
  void LoadIdentity() { SetDiagonal(1.0); }
  void SetDiagonal(double d);

  // returns matrix for rotation of theta radians about axis
  static TransformMatrix FromAxisAngle(double theta, const XYZ &axis);

  // returns matrix for rotation of theta radians about vector described by
  // matrices cross and tensor (more efficient if reusing an axis)
  static TransformMatrix FromAxisAngle(double theta,
                                       const TransformMatrix &cross,
                                       const TransformMatrix &tensor);

  // returns matrix such that (vector x vec2) = matrix * vec2
  // for efficient axis-angle rotation using Rodrigues' Formula
  static TransformMatrix CrossProduct(const XYZ &vector);

  // returns product of vector and vectorT
  // for efficient axis-angle rotation using Rodrigues' Formula
  static TransformMatrix TensorProduct(const XYZ &vector);

  // returns a uniformly random rotation on SO3 if arguments are uniformly
  // random on [0,1]. u1 and u3 can be uniformly random on [0,k] (0 < k < 1)
  // for a smaller but still uniformly random rotation
  static TransformMatrix UniformRandom(double u1, double u2, double u3);

  /*!Add a rotation transform to the matrix.
   */
  void AddRotationX(double angle);
  void AddRotationY(double angle);
  void AddRotationZ(double angle);
  // z-x-z Euler angles - beware ye gimbal lock
  void AddRotation(const XYZ &rotateBy);

  //! Set the matrix to a change-of-basis rotation from unity
  void BasisRotation(const XYZ &u, const XYZ &v, const XYZ &w);

  TransformMatrix operator*(const TransformMatrix &rhs) const;

  //! Apply the transformation to pos
  XYZ Apply(const XYZ &pos) const;

  // returns the inverse=transpose of this matrix
  TransformMatrix Inverse() const;

  // private:
  static const uint N = 3;
  double matrix[N][N];
};

inline void TransformMatrix::SetDiagonal(double d) {
  for (uint i = 0; i < N; ++i) {
    for (uint j = 0; j < N; ++j) {
      matrix[i][j] = 0.0;
    }
    matrix[i][i] = d;
  }
}

inline void TransformMatrix::AddRotation(const XYZ &rotateBy) {
  AddRotationZ(rotateBy.x);
  AddRotationX(rotateBy.y);
  AddRotationZ(rotateBy.z);
}

inline void TransformMatrix::AddRotationX(double angle) {
  TransformMatrix rot;
  double cosine = cos(angle);
  double sine = sin(angle);
  rot.matrix[1][1] = cosine;
  rot.matrix[1][2] = -sine;
  rot.matrix[2][1] = sine;
  rot.matrix[2][2] = cosine;
  *this = rot * *this;
}

inline void TransformMatrix::AddRotationY(double angle) {
  TransformMatrix rot;
  double cosine = cos(angle);
  double sine = sin(angle);
  rot.matrix[0][0] = cosine;
  rot.matrix[0][2] = sine;
  rot.matrix[2][0] = -sine;
  rot.matrix[2][2] = cosine;
  *this = rot * *this;
}

inline void TransformMatrix::AddRotationZ(double angle) {
  TransformMatrix rot;
  double cosine = cos(angle);
  double sine = sin(angle);
  rot.matrix[0][0] = cosine;
  rot.matrix[0][1] = -sine;
  rot.matrix[1][0] = sine;
  rot.matrix[1][1] = cosine;
  *this = rot * *this;
}

inline XYZ TransformMatrix::Apply(const XYZ &pos) const {
  XYZ result;
  result.x = matrix[0][0] * pos.x + matrix[0][1] * pos.y + matrix[0][2] * pos.z;
  result.y = matrix[1][0] * pos.x + matrix[1][1] * pos.y + matrix[1][2] * pos.z;
  result.z = matrix[2][0] * pos.x + matrix[2][1] * pos.y + matrix[2][2] * pos.z;
  return result;
}

inline void TransformMatrix::BasisRotation(const XYZ &u, const XYZ &v,
                                           const XYZ &w) {
  // column vectors are the basis vectors
  matrix[0][0] = u.x;
  matrix[1][0] = u.y;
  matrix[2][0] = u.z;

  matrix[0][1] = v.x;
  matrix[1][1] = v.y;
  matrix[2][1] = v.z;

  matrix[0][2] = w.x;
  matrix[1][2] = w.y;
  matrix[2][2] = w.z;
}

inline TransformMatrix
TransformMatrix::operator*(const TransformMatrix &o) const {
  TransformMatrix result(0.0);
  for (uint i = 0; i < N; ++i) {
    for (uint j = 0; j < N; ++j) {
      for (uint k = 0; k < N; ++k) {
        result.matrix[i][j] += matrix[i][k] * o.matrix[k][j];
      }
    }
  }
  return result;
}

inline TransformMatrix TransformMatrix::Inverse() const {
  TransformMatrix inverse;
  // this method does not work in the general case, but it works when we only
  // have rotations and translations transpose orthogonal rotation section
  for (uint i = 0; i < 3; ++i) {
    for (uint j = 0; j < 3; ++j) {
      inverse.matrix[i][j] = matrix[j][i];
    }
  }
  return inverse;
}

inline TransformMatrix TransformMatrix::CrossProduct(const XYZ &vector) {
  TransformMatrix cross;
  cross.matrix[0][0] = 0.0;
  cross.matrix[0][1] = -(vector.z);
  cross.matrix[0][2] = vector.y;

  cross.matrix[1][0] = vector.z;
  cross.matrix[1][1] = 0.0;
  cross.matrix[1][2] = -(vector.x);

  cross.matrix[2][0] = -(vector.y);
  cross.matrix[2][1] = (vector.x);
  cross.matrix[2][2] = 0.0;

  return cross;
}

inline TransformMatrix TransformMatrix::TensorProduct(const XYZ &u) {
  TransformMatrix tens;
  for (uint i = 0; i < 3; ++i) {
    tens.matrix[0][i] = u.x;
    tens.matrix[1][i] = u.y;
    tens.matrix[2][i] = u.z;
  }
  for (uint i = 0; i < 3; ++i) {
    tens.matrix[i][0] *= u.x;
    tens.matrix[i][1] *= u.y;
    tens.matrix[i][2] *= u.z;
  }
  return tens;
}

inline TransformMatrix
TransformMatrix::FromAxisAngle(double theta, const TransformMatrix &cross,
                               const TransformMatrix &tensor) {
  double c = cos(theta);
  TransformMatrix r(c);
  double s = sin(theta);
  for (uint i = 0; i < N; ++i) {
    for (uint j = 0; j < N; ++j) {
      r.matrix[i][j] += s * cross.matrix[i][j] + (1 - c) * tensor.matrix[i][j];
    }
  }
  return r;
}

inline TransformMatrix TransformMatrix::FromAxisAngle(double theta,
                                                      const XYZ &axis) {
  return FromAxisAngle(theta, CrossProduct(axis), TensorProduct(axis));
}

// Returns a rotation matrix that is uniformly random in S2 if u1/2/3 are
// uniformly random on [0, 1)
inline TransformMatrix TransformMatrix::UniformRandom(double u1, double u2,
                                                      double u3) {
  // 6 transcendentals, 17 double multiplies, 9 double adds

  // method from Arvo (1992)
  // rotate around the pole (0,0,1), then move the pole to a random
  // point in S2 via an inverted Householder reflection
  u1 *= 2.0 * M_PI; // theta - magnitude of polar rotation
  u1 -= M_PI;       // make restricted rotations symmetric for detailed balance
  u2 *= 2.0 * M_PI; // phi   - direction of pole shift
  u3 *= 2.0;        // z     - magnitude of pole shift

  double r = sqrt(u3);
  // Vector used for reflection
  // The reflection matrix is I - V * Transpose(V)
  XYZ V(sin(u2) * r, cos(u2) * r, sqrt(2.0 - u3));

  double sinTh = sin(u1);
  double cosTh = cos(u1);
  // Vector S = Transpose(V) * R
  // R = rotation of theta degrees about z axis
  double Sx = V.x * cosTh - V.y * sinTh;
  double Sy = V.x * sinTh + V.y * cosTh;
  //     Sz = V.z

  // Result Matrix = (V * Transpose(V) - I) R
  //              =  V * Transpose(V) * R - R
  //              =  V * S - R
  // R is a rotation about the pole
  // V * Transpose(V) - I is an inverted reflection to a random point on the
  // sphere
  TransformMatrix result;
  result.matrix[0][0] = V.x * Sx - cosTh;
  result.matrix[0][1] = V.x * Sy - sinTh;
  result.matrix[0][2] = V.x * V.z;

  result.matrix[1][0] = V.y * Sx + sinTh;
  result.matrix[1][1] = V.y * Sy - cosTh;
  result.matrix[1][2] = V.y * V.z;

  result.matrix[2][0] = V.z * Sx;
  result.matrix[2][1] = V.z * Sy;
  result.matrix[2][2] = 1.0 - u3;
  return result;
}
#endif /*TRANSFORM_MATRIX_H*/
