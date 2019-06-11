/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef XYZ_ARRAY_H
#define XYZ_ARRAY_H

#include "BasicTypes.h"
#include <string.h> //for memset, memcpy, etc.
#include <stdio.h> //for memset, memcpy, etc.
#include <cmath>
#include <utility>      //for swap (most modern compilers)
#include <algorithm>      //for swap pre-c++11 compilers
#ifdef _OPENMP
#include <omp.h>
#endif

//Forward declare to give accesss to internal arrays.
class BoxDimensions;



//XYZ structure of arrays for more efficient memory access
class XYZArray
{
public:
  friend class BoxDimensions;

  //Declare a set of coordinates with no data.
  XYZArray(void) : x(NULL), y(NULL), z(NULL), count(0),
    allocDone(false) {}

  //Declare a set of coordinates of an explicit size.
  explicit XYZArray(const uint n)
  {
    allocDone = false;
    Init(n);
  }

  XYZArray(const XYZArray& other);
  XYZArray& operator=(XYZArray other);

  friend void swap(XYZArray& a1, XYZArray& a2);
  ~XYZArray(void)
  {
    Uninit();
  }

  //Wipe out my memberarrays
  void Uninit();

  //Helper function to allocated memory to coordinate arrays.
  //Note: this has been laid out as subarrays of master array.
  void Init(const uint n);

  //Returns number of elements
  uint Count() const
  {
    return count;
  }

  //returns a XYZ copy of row i
  XYZ operator[](const uint i) const
  {
    return Get(i);
  }

  //returns a XYZ copy of row i
  XYZ Get(const uint i) const
  {
    return XYZ(x[i], y[i], z[i]);
  }

  //Returns smallest of x/y/z values in row i
  real Min(const uint i) const
  {
    real result = x[i];
    if (result > y[i])
      result = y[i];
    if (result > z[i])
      result = z[i];
    return result;
  }

  //Returns biggest of x/y/z values in row i
  real Max(const uint i) const
  {
    real result = x[i];
    if (y[i] > result)
      result = y[i];
    if (z[i] > result)
      result = z[i];
    return result;
  }


  //Set a specific coordinate to a value.
  void Set(const uint i, const real a, const real b, const real c)
  {
    x[i] = a;
    y[i] = b;
    z[i] = c;
  }

  void SetRange(const uint start, const uint stop,
                const real a, const real b, const real c);

  //Set a specific coordinate to a value.
  void Set(const uint i, XYZ const& val)
  {
    Set(i, val.x, val.y, val.z);
  }

  void SetRange(const uint start, const uint stop, XYZ const& val);

  ////////////////////////////////////////////////////////////////////

  //Add values in two diferent arrays.
  //Add the value of src[srcI] to this[i]
  void Add(XYZArray & src, const uint i, const uint srcI)
  {
    x[i] += src.x[srcI];
    y[i] += src.y[srcI];
    z[i] += src.z[srcI];
  }

  //return the sum of two rows in two XYZ arrays
  XYZ Sum(const XYZArray & src, const uint i, const uint srcI) const
  {
    return XYZ(x[i] + src.x[srcI], y[i] + src.y[srcI], z[i] + src.z[srcI]);
  }

  //Subtract the value of src[srcI] from this[i]
  void Sub(XYZArray & src, const uint i, const uint srcI)
  {
    x[i] -= src.x[srcI];
    y[i] -= src.y[srcI];
    z[i] -= src.z[srcI];
  }

  //returns row i - row j
  XYZ Difference(const uint i, const uint j) const
  {
    return XYZ(x[i] - x[j], y[i] - y[j], z[i] - z[j]);
  }

  real DifferenceX(const uint i, const uint j) const
  {
    return x[i] - x[j];
  }

  real DifferenceY(const uint i, const uint j) const
  {
    return y[i] - y[j];
  }

  real DifferenceZ(const uint i, const uint j) const
  {
    return z[i] - z[j];
  }

  real Length(const uint i) const
  {
    return sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]);
  }

  //calculate the adjoint and return the determinant
  real AdjointMatrix(XYZArray &Inv);

  //return the difference of two rows in two XYZ arrays
  XYZ Difference(const uint i, XYZArray const& other,
                 const uint otherI) const
  {
    return XYZ( x[i] - other.x[otherI], y[i] - other.y[otherI],
                z[i] - other.z[otherI]);
  }

  real DifferenceX(const uint i, XYZArray const& other,
                     const uint otherI) const
  {
    return x[i] - other.x[otherI];
  }

  real DifferenceY(const uint i, XYZArray const& other,
                     const uint otherI) const
  {
    return y[i] - other.y[otherI];
  }

  real DifferenceZ(const uint i, XYZArray const& other,
                     const uint otherI) const
  {
    return z[i] - other.z[otherI];
  }

  ////////////////////////////////////////////////////////////////////

  //Add row i to row j
  void Add(const uint i, const uint j)
  {
    x[i] += x[j];
    y[i] += y[j];
    z[i] += z[j];
  }

  //Subract row j from row i
  void Sub(const uint i, const uint j)
  {
    x[i] -= x[j];
    y[i] -= y[j];
    z[i] -= z[j];
  }

  //Add values to a single element of the array's values.
  void Add(const uint i, const real a, const real b, const real c)
  {
    x[i] += a;
    y[i] += b;
    z[i] += c;
  }

  //Adds XYZ value val to row i
  void Add(const uint i, XYZ const& val)
  {
    x[i] += val.x;
    y[i] += val.y;
    z[i] += val.z;
  }

  //Sub values from a single element of the array's values.
  void Sub(const uint i, const real a, const real b, const real c)
  {
    x[i] -= a;
    y[i] -= b;
    z[i] -= c;
  }

  //Scale values of a single element of the array's values.
  void Scale(const uint i, const real a, const real b, const real c)
  {
    x[i] *= a;
    y[i] *= b;
    z[i] *= c;
  }

  //Add values to a single array.
  void Add(const uint i, const real val)
  {
    x[i] += val;
    y[i] += val;
    z[i] += val;
  }

  //Sub to values within an array.
  void Sub(const uint i, const real val)
  {
    x[i] -= val;
    y[i] -= val;
    z[i] -= val;
  }

  //Multiply values within an array
  void Scale(const uint i, const real val)
  {
    x[i] *= val;
    y[i] *= val;
    z[i] *= val;
  }

  //Multiply values within an array
  void Scale(const uint i, XYZ const& val)
  {
    x[i] *= val.x;
    y[i] *= val.y;
    z[i] *= val.z;
  }

  //Add a constant transform x+a, y+b... to range of values.
  void AddRange(const uint start, const uint stop, XYZ const& val);

  //Subtract a constant transform x-a, y-b,... from a range of values.
  void SubRange(const uint start, const uint stop, XYZ const& val);

  //Multiply a constant transform x*a, y*b, ... from a range of values.
  void ScaleRange(const uint start, const uint stop, XYZ const& val);

  //Add a constant transform x+val, z+val,... to range of values.
  void AddRange(const uint start, const uint stop, const real val);

  //Subtract a constant transform x-val,y-val... from a range of values.
  void SubRange(const uint start, const uint stop, const real val);

  //Multiply a constant transform x*val,y*val... from a range of values.
  void ScaleRange(const uint start, const uint stop, const real val);

  //Copy one point
  void Set(XYZArray & dest, const uint srcIndex, const uint destIndex) const
  {
    dest.x[destIndex] = x[srcIndex];
    dest.y[destIndex] = y[srcIndex];
    dest.z[destIndex] = z[srcIndex];
  }

  //Add a constant transform x+a, y+b... to all values.
  void AddAll(XYZ const& val);

  //Subtract a constant transform x-a, y-b,... to all values.
  void SubAll(XYZ const& val);

  //Multiply a constant transform x*a, y*b, ... to all values.
  void ScaleAll(XYZ const& val);

  //Add a constant transform x+val, y+val... to all values.
  void AddAll(const real val);

  //Subtract a constant transform x-val, y-val,... to all values.
  void SubAll(const real val);

  //Multiply a constant transform x*val, y*val, ... to all values.
  void ScaleAll(const real val);

  ////////////////////////////////////////////////////////////////////

  //Copy range of points.
  void CopyRange(XYZArray & dest, const uint srcIndex, const uint destIndex,
                 const uint len) const;

  real * x, * y, * z;

protected:
  uint count;
  bool allocDone;
};


inline XYZArray::XYZArray(XYZArray const& other)
{
  count = other.count;
  x = new real[count];
  y = new real[count];
  z = new real[count];
  allocDone = true;
  memcpy(x, other.x, sizeof(real) * count);
  memcpy(y, other.y, sizeof(real) * count);
  memcpy(z, other.z, sizeof(real) * count);
}

//copy and swap assignment
inline XYZArray& XYZArray::operator=(XYZArray other)
{
  swap(*this, other);
  return *this;
}

inline void swap(XYZArray& a1, XYZArray& a2)
{
  using std::swap;
  swap(a1.x, a2.x);
  swap(a1.y, a2.y);
  swap(a1.z, a2.z);
  swap(a1.count, a2.count);
  swap(a1.allocDone, a2.allocDone);
}

inline void XYZArray::Uninit()
{
  if (x != NULL)
    delete[] x;
  if (y != NULL)
    delete[] y;
  if (z != NULL)
    delete[] z;
  allocDone = false;
}

inline void XYZArray::SetRange(const uint start, const uint stop,
                               const real a, const real b, const real c)
{
  SetRange(start, stop, XYZ(a, b, c));
}

inline void XYZArray::SetRange(const uint start, const uint stop,
                               XYZ const& val)
{
  for(uint i = start; i < stop ; ++i) {
    x[i] = val.x;
    y[i] = val.y;
    z[i] = val.z;
  }
}

inline void XYZArray::Init(const uint n)
{
  count = n;

  if (allocDone) {
    allocDone = false;
    if (x != NULL)
      delete[] x;
    if (y != NULL)
      delete[] y;
    if (z != NULL)
      delete[] z;
  }

  x = new real[n];
  y = new real[n];
  z = new real[n];

  allocDone = true;
}

inline void XYZArray::AddRange(const uint start, const uint stop,
                               XYZ const& val)
{
  for(uint i = start; i < stop ; ++i) {
    x[i] += val.x;
    y[i] += val.y;
    z[i] += val.z;
  }
}

inline void XYZArray::SubRange(const uint start, const uint stop,
                               XYZ const& val)
{
  for(uint i = start; i < stop ; ++i) {
    x[i] -= val.x;
    y[i] -= val.y;
    z[i] -= val.z;
  }
}

inline void XYZArray::ScaleRange(const uint start, const uint stop,
                                 XYZ const& val)
{
  for(uint i = start; i < stop ; ++i) {
    x[i] *= val.x;
    y[i] *= val.y;
    z[i] *= val.z;
  }
}

inline void XYZArray::AddRange(const uint start, const uint stop,
                               const real val)
{
  for(uint i = start; i < stop ; ++i) {
    x[i] += val;
    y[i] += val;
    z[i] += val;
  }
}

inline void XYZArray::SubRange(const uint start, const uint stop,
                               const real val)
{
  for(uint i = start; i < stop ; ++i) {
    x[i] -= val;
    y[i] -= val;
    z[i] -= val;
  }
}

inline void XYZArray::ScaleRange(const uint start, const uint stop,
                                 const real val)
{
  for(uint i = start; i < stop ; ++i) {
    x[i] *= val;
    y[i] *= val;
    z[i] *= val;
  }
}

//Add a constant transform x+a, y+b... to all values.
inline void XYZArray::AddAll(XYZ const& val)
{
  AddRange(0, count, val);
}

//Subtract a constant transform x-a, y-b,... to all values.
inline void XYZArray::SubAll(XYZ const& val)
{
  SubRange(0, count, val);
}

//Multiply a constant transform x*a, y*b, ... to all values.
inline void XYZArray::ScaleAll(XYZ const& val)
{
  ScaleRange(0, count, val);
}

//Add a constant transform x+val, y+val... to all values.
inline void XYZArray::AddAll(const real val)
{
  AddRange(0, count, val);
}

//Subtract a constant transform x-val, y-val,... to all values.
inline void XYZArray::SubAll(const real val)
{
  SubRange(0, count, val);
}

//Multiply a constant transform x*val, y*val, ... to all values.
inline void XYZArray::ScaleAll(const real val)
{
  ScaleRange(0, count, val);
}

//Copy range of points.
inline void XYZArray::CopyRange(XYZArray & dest, const uint srcIndex,
                                const uint destIndex, const uint len) const
{
#ifdef _OPENMP
  #pragma omp parallel default(shared)
#endif
  {
    memcpy(dest.x + destIndex, x + srcIndex, len * sizeof(real));
    memcpy(dest.y + destIndex, y + srcIndex, len * sizeof(real));
    memcpy(dest.z + destIndex, z + srcIndex, len * sizeof(real));
  }
}

inline real XYZArray::AdjointMatrix(XYZArray &Inv)
{
  Inv.x[0] = y[1] * z[2] - y[2] * z[1];
  Inv.y[0] = y[2] * z[0] - y[0] * z[2];
  Inv.z[0] = y[0] * z[1] - y[1] * z[0];

  Inv.x[1] = x[2] * z[1] - x[1] * z[2];
  Inv.y[1] = x[0] * z[2] - x[2] * z[0];
  Inv.z[1] = x[1] * z[0] - x[0] * z[1];

  Inv.x[2] = x[1] * y[2] - x[2] * y[1];
  Inv.y[2] = x[2] * y[0] - x[0] * y[2];
  Inv.z[2] = x[0] * y[1] - x[1] * y[0];

  real det = x[0] * Inv.x[0] + x[1] * Inv.y[0] + x[2] * Inv.z[0];
  return det;
}

#endif /*XYZ_ARRAY_H*/
