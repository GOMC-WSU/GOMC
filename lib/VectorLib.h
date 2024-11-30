/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef VECTOR_LIB_H
#define VECTOR_LIB_H

#include <algorithm> //for vector copying
#include <sstream>   //for string parsing.
#include <string>    //for string splitting functions
#include <vector>    //main type this lib deals with

namespace vect {
inline std::vector<std::string> split(std::vector<std::string> &elems,
                                      std::string const &s, const char delim) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim))
    elems.push_back(item);
  return elems;
}
inline std::vector<std::string> split(const std::string &s, const char delim) {
  std::vector<std::string> elems;
  split(elems, s, delim);
  return elems;
}

// returns a pointer to an array containing a new copy of vec
template <typename T> T *transfer(const std::vector<T> &vec) {
  T *array = new T[vec.size()];
  std::copy(vec.begin(), vec.end(), array);
  return array;
}

// returns a pointer to an array containing a new copy of vec
template <typename T> T *TransferInto(T *array, const std::vector<T> &vec) {
  std::copy(vec.begin(), vec.end(), array);
  return array;
}
// overloaded because vector<bool> is a lie and copy may not be supported
inline bool *transfer(const std::vector<bool> &vec) {
  bool *array = new bool[vec.size()];
  for (size_t i = 0; i < vec.size(); ++i) {
    array[i] = vec[i];
  }
  return array;
}

inline void transfer(char *&arr, std::vector<std::string> const &v,
                     const uint len) {
  if (arr == NULL)
    arr = new char[v.size() * len];
  for (unsigned int i = 0; i < v.size(); i += len)
    for (unsigned int j = 0; j < len; j++)
      if (j < v[i].size())
        *(arr + i + j) = v[i][j];
      else
        *(arr + i + j) = ' ';
}
} // namespace vect

#endif /*VECTOR_LIB_H*/
