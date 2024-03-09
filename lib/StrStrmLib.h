/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef STR_STRM_LIB_H
#define STR_STRM_LIB_H

#include <sstream> //The obj. this set of util. functions deals primarily with.
#include <string>  //For functions involving strings

#include "BasicTypes.h" //For ulong/uint
#include "ConstField.h" //For const field type.
#include "StrLib.h"     //For str::GetDat

namespace align {
static const char LEFT = 'l';
static const char RIGHT = 'r';
} // namespace align

namespace sstrm {
inline std::string StripWS(std::string const &str) {
  std::string val;
  std::stringstream strm(str);
  strm >> val;
  return val;
}

// Count words, ignoring white space.
// See:
// http://stackoverflow.com/questions/3672234/
//     c-function-to-count-all-the-words-in-a-string
inline uint CountWords(std::string &str) {
  std::stringstream strm;
  // sneaky way to use the string as the buffer to avoid copy.
  strm.rdbuf()->pubsetbuf(&str[0], str.length());
  return std::distance(std::istream_iterator<std::string>(strm),
                       std::istream_iterator<std::string>());
}

template <typename T> inline T FromStr(T &var, std::string const &str) {
  std::stringstream strm(str);
  strm >> var;
  return var;
}

inline bool FromStr(bool &b, std::string const &str) {
  std::stringstream strm(str);
  if (strm.str().length() == 1)
    strm << std::noboolalpha;
  strm >> b;
  return b;
}
inline char *&FromStr(char *&c, std::string const &str) {
  for (unsigned int i = 0; i < str.length(); i++)
    c[i] = str[i];
  return c;
}

struct Converter {
  std::stringstream strm;

  Converter() { AutoFmt(); }

  // No width or precision set, by default.
  Converter &AutoFmt() {
    Align(align::RIGHT).Fixed();
    return *this;
  }

  Converter &Align(const char align) {
    if (align == align::LEFT)
      strm.setf(std::ios_base::left, std::ios_base::adjustfield);
    else // if (align == align::RIGHT)
      strm.setf(std::ios_base::right, std::ios_base::adjustfield);
    return *this;
  }
  Converter &Fixed() {
    strm.setf(std::ios_base::fixed, std::ios_base::floatfield);
    return *this;
  }

  Converter &Precision(const int precision) {
    if (precision != -1)
      strm.precision(precision);
    return *this;
  }

  Converter &Width(const int len) {
    if (len != -1)
      strm.width(len);
    return *this;
  }

  template <typename T> Converter &operator<<(const T &var) {
    strm << var;
    return *this;
  }

  Converter &operator<<(const bool b) {
    if (strm.precision() == 1)
      strm << std::noboolalpha;
    strm << b;
    return *this;
  }
  Converter &operator>>(std::string &str) {
    str = strm.str();
    Flush();
    return *this;
  }

  template <typename T>
  void Replace(std::string &str, const T &var, ConstField const &field) {
    std::string temp;
    Width(field.LENGTH);
    *this << var;
    *this >> temp;
    str.replace(field.START, field.LENGTH, temp);
  }

  void Replace(std::string &str, char const *const cStr,
               ConstField const &field) {
    for (uint c = 0; c < field.LENGTH; ++c)
      str[c + field.START] = cStr[c];
  }

  Converter &Flush() {
    strm.clear();
    strm.str("");
    AutoFmt();
    return *this;
  }
};

inline void Flush(std::stringstream &strm) {
  strm.clear();
  strm.str("");
}
} // namespace sstrm

#endif /*STR_STRM_LIB_H*/
