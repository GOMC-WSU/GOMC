/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef STR_STRM_LIB_H
#define STR_STRM_LIB_H

#include <sstream> //The obj. this set of util. functions deals primarily with.
#include <string> //For functions involving strings

#include "BasicTypes.h" //For ulong/uint
#include "StrLib.h" //For str::GetDat
#include "ConstField.h" //For const field type.

namespace align
{
static const char LEFT = 'l';
static const char RIGHT = 'r';
}

namespace sstrm
{
inline std::string StripWS(std::string const& str)
{
  std::string val;
  std::stringstream strm(str);
  strm >> val;
  return val;
}

// Count words, ignoring white space.
// See:
// http://stackoverflow.com/questions/3672234/
//     c-function-to-count-all-the-words-in-a-string
inline uint CountWords(std::string & str)
{
  std::stringstream strm;
  // sneaky way to use the string as the buffer to avoid copy.
  strm.rdbuf()->pubsetbuf(&str[0], str.length() );
  return std::distance(std::istream_iterator<std::string>(strm),
                       std::istream_iterator<std::string>());
}

inline real FromStr(real & d, std::string const& str)
{
  std::stringstream strm(str);
  strm >> d;
  return d;
}
inline ulong FromStr(ulong & ul, std::string const& str)
{
  std::stringstream strm(str);
  strm >> ul;
  return ul;
}
inline long FromStr(long & l, std::string const& str)
{
  std::stringstream strm(str);
  strm >> l;
  return l;
}
inline int FromStr(int & i, std::string const& str)
{
  std::stringstream strm(str);
  strm >> i;
  return i;
}
inline uint FromStr(uint & ui, std::string const& str)
{
  std::stringstream strm(str);
  strm >> ui;
  return ui;
}
inline bool FromStr(bool & b, std::string const& str)
{
  std::stringstream strm(str);
  if (strm.str().length() == 1) strm << std::noboolalpha;
  strm >> b;
  return b;
}
inline char *& FromStr(char *& c, std::string const& str)
{
  for (unsigned int i = 0; i < str.length(); i++)
    c[i] = str[i];
  return c;
}


struct Converter {
  std::stringstream strm;

  Converter()
  {
    AutoFmt();
  }

  //No width or precision set, by default.
  Converter & AutoFmt()
  {
    Align(align::RIGHT).Fixed();
    return *this;
  }

  Converter & Align(const char align)
  {
    if (align == align::LEFT)
      strm.setf(std::ios_base::left, std::ios_base::adjustfield);
    else //if (align == align::RIGHT)
      strm.setf(std::ios_base::right, std::ios_base::adjustfield);
    return *this;
  }
  Converter & Fixed()
  {
    strm.setf(std::ios_base::fixed, std::ios_base::floatfield);
    return *this;
  }

  Converter & Precision(const int precision)
  {
    if (precision != -1)
      strm.precision(precision);
    return *this;
  }

  Converter & Width(const int len)
  {
    if (len != -1)
      strm.width(len);
    return *this;
  }

  Converter & operator<<(const real d)
  {
    strm << d;
    return *this;
  }
  Converter & operator<<(const long l)
  {
    strm << l;
    return *this;
  }
  Converter & operator<<(const int i)
  {
    strm << i;
    return *this;
  }
  Converter & operator<<(const uint ui)
  {
    strm << ui;
    return *this;
  }
  Converter & operator<<(std::string const& s)
  {
    strm << s;
    return * this;
  }
  Converter & operator<<(const bool b)
  {
    if (strm.precision() == 1)
      strm << std::noboolalpha;
    strm << b;
    return *this;
  }
  Converter & operator>>(std::string & str)
  {
    str = strm.str();
    Flush();
    return *this;
  }

  void Replace(std::string & str, const real d, ConstField const& field)
  {
    std::string temp;
    Width(field.LENGTH);
    *this << d;
    *this >> temp;
    str.replace(field.START, field.LENGTH, temp);
  }
  void Replace(std::string & str, const long l, ConstField const& field)
  {
    std::string temp;
    Width(field.LENGTH);
    *this << l;
    *this >> temp;
    str.replace(field.START, field.LENGTH, temp);
  }
  void Replace(std::string & str, const int i, ConstField const& field)
  {
    std::string temp;
    Width(field.LENGTH);
    *this << i;
    *this >> temp;
    str.replace(field.START, field.LENGTH, temp);
  }
  void Replace(std::string & str, const uint ui, ConstField const& field)
  {
    std::string temp;
    Width(field.LENGTH);
    *this << ui;
    *this >> temp;
    str.replace(field.START, field.LENGTH, temp);
  }
  void Replace(std::string & str, const bool b, ConstField const& field)
  {
    std::string temp;
    Width(field.LENGTH);
    *this << b;
    *this >> temp;
    str.replace(field.START, field.LENGTH, temp);
  }
  void Replace(std::string & str, char const*const cStr,
               ConstField const& field)
  {
    for (uint c = 0; c < field.LENGTH; ++c)
      str[c + field.START] = cStr[c];
  }
  void Replace(std::string & str, std::string const& subStr,
               ConstField const& field)
  {
    std::string temp;
    Width(field.LENGTH);
    *this << subStr;
    *this >> temp;
    str.replace(field.START, field.LENGTH, temp);
  }

  Converter & Flush()
  {
    strm.clear();
    strm.str("");
    AutoFmt();
    return *this;
  }
};

inline void Flush(std::stringstream & strm)
{
  strm.clear();
  strm.str("");
}
}

#endif /*STR_STRM_LIB_H*/
