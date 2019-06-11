/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef FIXED_WIDTH_READER_H
#define FIXED_WIDTH_READER_H

#include "BasicTypes.h" //For "uint"
#include "Reader.h" //Parent class
#include "ConstField.h" //For ConstField kind.
#include "StrLib.h" //FromStr, StripWS
#include "StrStrmLib.h" //For stringstream operators


class FixedWidthReader : public Reader
{
public:
  FixedWidthReader(std::string const& name, std::string const& alias,
                   const bool crit = true, const bool note = true):
    Reader(name, alias, false, NULL, false, NULL, crit, note), line("") {}

  //Functions to get values from file, using fields.
  FixedWidthReader & Get(real & d, ConstField const& field)
  {
    sstrm::FromStr(d, Str(field));
    return *this;
  }
  FixedWidthReader & Get(uint & ui, ConstField const& field)
  {
    sstrm::FromStr(ui, Str(field));
    return *this;
  }
  FixedWidthReader & Get(ulong & ul, ConstField const& field)
  {
    sstrm::FromStr(ul, Str(field));
    return *this;
  }
  FixedWidthReader & Get(std::string & s, ConstField const& field)
  {
    s = sstrm::StripWS(Str(field));
    return *this;
  }
  FixedWidthReader & Get(char & c, ConstField const& field)
  {
    c = line[field.START];
    return *this;
  }

  std::string GetLineCopy() const
  {
    return line;
  }

  //Gets line.
  bool Read(std::string & str, ConstField const& field)
  {
    if (GoodFileWData()) {
      std::getline(file, line);
      str = Str(field);
#ifndef NDEBUG
      //big ol' waste of lines
      //std::cout << line << std::endl;
#endif
    }
    return GoodFileWData();
  }

protected:
  std::string Str(ConstField const& field)
  {
    return line.substr(field.START, field.LENGTH);
  }
  std::string line;
};

#endif /*FIXED_WIDTH_READER_H*/
