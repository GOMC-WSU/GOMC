/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef FIXED_WIDTH_READER_H
#define FIXED_WIDTH_READER_H

#include "BasicTypes.h" //For "uint"
#include "ConstField.h" //For ConstField kind.
#include "Reader.h"     //Parent class
#include "StrLib.h"     //FromStr, StripWS
#include "StrStrmLib.h" //For stringstream operators

class FixedWidthReader : public Reader {
public:
  FixedWidthReader(std::string const &name, std::string const &alias,
                   const bool crit = true, const bool note = true)
      : Reader(name, alias, false, NULL, false, NULL, crit, note), line("") {}

  // this default constructor was defined for PDBSetup class to be able to
  // create objects without initializing right away.
  FixedWidthReader()
      : Reader("", "", false, NULL, false, NULL, true, true), line("") {}

  // Functions to get values from file, using fields.
  FixedWidthReader &Get(double &d, ConstField const &field) {
    sstrm::FromStr(d, Str(field));
    return *this;
  }
  FixedWidthReader &Get(float &f, ConstField const &field) {
    sstrm::FromStr(f, Str(field));
    return *this;
  }
  FixedWidthReader &Get(uint &ui, ConstField const &field) {
    sstrm::FromStr(ui, Str(field));
    return *this;
  }
  FixedWidthReader &Get(ulong &ul, ConstField const &field) {
    sstrm::FromStr(ul, Str(field));
    return *this;
  }
  FixedWidthReader &Get(std::string &s, ConstField const &field) {
    s = sstrm::StripWS(Str(field));
    return *this;
  }
  FixedWidthReader &Get(char &c, ConstField const &field) {
    c = line[field.START];
    return *this;
  }

  std::string GetLineCopy() const { return line; }

  // Gets line.
  bool Read(std::string &str, ConstField const &field) {
    if (GoodFileWData()) {
      std::getline(file, line);
      str = Str(field);
    }
    return GoodFileWData();
  }

  void SetData(std::string const &name, std::string const &alias) {
    fileName = name;
    fileAlias = alias;
    nameWAlias = fileAlias + ":  \t" + fileName;
  }

protected:
  std::string Str(ConstField const &field) {
    return line.substr(field.START, field.LENGTH);
  }
  std::string line;
};

#endif /*FIXED_WIDTH_READER_H*/
