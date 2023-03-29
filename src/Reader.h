/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef READER_H
#define READER_H

#include <stdio.h> //for exit

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> //for skipping to end of line.
#include <string> //for skipping comparisons, file names etc.

#include "BasicTypes.h" //for uint, ulong
#include "StrLib.h"     //for str::compare(,) wrapper
#include "VectorLib.h"  //For vect::split() string sentence parser

inline std::ifstream &operator>>(std::ifstream &is, bool &val) {
  is >> std::noboolalpha >> val;
  if (is.fail()) {
    is.clear();
    is >> std::boolalpha >> val;
  }
  return is;
}

class Reader {
public:
  Reader(std::string const &name, std::string const &alias,
         bool useSkipW = false, std::string const *const skipW = NULL,
         bool useSkipC = false, std::string const *const skipC = NULL,
         const bool crit = true, const bool note = true)
      : firstVal(""), skipChars(""), skipWordsEnable(false),
        skipCharsEnable(false) {
    Init(name, alias, crit, note);
    if (useSkipW && skipW != NULL) {
      std::string s = *skipW;
      SetSkipWords(s);
    }
    if (useSkipC && skipC != NULL) {
      std::string s = *skipC;
      SetSkipChars(s);
    }
  }

  ~Reader(void) {
    if (isOpen)
      close();
  }

  // Set main class vars.
  void Init(std::string const &name, std::string const &alias, const bool crit,
            const bool note) {
    fileName = name;
    fileAlias = alias;
    critical = crit;
    notify = note;
    comp = true;
    isOpen = false;
    nameWAlias = fileAlias + ":  \t" + fileName;
  }

  // Open or close a file, with basic protections
  void open(void) {
    if (isOpen)
      return;
    file.open(fileName.c_str(), std::ios::in);
    CheckFileState(true, "...could not be opened.", "Reading from ");
  }

  void close(void) {
    if (!isOpen)
      return;
    file.close();
  }

  // Go to start of file
  void FileStartOver(void) {
    file.clear();
    file.seekg(0, std::ios::beg);
  }

  void SkipLine(void) {
    file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  void Skip(void) {
    if (GoodFileWData()) {
      std::string dummy;
      file >> dummy;
    }
  }

  bool Read(std::string &firstItem);

  // Make public to expose object.
  std::ifstream file;

protected:
  bool GoodFileWData(void) {
    return file.is_open() && file.good() && !file.eof();
  }

  void HandleError(std::string const &msg) {
    std::cerr << ((critical) ? "Error " : "Warning ") << nameWAlias << std::endl
              << msg << std::endl;
    if (critical)
      exit(1);
  }

  void HandleNote(std::string const &msg) {
    std::cout << msg << nameWAlias << std::endl;
  }
  void CheckFileState(const bool expected, std::string const &errMessage,
                      std::string const &note) {
    isOpen = GoodFileWData();
    if (isOpen == expected && notify)
      HandleNote(note);
    else if (isOpen != expected)
      HandleError(errMessage);
  }

  // Skip unwanted items -- comments, etc.
  void SetSkipWords(std::string const &skipW) {
    skipWordsEnable = (skipW.length() > 0);
    vect::split(skipWords, skipW, ' ');
  }

  void SetSkipChars(std::string const &skipC) {
    skipCharsEnable = (skipC.length() > 0);
    skipChars = skipC;
  }

  bool CheckSkipChars(std::string const &s) {
    bool skip = false;
    if (skipCharsEnable && !skip)
      for (uint c = 0; c < skipChars.length() && !skip; c++)
        skip |= (s[0] == skipChars[c]);
    return skip;
  }

  bool CheckSkipWords(std::string const &s) {
    bool skip = false;
    if (skipWordsEnable)
      for (uint w = 0; w < skipWords.size() && !skip; w++)
        skip |= str::compare(s, skipWords[w]);
    return skip;
  }

  std::string firstVal;

  // For skipping unwanted items
  std::string skipChars;
  std::vector<std::string> skipWords;
  bool skipWordsEnable, skipCharsEnable;

  std::string fileName, fileAlias, nameWAlias;
  bool critical, notify, isOpen, comp;
};

#endif /*READER_H*/
