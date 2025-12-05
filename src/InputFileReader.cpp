/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#include "InputFileReader.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

/*****************************
**
**    Author: Younes Nejahi
**
*****************************/

InputFileReader::InputFileReader(std::string inputFileName) {
  fs.open(inputFileName.c_str(), std::fstream::in);
}

InputFileReader::InputFileReader() {}

InputFileReader::~InputFileReader() {
  // fs.close();
}

void InputFileReader::Open(std::string inputFileName) {
  fs.open(inputFileName.c_str(), std::fstream::in);
  if (!fs.is_open()) {
    std::cout << "Cannot open input file!" << std::endl;
    exit(EXIT_FAILURE);
  }
}

/*
** Read one line from the input and seperate it into strings.
** ***** str -> components of each line will be stored here
*/
bool InputFileReader::readNextLine(std::vector<std::string> &str) {
  std::string line;
  do {
    if (fs.eof() || fs.bad() || fs.fail()) {
      return false;
    }
    std::getline(fs, line);
    if (!line.size())
      line = "#";
  } while (line[0] == '#' || line[0] == '\0');

  std::istringstream iss(line);
  copy(std::istream_iterator<std::string>(iss),
       std::istream_iterator<std::string>(), back_inserter(str));
  return true;
}
