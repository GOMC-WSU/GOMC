/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "Reader.h"

#include <iostream> //for cout
#include <string>   //value read

bool Reader::Read(std::string &firstItem) {
  while (GoodFileWData() && (file >> firstVal))
    if (CheckSkipChars(firstVal) || CheckSkipWords(firstVal))
      SkipLine();
    else
      break;
  // commented out debug because it only tells us we have successfully
  // ignored comments and prints a lot of text to do so
  /*
  #ifndef NDEBUG
  std::streampos pos = file.tellg();
  std::string currLine;
  if ( std::getline(file, currLine) )
  {
     if (file.eof())
   file.clear();
     std::cout << firstVal << currLine << std::endl;
     file.seekg(-(currLine.size()+1), std::ios_base::cur);
  }
  #endif
  */
  firstItem = firstVal;
  return GoodFileWData();
}
