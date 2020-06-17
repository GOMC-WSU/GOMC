/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include <string> //value read
#include <iostream> //for cout

#include "Reader.h"

bool Reader::Read(std::string & firstItem)
{
  while(GoodFileWData() && (file >> firstVal) )
    if ( CheckSkipChars(firstVal) || CheckSkipWords(firstVal) )
      SkipLine();
    else
      break;
  //commented out debug because it only tells us we have successfully
  //ignored comments and prints a lot of text to do so
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
