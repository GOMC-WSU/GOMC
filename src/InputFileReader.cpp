/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "InputFileReader.h"
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

/*****************************
**
**	Author: Younes Nejahi
**
*****************************/

InputFileReader::InputFileReader(std::string inputFileName)
{
  fs.open(inputFileName.c_str(), std::fstream::in);

}

InputFileReader::InputFileReader()
{
}

InputFileReader::~InputFileReader()
{
  //fs.close();
}

void InputFileReader::Open(std::string inputFileName)
{
  fs.open(inputFileName.c_str(), std::fstream::in);
  if(!fs.is_open()) {
    std::cout << "Cannot open input file!" << std::endl;
    exit(EXIT_FAILURE);
  }
}


/*
** Read one line from the input and seperate it into strings.
** ***** str -> components of each line will be stored here
*/
bool InputFileReader::readNextLine(std::vector<std::string> & str)
{
  std::string line;
  do {
    if (fs.eof() || fs.bad() || fs.fail()) {
      return false;
    }
    std::getline(fs, line);
    if(!line.size())
      line = "#";
  } while (line[0] == '#' || line[0] == '\0');

  std::istringstream iss(line);
  copy(std::istream_iterator<std::string>(iss),
       std::istream_iterator<std::string>(),
       back_inserter(str));
  return true;
}
