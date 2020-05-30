/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#include <iostream>
#include <fstream>
#include <vector>

class InputFileReader
{
private:
  std::fstream fs;
  std::vector<string> & split(const std::string &s, char delim, std::vector<std::string> &elems);

public:
  bool readNextLine(std::vector<std::string> & str);
  void Open(std::string fileName);
  InputFileReader(std::string fileName);
  InputFileReader(void);
  ~InputFileReader();
};
