/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.11
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#pragma once
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class InputFileReader
{
private:
  fstream fs;
  vector<string> & split(const string &s, char delim, vector<string> &elems);

public:
  bool readNextLine(vector<string> & str);
  void Open(string fileName);
  InputFileReader(string fileName);
  InputFileReader(void);
  ~InputFileReader();
};
