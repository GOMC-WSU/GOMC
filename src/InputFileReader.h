/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
#ifndef INPUT_FILE_READER_H
#define INPUT_FILE_READER_H

#if GOMC_LIB_MPI
#include <mpi.h>
#endif

#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

class InputFileReader {
public:
  InputFileReader(const std::string& fileName);
  InputFileReader() : fs{} {}
  ~InputFileReader();
  bool readNextLine(std::vector<std::string> &str);
  void Open(const std::string& fileName);
  void Test(const std::string& fileName);

private:
  std::ifstream fs;
};

#endif /*INPUT_FILE_READER_H*/
