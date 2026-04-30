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

InputFileReader::InputFileReader(const std::string& inputFileName) {
  fs.open(inputFileName.c_str(), std::fstream::in);
}

InputFileReader::~InputFileReader() {
  if (fs.is_open()) {
    fs.close();
  }
}

void InputFileReader::Open(const std::string& inputFileName) {
  fs.open(inputFileName.c_str(), std::fstream::in);
  if (!fs.is_open()) {
    std::cout << "Cannot open input file!" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// Used to verify the existence of a file, ensure that it can be opened for
// input, and generate an appropriate error message if it can't be. If the file
// is opened, then we close the file, since it's not being used yet, just
// verified.
void InputFileReader::Test(const std::string& inputFileName) {
  fs.open(inputFileName.c_str(), std::fstream::in);

  // Check if file has been opened. If not, output an error message and exit
  if (!fs.is_open()) {
    // First check if the file exists
    if (!std::filesystem::exists(inputFileName.c_str())) {
#if GOMC_LIB_MPI
      MPI_Finalize();
#endif
      int error_num = static_cast<int>(std::errc::no_such_file_or_directory);
      std::error_code ec =
          std::make_error_code(std::errc::no_such_file_or_directory);
      std::cout << "Problem opening " << inputFileName << ": " << ec.message()
                << "\n";
      exit(error_num);
    }
    // Since the file exists, check if there is a problem with the directory
    // or file permissions
    std::ifstream tempFile(inputFileName.c_str());
    if (!tempFile.is_open()) {
      int error_num = static_cast<int>(std::errc::permission_denied);
      std::error_code ec = std::make_error_code(std::errc::permission_denied);
      std::cout << "Problem opening " << inputFileName << ": File found but "
                << ec.message() << "\n";
      exit(error_num);
    }
    // Failed for some unknown reason. Shouldn't reach here but need to
    // terminate if we do
    std::cout << "Problem opening " << inputFileName << ": Unexpected error!\n";
    exit(2);
  }

  // Close file to be reopened later when the file is ready to be used
  fs.close();
}

/*
** Read one line from the input and separate it into strings.
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
