/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef COORDINATE_SETUP_H
#define COORDINATE_SETUP_H

#include <string>
#include <vector>

struct CoordinateSetup {
  std::vector<double> partX;
  std::vector<double> partY;
  std::vector<double> partZ;
  std::vector<unsigned int> partBox;

  void Init(const std::string &pdbFilename);
  void SetCOM(const MolSetupData &molData);
}

#endif /*COORDINATE_SETUP_H*/
