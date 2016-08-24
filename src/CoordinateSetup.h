#ifndef COORDINATESETUP_H
#define COORDINATESETUP_H

#include <vector>
#include <string>

struct CoordinateSetup
{
  std::vector<double> partX;
  std::vector<double> partY;
  std::vector<double> partZ;
  std::vector<unsigned int> partBox;

  void Init(const std::string& pdbFilename);
  void SetCOM(const MolSetupData& molData);
}

#endif
