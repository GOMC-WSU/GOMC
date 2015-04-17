/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
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

