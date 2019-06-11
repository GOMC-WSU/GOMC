/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef COORDINATESETUP_H
#define COORDINATESETUP_H

#include <vector>
#include <string>

struct CoordinateSetup {
  std::vector<real> partX;
  std::vector<real> partY;
  std::vector<real> partZ;
  std::vector<unsigned int> partBox;

  void Init(const std::string& pdbFilename);
  void SetCOM(const MolSetupData& molData);
}

#endif
