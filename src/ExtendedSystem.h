/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.70
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef EXTENDED_SYSTEM
#define EXTENDED_SYSTEM

#include "InputAbstracts.h" //For FWReadableBase
#include "BasicTypes.h" //For uint
#include "EnsemblePreprocessor.h" //For BOX_TOTAL, etc.
#include "XYZArray.h" //For box dimensions.
#include "PDBSetup.h"
#include "GeomLib.h" // for PI and Dot product
#include <vector>

namespace config_setup
{
  struct RestartSettings;
  struct InFiles;
  struct Input;
}

namespace pdb_setup
{
  struct Cryst1;
}


class ExtendedSystem  {
  public:
  ExtendedSystem();
  ~ExtendedSystem() {};
  void Init(PDBSetup &pdb, config_setup::Input inputFiles);
  private:
    void ReadExtendedSystem(const char *filename, const int box);
    void UpdateCellBasis(PDBSetup &pdb, const int box);
    ulong firstStep;
    XYZ center[BOX_TOTAL];
    XYZArray cellBasis[BOX_TOTAL];
    XYZArray axis;
    double cosAngle[BOX_TOTAL][3];
    double cellAngle[BOX_TOTAL][3];
    bool hasCellBasis[BOX_TOTAL];
};


#endif /*EXTENDED_SYSTEM*/
