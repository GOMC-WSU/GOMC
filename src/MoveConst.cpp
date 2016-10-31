/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.8
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "MoveConst.h" //Spec file.

std::vector<std::string> mv::MoveNames()
{
   std::vector<std::string> v;
#if ENSEMBLE == NVT
   v.push_back("Displace (Box 0)");
   v.push_back("Rotate (Box 0)");
   v.push_back("IntraSwap (Box 0)");
#elif ENSEMBLE == GCMC
   v.push_back("Displace (Box 0)");
   v.push_back("Rotate (Box 0)");
   v.push_back("Deletion (from Box 0)");
   v.push_back("Insertion (into Box 0)");
   v.push_back("IntraSwap (Box 0)");
#elif ENSEMBLE == GEMC
   v.push_back("Displace (Box 0)");
   v.push_back("Displace (Box 1)");
   v.push_back("Rotate (Box 0)");
   v.push_back("Rotate (Box 1)");
   v.push_back("Volume Transfer (Box 0 -> Box 1)");
   v.push_back("Volume Transfer (Box 1 -> Box 0)");
   v.push_back("Molecule Transfer (Box 0 -> Box 1)");
   v.push_back("Molecule Transfer (Box 1 -> Box 0)");
   v.push_back("IntraSwap (Box 0)");
   v.push_back("IntraSwap (Box 1)");
#endif
   return v;
}

std::vector<std::string> mv::ScaleMoveNames()
{
#if ENSEMBLE == GEMC
   std::vector<std::string> v(VOL_TRANSFER+1);
#else
   std::vector<std::string> v(ROTATE+1);
#endif
   v[DISPLACE] = "DISPLACE";
   v[ROTATE] = "ROTATE";
#if ENSEMBLE == GEMC
   v[VOL_TRANSFER] = "VOLUME";
#endif
   return v;
}
