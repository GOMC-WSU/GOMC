/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef PDB_CONST_H
#define PDB_CONST_H

#include <string> //For constant labels

#include "ConstField.h" //For ConstField kind.
#include "StrStrmLib.h"

namespace pdb_entry
{
static const uint LINE_WIDTH = 80;
namespace end
{
static const std::string STR = "END";
static const ConstField POS(0, 3);
}
namespace label
{
static const std::string REMARK = "REMARK";
static const std::string CRYST1 = "CRYST1";
static const std::string ATOM = "ATOM  ";
static const ConstField POS(0, 6);
}
namespace remark
{
namespace field
{
static const uint NON_SCALE_TAGS = 2;
namespace rem_num
{
static const ConstField POS(7, 3);
}
namespace name
{
static const std::string STR_GOMC = "GOMC";
static const std::string STR_STEP = "STEP";
static const ConstField POS(11, 15);
}
namespace data
{
static const ConstField POS(30, 15);
}
}
}
namespace cryst1
{
namespace field
{
namespace x
{
static const ConstField POS(6, 9);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace y
{
static const ConstField POS(15, 9);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace z
{
static const ConstField POS(24, 9);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace ang_alpha
{
static const ConstField POS(33, 7);
static const real DEFAULT = 90.0;
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
}
namespace ang_beta
{
static const ConstField POS(40, 7);
static const real DEFAULT = 90.0;
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
}
namespace ang_gamma
{
static const ConstField POS(47, 7);
static const real DEFAULT = 90.0;
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
}
namespace space
{
static const ConstField POS(55, 11);
static const std::string DEFAULT = "P 1";
}
namespace zvalue
{
static const ConstField POS(66, 4);
static const std::string DEFAULT = "1";
}
namespace dis
{
static const ConstField POS(15, 10);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace rot
{
static const ConstField POS(25, 10);
static const uint PRECISION = 5;
static const char ALIGN = align::RIGHT;
}
namespace vol
{
static const ConstField POS(35, 12);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace frameNum
{
static const ConstField POS(15, 10);
static const uint PRECISION = 0;
static const char ALIGN = align::RIGHT;
}
namespace stepsNum
{
static const ConstField POS(50, 15);
static const uint PRECISION = 0;
static const char ALIGN = align::RIGHT;
}
}
}
namespace atom
{
namespace field
{
namespace atom_num
{
static const ConstField POS(6, 5);
static const char ALIGN = align::RIGHT;
}
namespace alias
{
static const ConstField POS(12, 4);
}
namespace alt_location
{
static const ConstField POS(16, 1);
}
namespace res_name
{
static const ConstField POS(17, 4);
}
namespace chain
{
static const ConstField POS(21, 1);
}
namespace res_num
{
static const ConstField POS(22, 4);
static const char ALIGN = align::RIGHT;
}
namespace icode
{
static const ConstField POS(26, 1);
}
namespace x
{
static const ConstField POS(30, 8);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace y
{
static const ConstField POS(38, 8);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace z
{
static const ConstField POS(46, 8);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
}
namespace occupancy
{
static const ConstField POS(54, 6);
static const std::string BOX[] = {"0.00", "1.00"};
static const char ALIGN = align::RIGHT;
}
namespace beta
{
static const ConstField POS(60, 6);
static const real DEFAULT = 0.00;
static const std::string FIX[] = {"0.00", "1.00", "2.00"};
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
}
namespace el
{
static const ConstField POS(76, 2);
}
namespace charge
{
static const ConstField POS(78, 2);
}
}
}
}

#endif /*PDB_CONST_H*/
