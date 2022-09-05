/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef PDB_CONST_H
#define PDB_CONST_H

#include <string> //For constant labels

#include "ConstField.h" //For ConstField kind.
#include "StrStrmLib.h"

namespace pdb_entry {
static const uint LINE_WIDTH = 80;
namespace end {
static const std::string STR = "END";
static const ConstField POS(0, 3);
} // namespace end
namespace label {
static const std::string REMARK = "REMARK";
static const std::string CRYST1 = "CRYST1";
static const std::string ATOM = "ATOM  ";
static const std::string HETATM = "HETATM";
static const ConstField POS(0, 6);
} // namespace label
namespace remark {
namespace field {
static const uint NON_SCALE_TAGS = 2;
namespace rem_num {
static const ConstField POS(7, 3);
}
namespace name {
static const std::string STR_GOMC = "GOMC";
static const std::string STR_STEP = "STEP";
static const ConstField POS(11, 15);
} // namespace name
namespace data {
static const ConstField POS(30, 15);
}
} // namespace field
} // namespace remark
namespace cryst1 {
namespace field {
namespace x {
static const ConstField POS(6, 9);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace x
namespace y {
static const ConstField POS(15, 9);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace y
namespace z {
static const ConstField POS(24, 9);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace z
namespace ang_alpha {
static const ConstField POS(33, 7);
static const double DEFAULT = 90.0;
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
} // namespace ang_alpha
namespace ang_beta {
static const ConstField POS(40, 7);
static const double DEFAULT = 90.0;
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
} // namespace ang_beta
namespace ang_gamma {
static const ConstField POS(47, 7);
static const double DEFAULT = 90.0;
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
} // namespace ang_gamma
namespace space {
static const ConstField POS(55, 11);
static const std::string DEFAULT = "P 1";
} // namespace space
namespace zvalue {
static const ConstField POS(66, 4);
static const std::string DEFAULT = "1";
} // namespace zvalue
namespace dis {
static const ConstField POS(15, 10);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace dis
namespace rot {
static const ConstField POS(25, 10);
static const uint PRECISION = 5;
static const char ALIGN = align::RIGHT;
} // namespace rot
namespace vol {
static const ConstField POS(35, 12);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace vol
namespace frameNum {
static const ConstField POS(15, 10);
static const uint PRECISION = 0;
static const char ALIGN = align::RIGHT;
} // namespace frameNum
namespace stepsNum {
static const ConstField POS(50, 15);
static const uint PRECISION = 0;
static const char ALIGN = align::RIGHT;
} // namespace stepsNum
} // namespace field
} // namespace cryst1
namespace atom {
namespace field {
namespace atom_num {
static const ConstField POS(6, 5);
static const char ALIGN = align::RIGHT;
} // namespace atom_num
namespace alias {
static const ConstField POS(12, 4);
}
namespace alt_location {
static const ConstField POS(16, 1);
}
namespace res_name {
static const ConstField POS(17, 4);
}
namespace chain {
static const ConstField POS(21, 1);
}
namespace res_num {
static const ConstField POS(22, 4);
static const char ALIGN = align::RIGHT;
} // namespace res_num
namespace icode {
static const ConstField POS(26, 1);
}
namespace x {
static const ConstField POS(30, 8);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace x
namespace y {
static const ConstField POS(38, 8);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace y
namespace z {
static const ConstField POS(46, 8);
static const uint PRECISION = 3;
static const char ALIGN = align::RIGHT;
} // namespace z
namespace occupancy {
static const ConstField POS(54, 6);
static const std::string BOX[] = {"0.00", "1.00"};
static const char ALIGN = align::RIGHT;
} // namespace occupancy
namespace beta {
static const ConstField POS(60, 6);
static const double DEFAULT = 0.00;
static const std::string FIX[] = {"0.00", "1.00", "2.00"};
static const uint PRECISION = 2;
static const char ALIGN = align::RIGHT;
} // namespace beta
namespace el {
static const ConstField POS(76, 2);
}
namespace charge {
static const ConstField POS(78, 2);
}
} // namespace field
} // namespace atom
} // namespace pdb_entry

#endif /*PDB_CONST_H*/
