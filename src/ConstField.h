/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.60
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef CONST_FIELD_H
#define CONST_FIELD_H

#include "BasicTypes.h" //For "uint"

struct ConstField {
  ConstField(const uint st, const uint len) : START(st), LENGTH(len) {}
  const uint START;
  const uint LENGTH;
};

#endif /*CONST_FIELD_H*/
