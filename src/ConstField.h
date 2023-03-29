/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
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
