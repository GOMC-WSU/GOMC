#ifndef CONST_FIELD_H
#define CONST_FIELD_H

#include "BasicTypes.h" //For "uint"

struct ConstField
{
  ConstField(const uint st, const uint len) : START(st), LENGTH(len) {}
  const uint START;
  const uint LENGTH;
};

#endif /*CONST_FIELD_H*/
