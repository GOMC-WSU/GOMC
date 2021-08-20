#include <vector>
#include "BasicTypes.h"
#include "Endian.h"
#include <stdint.h>
#include <iostream>

namespace
{
union dbl_input_union {
  char bin_value[8];
  double dbl_value;
};

union uint64_input_union {
  char bin_value[8];
  uint64_t uint_value;
};

union uint32_input_union {
  char bin_value[4];
  uint32_t uint_value;
};

union int8_input_union {
  char bin_value[1];
  int8_t int_value;
};
}

class CheckpointUtilities {

static void readVector3DDouble(FILE * stream, std::vector< std::vector< std::vector <double> > > & data);
static void readVector3DUint(FILE * stream, std::vector< std::vector< std::vector <uint> > > & data);
static void readVector2DUint(FILE * stream, std::vector< std::vector< uint > > & data);
static void readVector1DDouble(FILE * stream, std::vector< double > & data);
static double read_double_binary(FILE * stream);
static int8_t read_uint8_binary(FILE * stream);
static uint32_t read_uint32_binary(FILE * stream);
static uint64_t read_uint64_binary(FILE * stream);

};
