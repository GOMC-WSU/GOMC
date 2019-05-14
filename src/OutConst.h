/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef OUT_CONST_H
#define OUT_CONST_H

#include <string>

#include "EnsemblePreprocessor.h" //For VARIABLE_x, etc.
#include "BasicTypes.h" //For uint

namespace out
{
extern const std::string ENERGY_TOTAL;
static const uint ENERGY_TOTAL_IDX = 0;
extern const std::string ENERGY_INTER;
static const uint ENERGY_INTER_IDX = 1;
extern const std::string ENERGY_TC;
static const uint ENERGY_TC_IDX = 2;
extern const std::string ENERGY_INTRA_B;
static const uint ENERGY_INTRA_B_IDX = 3;
extern const std::string ENERGY_INTRA_NB;
static const uint ENERGY_INTRA_NB_IDX = 4;
extern const std::string ENERGY_ELECT;
static const uint ENERGY_ELECT_IDX = 5;
extern const std::string ENERGY_REAL;
static const uint ENERGY_REAL_IDX = 6;
extern const std::string ENERGY_RECIP;
static const uint ENERGY_RECIP_IDX = 7;
extern const std::string VIRIAL_TOTAL;
static const uint VIRIAL_TOTAL_IDX = 8;
extern const std::string PRESSURE;
static const uint PRESSURE_IDX = 9;
extern const std::string MOL_NUM;
static const uint MOL_NUM_IDX = 10;
extern const std::string DENSITY;
static const uint DENSITY_IDX = 11;
extern const std::string SURF_TENSION;
static const uint SURF_TENSION_IDX = 12;
#if ENSEMBLE == NVT || ENSEMBLE == GCMC
static const uint TOTAL_SINGLE = 13;
#elif ENSEMBLE == NPT
extern const std::string VOLUME;
static const uint VOLUME_IDX = 13;
static const uint TOTAL_SINGLE = 14;
#else
extern const std::string VOLUME;
extern const std::string HEAT_OF_VAP;
static const uint VOLUME_IDX = 13;
static const uint HEAT_OF_VAP_IDX = 14;
static const uint TOTAL_SINGLE = 15;
#endif

//MULTI
extern const std::string MOL_FRACTION;
extern const std::string MOL_DENSITY;
static const uint MOL_FRACTION_IDX = 0;
static const uint MOL_DENSITY_IDX = 1;
static const uint TOTAL_K = 2;
}

#endif /*OUT_CONST_H*/
