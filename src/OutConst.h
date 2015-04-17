/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) BETA 0.97 (Serial version)
Copyright (C) 2015  GOMC Group

A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef OUT_CONST_H
#define OUT_CONST_H

#include <string>

#include "EnsemblePreprocessor.h" //For VARIABLE_x, etc.
#include "../lib/BasicTypes.h" //For uint

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
   extern const std::string VIRIAL_TOTAL;
   static const uint VIRIAL_TOTAL_IDX = 5;
   extern const std::string VIRIAL_INTER;
   static const uint VIRIAL_INTER_IDX = 6;
   extern const std::string VIRIAL_TC;
   static const uint VIRIAL_TC_IDX = 7;
   extern const std::string PRESSURE;
   static const uint PRESSURE_IDX = 8;
#if ENSEMBLE == NVT || ENSEMBLE == GCMC
   static const uint TOTAL_SINGLE = 9;
#elif ENSEMBLE == GEMC
   extern const std::string VOLUME;
   static const uint VOLUME_IDX = 9;
   extern const std::string HEAT_OF_VAP;
   static const uint HEAT_OF_VAP_IDX = 10;
   static const uint TOTAL_SINGLE = 11; 
#endif

   //MULTI
#if ENSEMBLE == GCMC || ENSEMBLE == GEMC
   extern const std::string MOL_NUM;
   static const uint MOL_NUM_IDX = 0;
   extern const std::string DENSITY;
   static const uint DENSITY_IDX = 1;
   extern const std::string MOL_FRACTION;
   static const uint MOL_FRACTION_IDX = 2;
   static const uint TOTAL_K = 3;
#elif ENSEMBLE == NVT
   static const uint TOTAL_K = 0;
#endif
}

#endif /*OUT_CONST_H*/

