/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#include "OutConst.h" //For namespace spec;

namespace out {
const std::string ENERGY_TOTAL = "TOT_EN";
const std::string ENERGY_INTER = "EN_INTER";
const std::string ENERGY_LRC = "EN_LRC";
const std::string ENERGY_INTRA_B = "EN_INTRA(B)";
const std::string ENERGY_INTRA_NB = "EN_INTRA(NB)";
const std::string ENERGY_ELECT = "EN_ELECT";
const std::string ENERGY_REAL = "EN_REAL";
const std::string ENERGY_RECIP = "EN_RECIP";
const std::string VIRIAL_TOTAL = "TOTAL_VIR";
const std::string PRESSURE = "PRESSURE";
const std::string COMPRESSIBILITY = "COMPRESSIBILITY";
const std::string ENTHALPY = "ENTHALPY";
#if ENSEMBLE == GEMC
const std::string HEAT_OF_VAP = "HEAT_VAP";
const std::string VOLUME = "VOLUME";
#endif
#if ENSEMBLE == NPT
const std::string VOLUME = "VOLUME";
#endif
const std::string DENSITY = "TOT_DENS";
const std::string MOL_NUM = "TOT_MOL";
const std::string MOL_FRACTION = "MOLFRACT";
const std::string MOL_DENSITY = "MOLDENS";
const std::string SURF_TENSION = "SURF_TENSION";
} // namespace out
