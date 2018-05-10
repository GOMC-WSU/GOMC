/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.30
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "OutConst.h" //For namespace spec;

namespace out
{
const std::string ENERGY_TOTAL = "TOT_EN";
const std::string ENERGY_INTER = "EN_INTER";
const std::string ENERGY_TC = "EN_TC";
const std::string ENERGY_INTRA_B  = "EN_INTRA(B)";
const std::string ENERGY_INTRA_NB = "EN_INTRA(NB)";
const std::string ENERGY_ELECT = "EN_ELECT";
const std::string ENERGY_REAL = "EN_REAL";
const std::string ENERGY_RECIP = "EN_RECIP";
const std::string VIRIAL_TOTAL = "TOTAL_VIR";
const std::string PRESSURE = "PRESSURE";
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
const std::string SURF_TENSION = "SURF_TENSION";
}
