/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.70 (Serial version)
Copyright (C) 2015  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#include "OutConst.h" //For namespace spec;

namespace out
{
   const std::string ENERGY_TOTAL = "EnergyTotal";
#ifdef EN_SUBCAT_OUT
   const std::string ENERGY_INTER = "EnergyInter";
   const std::string ENERGY_TC = "EnergyTC";
   const std::string ENERGY_INTRA_B  = "EnergyIntraBond";
   const std::string ENERGY_INTRA_NB = "EnergyIntraNonbond";
   const std::string ENERGY_ELECT = "EnergyElect";
   const std::string ENERGY_REAL = "EnergyReal";
   const std::string ENERGY_RECIP = "EnergyRecip";
#endif
   const std::string VIRIAL_TOTAL = "VirialTotal";
#ifdef VIR_SUBCAT_OUT
   const std::string VIRIAL_INTER = "VirialInter";
   const std::string VIRIAL_TC = "VirialTC";
#endif
   const std::string PRESSURE = "Pressure";
#ifdef VARIABLE_VOLUME
   const std::string VOLUME = "Volume";
#endif
#if ENSEMBLE == GEMC
   const std::string HEAT_OF_VAP = "HeatOfVaporization";
#endif
#ifdef VARIABLE_DENSITY
   const std::string DENSITY = "Density";
#endif
#ifdef VARIABLE_PARTICLE_NUMBER
   const std::string MOL_NUM = "MolNum";
   const std::string MOL_FRACTION = "MolFraction";
#endif
}
