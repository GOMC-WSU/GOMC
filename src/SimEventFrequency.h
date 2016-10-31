/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.8
Copyright (C) 2016  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef SIM_EVENT_FREQUENCY_H
#define SIM_EVENT_FREQUENCY_H

#include "../lib/BasicTypes.h" //for ulong
#include "ConfigSetup.h" //for event frequencies from config file.

struct SimEventFrequency
{
   ulong total, perAdjust, tillEquil;
   void Init(config_setup::Step const& s)
   { total = s.total; perAdjust = s.adjustment; tillEquil = s.equil; }
};

#endif /*SIM_EVENT_FREQUENCY_H*/
