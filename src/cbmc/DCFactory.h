/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.40
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCFACTORY_H
#define DCFACTORY_H

#include "BasicTypes.h"

namespace cbmc
{
class DCComponent;

class DCFactory
{
public:
  virtual DCComponent* MakeComponent(uint previous) = 0;
  virtual ~DCFactory() {};
};
}

#endif
