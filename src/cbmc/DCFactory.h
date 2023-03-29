/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.75
Copyright (C) 2022 GOMC Group
A copy of the MIT License can be found in License.txt
along with this program, also can be found at
<https://opensource.org/licenses/MIT>.
********************************************************************************/
#ifndef DCFACTORY_H
#define DCFACTORY_H

#include "BasicTypes.h"

namespace cbmc {
class DCComponent;

class DCFactory {
public:
  virtual DCComponent *MakeComponent(uint previous) = 0;
  virtual ~DCFactory(){};
};
} // namespace cbmc

#endif
