/******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) Copyright (C) GOMC Group
A copy of the MIT License can be found in License.txt with this program or at
<https://opensource.org/licenses/MIT>.
******************************************************************************/
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

#endif /*DCFACTORY_H*/
