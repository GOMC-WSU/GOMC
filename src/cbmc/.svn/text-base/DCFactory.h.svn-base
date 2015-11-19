#ifndef DCFACTORY_H
#define DCFACTORY_H

#include "../../lib/BasicTypes.h"

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
