#ifndef ReplicaCommunicator_H
#define ReplicaCommunicator_H

#include "GOMC_Config.h"    //For version number

#if GOMC_LIB_MPI
#include <mpi.h>
#endif
#include "XYZArray.h"
#include "Coordinates.h"
#include "COM.h"
#include <cstdint>
#include <cstring>

class ReplicaCommunicator{
  public:
  #if GOMC_LIB_MPI
    ReplicaCommunicator();
    void exchangeXYZArrayNonBlocking(XYZArray * myXYZArray, int exchangePartner);
    void exchangeBoxDimensionsNonBlocking(BoxDimensions * myBoxDimensions, int exchangePartner);

  #endif
};

#endif