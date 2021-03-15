#ifndef ReplicaCommunicator_H
#define ReplicaCommunicator_H

#include "GOMC_Config.h"    //For version number

#if GOMC_LIB_MPI
#include <mpi.h>
#endif
#include "XYZArray.h"

class ReplicaCommunicator{
  public:
    ReplicaCommunicator();
    void exchangePositionsNonBlocking(XYZArray * myPos, int exchangePartner);
    void exchangeCOMsNonBlocking(XYZArray * myCOMs, int exchangePartner);
};

#endif