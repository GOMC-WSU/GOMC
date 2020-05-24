/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "ParallelTemperingUtilities.h"

#if GOMC_LIB_MPI

ParallelTemperingUtilities::ParallelTemperingUtilities(MultiSim const*const& multisim, System & sys):
ms(multisim), sysPotRef(sys.potential){

    #if BOX_TOTAL == 1
        global_energies.resize(multisim->worldSize, 0.0);
    #else
        global_energies.resize(2, std::vector<double>(multisim->worldSize, 0.0));
    #endif

}

bool ParallelTemperingUtilities::evaluateExchangeCriteria(uint step){

    for (int i = 0; i < 2; i++){
        std::cout << "Before fill : energy[" << i << "] : " << global_energies[i] << std::endl;
    }
    #if BOX_TOTAL == 1
    std::memset(&global_energies[0], 0, ms->worldSize * sizeof(double));
    #else
    std::memset(&global_energies[0], 0, ms->worldSize * sizeof(double));
    std::memset(&global_energies[1], 0, ms->worldSize * sizeof(double));
    #endif

    for (int i = 0; i < 2; i++){
        std::cout << "After fill : energy[" << i << "] : " << global_energies[i] << std::endl;
    }

    #if BOX_TOTAL == 1
        global_energies[ms->worldRank] = sysPotRef.boxEnergy[0].total;

    for (int i = 0; i < 2; i++){
        std::cout << "After set local energy : energy[" << i << "] : " << global_energies[i] << std::endl;
    }

    MPI_Allreduce(MPI_IN_PLACE, &global_energies[0], ms->worldSize, MPI_DOUBLE, MPI_SUM,
            MPI_COMM_WORLD);

    //energies = global_energies;

    for (int i = 0; i < 2; i++){
        std::cout << "After allreduce : energy[" << i << "] : " << global_energies[i] << std::endl;
    }

    #else

        global_energies[0][ms->worldRank] = sysPotRef.boxEnergy[0].total;
        global_energies[1][ms->worldRank] = sysPotRef.boxEnergy[1].total;

        MPI_Allreduce(MPI_IN_PLACE, &global_energies[0], ms->worldSize*2, MPI_DOUBLE, MPI_SUM,
              MPI_COMM_WORLD);
        //MPI_Allreduce(MPI_IN_PLACE, &energies[1], ms->worldSize, MPI_DOUBLE, MPI_SUM,
          //    MPI_COMM_WORLD);

    #endif

    //exp(-BETA * (sysPotNew.Total() - sysPotRef.Total()));


}

void ParallelTemperingUtilities::exchangePositions(XYZArray & myPos, MultiSim const*const& multisim){
    
    XYZArray buffer(myPos);

    if (multisim->worldRank == 0){
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myPos.x, myPos.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myPos.x, myPos.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }
    if (multisim->worldRank == 0){
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myPos.y, myPos.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myPos.y, myPos.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }
    if (multisim->worldRank == 0){
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myPos.z, myPos.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myPos.z, myPos.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }

    for (int i = 0; i < 3; i++){
        std::cout << "My coords : " << myPos.x[i] << " x " << myPos.y[i] << " x " << myPos.z[i] << std::endl;
        std::cout << "Send buffer : " << buffer.x[i] << " x " << buffer.y[i] << " x " << buffer.z[i] << std::endl;
    }
}

void ParallelTemperingUtilities::exchangeCOMs(XYZArray & myCOMs, MultiSim const*const& multisim){
    
    XYZArray buffer(myCOMs);

    if (multisim->worldRank == 0){
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myCOMs.x, myCOMs.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myCOMs.x, myCOMs.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }
    if (multisim->worldRank == 0){
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myCOMs.y, myCOMs.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myCOMs.y, myCOMs.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }
    if (multisim->worldRank == 0){
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myCOMs.z, myCOMs.Count(), MPI_DOUBLE, 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myCOMs.z, myCOMs.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }

    for (int i = 0; i < 3; i++){
        std::cout << "My COMs : " << myCOMs.x[i] << " x " << myCOMs.y[i] << " x " << myCOMs.z[i] << std::endl;
        std::cout << "Send buffer : " << buffer.x[i] << " x " << buffer.y[i] << " x " << buffer.z[i] << std::endl;
    }
}


#endif