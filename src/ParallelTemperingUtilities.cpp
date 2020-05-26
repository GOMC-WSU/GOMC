/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "ParallelTemperingUtilities.h"

#if GOMC_LIB_MPI

ParallelTemperingUtilities::ParallelTemperingUtilities(MultiSim const*const& multisim, System & sys, StaticVals const& statV, ulong parallelTempFreq):
ms(multisim), sysPotRef(sys.potential), parallelTempFreq(parallelTempFreq), prng(sys.prngParallelTemp){

    #if BOX_TOTAL == 1
        global_energies.resize(ms->worldSize, 0.0);
    #else
        global_energies.resize(2, std::vector<double>(ms->worldSize, 0.0));
    #endif
    global_betas.resize(ms->worldSize, 0.0);
    exchangeProbabilities.resize(ms->worldSize, 0.0);
    exchangeResults.resize(ms->worldSize, false);
    global_betas[ms->worldRank] = statV.forcefield.beta;
    MPI_Allreduce(MPI_IN_PLACE, &global_betas[0], ms->worldSize, MPI_DOUBLE, MPI_SUM,
        MPI_COMM_WORLD);
}

vector<bool> ParallelTemperingUtilities::evaluateExchangeCriteria(ulong step){

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

    int parity = step / parallelTempFreq % 2;
    double uBoltz;
    bool bPrint = false;
    double printRecord;
    for (int i = 1; i < ms->worldSize; i++){
        bPrint = ms->worldRank == i || ms->worldRank == i -1;
        if (i % 2 == parity){
            uBoltz = exp((global_betas[i-1] - global_betas[i]) * (global_energies[i-1] - global_energies[i]));
            exchangeProbabilities[i] = min(uBoltz, 1.0);
            exchangeResults[i] = (printRecord = prng()) < uBoltz;
            std::cout << "Swapping repl " << i-1 << " and repl " << i << " uBoltz :" << uBoltz << "prng : " << printRecord << std::endl;
        } else {
            exchangeResults[i] = false;
            exchangeProbabilities[i] = 0.0;
        }
    }
    return exchangeResults; 

}

void ParallelTemperingUtilities::exchangePositions(XYZArray & myPos, MultiSim const*const& multisim, int exchangePartner, bool leader){
    
    XYZArray buffer(myPos);
// if im 0, im the follower and i get 1 as a
    if (leader){
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myPos.x, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myPos.x, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    if (leader){
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myPos.y, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myPos.y, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    if (leader){
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myPos.z, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myPos.z, myPos.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

    for (int i = 0; i < 3; i++){
        std::cout << "My coords : " << myPos.x[i] << " x " << myPos.y[i] << " x " << myPos.z[i] << std::endl;
        std::cout << "Send buffer : " << buffer.x[i] << " x " << buffer.y[i] << " x " << buffer.z[i] << std::endl;
    }
}

void ParallelTemperingUtilities::exchangeCOMs(XYZArray & myCOMs, MultiSim const*const& multisim, int exchangePartner, bool leader){
    
    XYZArray buffer(myCOMs);

    if (leader){
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myCOMs.x, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myCOMs.x, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.x, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    if (leader){
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myCOMs.y, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myCOMs.y, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.y, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    if (leader){
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(myCOMs.z, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(myCOMs.z, myCOMs.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(buffer.z, buffer.Count(), MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

    for (int i = 0; i < 3; i++){
        std::cout << "My COMs : " << myCOMs.x[i] << " x " << myCOMs.y[i] << " x " << myCOMs.z[i] << std::endl;
        std::cout << "Send buffer : " << buffer.x[i] << " x " << buffer.y[i] << " x " << buffer.z[i] << std::endl;
    }
}

#endif