/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 2.51
Copyright (C) 2018  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/

#include "ParallelTemperingUtilities.h"

#if GOMC_LIB_MPI

ParallelTemperingUtilities::ParallelTemperingUtilities(MultiSim const*const& multisim, System & sys, StaticVals const& statV, ulong parallelTempFreq):
ms(multisim), sysPotRef(sys.potential), parallelTempFreq(parallelTempFreq), prng(sys.prngParallelTemp), newMolsPos(sys.boxDimRef, newCOMs, sys.molLookupRef, sys.prng, statV.mol),
  newCOMs(sys.boxDimRef, newMolsPos, sys.molLookupRef, statV.mol){

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

    //for (int i = 0; i < 2; i++){
      //  std::cout << "Before fill : energy[" << i << "] : " << global_energies[i] << std::endl;
    //}
    #if BOX_TOTAL == 1
    std::memset(&global_energies[0], 0, ms->worldSize * sizeof(double));
    #else
    std::memset(&global_energies[0], 0, ms->worldSize * sizeof(double));
    std::memset(&global_energies[1], 0, ms->worldSize * sizeof(double));
    #endif

    //for (int i = 0; i < 2; i++){
      //  std::cout << "After fill : energy[" << i << "] : " << global_energies[i] << std::endl;
    //}

    #if BOX_TOTAL == 1
//        global_energies[ms->worldRank] = sysPotRef.boxEnergy[0].total;
        global_energies[ms->worldRank] = sysPotRef.totalEnergy.total;

    //for (int i = 0; i < 2; i++){
      //  std::cout << "After set local energy : energy[" << i << "] : " << global_energies[i] << std::endl;
    //}

    MPI_Allreduce(MPI_IN_PLACE, &global_energies[0], ms->worldSize, MPI_DOUBLE, MPI_SUM,
            MPI_COMM_WORLD);

    //for (int i = 0; i < 2; i++){
     //   std::cout << "After allreduce : energy[" << i << "] : " << global_energies[i] << std::endl;
   // }

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
            #if ENSEMBLE == NVT
            uBoltz = exp((global_betas[i] - global_betas[i-1]) * (global_energies[i] - global_energies[i-1]));
            #endif
            exchangeProbabilities[i] = min(uBoltz, 1.0);
            exchangeResults[i] = (printRecord = prng()) < uBoltz;
            //std::cout << "Swapping repl " << i-1 << " and repl " << i << " uBoltz :" << uBoltz << "prng : " << printRecord << std::endl;
        } else {
            exchangeResults[i] = false;
            exchangeProbabilities[i] = 0.0;
        }
    }
    return exchangeResults; 

}

void ParallelTemperingUtilities::conductExchanges(Coordinates & coordCurrRef, COM & comCurrRef, MultiSim const*const& ms, vector<bool> & resultsOfExchangeCriteria){


        if (resultsOfExchangeCriteria[ms->worldRank] == true){
/*

When you wake up in the morning, if the 50m sim returns that flat line,

try building a coord and com object, and then call

    swap(coordCurrRef, newMolsPos);
    swap(comCurrRef, newCOMs);

*/
                newMolsPos = coordCurrRef;
                newCOMs = comCurrRef;

                std::cout << "Reassigned new pos and coms" << std::endl;

                exchangePositions(newMolsPos, ms, ms->worldRank+1, false);
                exchangeCOMs(newCOMs, ms, ms->worldRank+1, false);

                std::cout << "Exchanged new pos and coms" << std::endl;


                swap(coordCurrRef, newMolsPos);
                swap(comCurrRef, newCOMs);

                std::cout << "Swapped new pos and coms" << std::endl;

/*

Eventually add this back, but I am isolating the swapping from regrid and repot to see if this is the src of the issue

                myCellList.GridAll(system->boxDimRef, newMolsPos, system->molLookup);

                sysPotNew = system->calcEnergy.SystemTotal();

                sysPotRef = sysPotNew;
*/
      } else if(ms->worldRank+1 != ms->worldSize && resultsOfExchangeCriteria[ms->worldRank+1] == true) {

                newMolsPos = coordCurrRef;
                newCOMs = comCurrRef;

                std::cout << "Reassigned new pos and coms" << std::endl;

                exchangePositions(newMolsPos, ms, ms->worldRank+1, false);
                exchangeCOMs(newCOMs, ms, ms->worldRank+1, false);
                std::cout << "Exchanged new pos and coms" << std::endl;


                swap(coordCurrRef, newMolsPos);
                swap(comCurrRef, newCOMs);
                std::cout << "Swapped new pos and coms" << std::endl;

/*

                myCellList.GridAll(system->boxDimRef, newMolsPos, system->molLookup);

                sysPotNew = system->calcEnergy.SystemTotal();

                sysPotRef = sysPotNew;
*/

      }
}


void ParallelTemperingUtilities::exchangePositions(Coordinates & myPos, MultiSim const*const& multisim, int exchangePartner, bool leader){
    
    Coordinates buffer = myPos;
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
}

void ParallelTemperingUtilities::exchangeCOMs(COM & myCOMs, MultiSim const*const& multisim, int exchangePartner, bool leader){
    
    COM buffer = myCOMs;

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
}

void ParallelTemperingUtilities::exchangeCellLists(CellList & myCellList, MultiSim const*const& multisim, int exchangePartner, bool leader){
    
    CellList buffer(myCellList);
    std::cout << " buffer list size : " << buffer.list.size() << " myCellList list Size : " << myCellList.list.size() << std::endl; 
    std::cout << " buffer head size : " << buffer.head[0].size() << " myCellList head Size : " << myCellList.head[0].size() << std::endl; 

// if im 0, im the follower and i get 1 as a

    if (leader){
        MPI_Send(&buffer.list[0], buffer.list.size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&myCellList.list[0], myCellList.list.size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&myCellList.list[0], myCellList.list.size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.list[0], buffer.list.size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    
    if (leader){
        MPI_Send(&buffer.head[0][0], buffer.head[0].size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&myCellList.head[0][0], myCellList.head[0].size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&myCellList.head[0][0], myCellList.head[0].size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.head[0][0], buffer.head[0].size(), MPI_INT, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
}

void ParallelTemperingUtilities::exchangePotentials(SystemPotential & mySystemPotential, MultiSim const*const& multisim, int exchangePartner, bool leader){
    
    SystemPotential buffer(mySystemPotential);

// if im 0, im the follower and i get 1 as a

    if (leader){
        MPI_Send(&buffer.totalEnergy.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    
    if (leader){
        MPI_Send(&buffer.totalEnergy.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalEnergy.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalEnergy.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalEnergy.intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalEnergy.intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalEnergy.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalEnergy.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalEnergy.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

    if (leader){
        MPI_Send(&buffer.totalEnergy.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalEnergy.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalEnergy.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalEnergy.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    for (int b = 0; b < BOX_TOTAL; b++){
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }
        
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].intraBond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].intraNonbond, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxEnergy[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxEnergy[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

        if (leader){
            MPI_Send(&buffer.boxEnergy[b].tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxEnergy[b].tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxEnergy[b].tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.totalEnergy.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }
    }
}

void ParallelTemperingUtilities::exchangeVirials(SystemPotential & mySystemPotential, MultiSim const*const& multisim, int exchangePartner, bool leader){
    SystemPotential buffer(mySystemPotential);

// if im 0, im the follower and i get 1 as a

    if (leader){
        MPI_Send(&buffer.totalVirial.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.total, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
    
    if (leader){
        MPI_Send(&buffer.totalVirial.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalVirial.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.correction, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalVirial.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.inter, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }
        
    if (leader){
        MPI_Send(&buffer.totalVirial.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.real, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalVirial.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.recip, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

        
    if (leader){
        MPI_Send(&buffer.totalVirial.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.self, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

    if (leader){
        MPI_Send(&buffer.totalVirial.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
        MPI_Recv(&mySystemPotential.totalVirial.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(&mySystemPotential.totalVirial.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(&buffer.totalVirial.tc, 1, MPI_DOUBLE, exchangePartner, 0,
                 MPI_COMM_WORLD);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (leader){
                MPI_Send(&buffer.totalVirial.totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
                MPI_Recv(&mySystemPotential.totalVirial.totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(&mySystemPotential.totalVirial.totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&buffer.totalVirial.totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
            }

            if (leader){
                MPI_Send(&buffer.totalVirial.interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
                MPI_Recv(&mySystemPotential.totalVirial.interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(&mySystemPotential.totalVirial.interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&buffer.totalVirial.interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
            }

            if (leader){
                MPI_Send(&buffer.totalVirial.realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
                MPI_Recv(&mySystemPotential.totalVirial.realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(&mySystemPotential.totalVirial.realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&buffer.totalVirial.realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
            }

            if (leader){
                MPI_Send(&buffer.totalVirial.recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
                MPI_Recv(&mySystemPotential.totalVirial.recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(&mySystemPotential.totalVirial.recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&buffer.totalVirial.recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
            }

            if (leader){
                MPI_Send(&buffer.totalVirial.corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
                MPI_Recv(&mySystemPotential.totalVirial.corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                MPI_Recv(&mySystemPotential.totalVirial.corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Send(&buffer.totalVirial.corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                        MPI_COMM_WORLD);
            }
        }
    }

    for (int b = 0; b < BOX_TOTAL; b++){
        if (leader){
            MPI_Send(&buffer.boxVirial[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].total, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }
        
        if (leader){
            MPI_Send(&buffer.boxVirial[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].totalElect, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxVirial[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].correction, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxVirial[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].inter, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }
            
        if (leader){
            MPI_Send(&buffer.boxVirial[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].real, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxVirial[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].recip, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

            
        if (leader){
            MPI_Send(&buffer.boxVirial[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].self, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }

        if (leader){
            MPI_Send(&buffer.boxVirial[b].tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
            MPI_Recv(&mySystemPotential.boxVirial[b].tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&mySystemPotential.boxVirial[b].tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&buffer.boxVirial[b].tc, 1, MPI_DOUBLE, exchangePartner, 0,
                    MPI_COMM_WORLD);
        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (leader){
                    MPI_Send(&buffer.boxVirial[b].totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                    MPI_Recv(&mySystemPotential.boxVirial[b].totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(&mySystemPotential.boxVirial[b].totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&buffer.boxVirial[b].totalTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                }

                if (leader){
                    MPI_Send(&buffer.boxVirial[b].interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                    MPI_Recv(&mySystemPotential.boxVirial[b].interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(&mySystemPotential.boxVirial[b].interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&buffer.boxVirial[b].interTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                }

                if (leader){
                    MPI_Send(&buffer.boxVirial[b].realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                    MPI_Recv(&mySystemPotential.boxVirial[b].realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(&mySystemPotential.boxVirial[b].realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&buffer.boxVirial[b].realTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                }

                if (leader){
                    MPI_Send(&buffer.boxVirial[b].recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                    MPI_Recv(&mySystemPotential.boxVirial[b].recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(&mySystemPotential.boxVirial[b].recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&buffer.boxVirial[b].recipTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                }

                if (leader){
                    MPI_Send(&buffer.boxVirial[b].corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                    MPI_Recv(&mySystemPotential.boxVirial[b].corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(&mySystemPotential.boxVirial[b].corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&buffer.boxVirial[b].corrTens[i][j], 1, MPI_DOUBLE, exchangePartner, 0,
                            MPI_COMM_WORLD);
                }
            }
        }
    }
}


#endif