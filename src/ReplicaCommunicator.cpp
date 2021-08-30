#include "ReplicaCommunicator.h"

#if GOMC_LIB_MPI
ReplicaCommunicator::ReplicaCommunicator(){}

void ReplicaCommunicator::exchangeXYZArrayNonBlocking(XYZArray * myXYZArray, int exchangePartner)
{
    XYZArray outBuffer(*myXYZArray);
    MPI_Request mpi_req;

    int myXYZArrayCount = outBuffer.Count();
    int otherXYZCount;

    MPI_Isend(&myXYZArrayCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(&otherXYZCount, 1 * sizeof(int), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    XYZArray inBuffer(otherXYZCount);

    MPI_Isend(outBuffer.x, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.x, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.y, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.y, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    MPI_Isend(outBuffer.z, outBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req);
    MPI_Recv(inBuffer.z, inBuffer.Count() * sizeof(double), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);

    swap(*myXYZArray, inBuffer);
}

void ReplicaCommunicator::exchangeMolLookupNonBlocking(MoleculeLookup * newMolLookup, MoleculeLookup * oldMolLookup, int exchangePartner)
{
    MPI_Request mpi_req0;
    MPI_Request mpi_req1;
    MPI_Request mpi_req2;
    MPI_Request mpi_req3;
    MPI_Request mpi_req4;
    MPI_Request mpi_req5;
    MPI_Request mpi_req6;
    MPI_Request mpi_req7;
    MPI_Request mpi_req8;
    MPI_Request mpi_req9;
    MPI_Request mpi_req10;
    MPI_Request mpi_req11;
    MPI_Request mpi_req12;
    MPI_Request mpi_req13;

    MPI_Isend(oldMolLookup->atomCharge.data(), oldMolLookup->atomCharge.size() * sizeof(decltype(oldMolLookup->atomCharge.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req0);
    MPI_Isend(oldMolLookup->atomIndex.data(), oldMolLookup->atomIndex.size() * sizeof(decltype(oldMolLookup->atomIndex.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req1);
    MPI_Isend(oldMolLookup->atomKind.data(), oldMolLookup->atomKind.size() * sizeof(decltype(oldMolLookup->atomKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req2);
    MPI_Isend(oldMolLookup->boxAndKindStart.data(), oldMolLookup->boxAndKindStart.size() * sizeof(decltype(oldMolLookup->boxAndKindStart.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req3);
    MPI_Isend(&oldMolLookup->boxAndKindStartCount, 1 * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req4);
    MPI_Isend(oldMolLookup->boxAndKindSwappableCounts.data(), oldMolLookup->boxAndKindSwappableCounts.size() * sizeof(decltype(oldMolLookup->boxAndKindSwappableCounts.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req5);
    MPI_Isend(oldMolLookup->canMoveKind.data(), oldMolLookup->canMoveKind.size() * sizeof(decltype(oldMolLookup->canMoveKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req6);
    MPI_Isend(oldMolLookup->canSwapKind.data(), oldMolLookup->canSwapKind.size() * sizeof(decltype(oldMolLookup->canSwapKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req7);
    MPI_Isend(oldMolLookup->fixedMolecule.data(), oldMolLookup->fixedMolecule.size() * sizeof(decltype(oldMolLookup->fixedMolecule.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req8);
    MPI_Isend(oldMolLookup->molIndex.data(), oldMolLookup->molIndex.size() * sizeof(decltype(oldMolLookup->molIndex.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req9);
    MPI_Isend(oldMolLookup->molKind.data(), oldMolLookup->molKind.size() * sizeof(decltype(oldMolLookup->molKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req10);
    MPI_Isend(oldMolLookup->molLookup.data(), oldMolLookup->molLookup.size() * sizeof(decltype(oldMolLookup->molLookup.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req11);
    MPI_Isend(&oldMolLookup->molLookupCount, 1 * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req12);
    MPI_Isend(&oldMolLookup->numKinds, 1 * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, &mpi_req13);

    MPI_Recv(newMolLookup->atomCharge.data(), newMolLookup->atomCharge.size() * sizeof(decltype(newMolLookup->atomCharge.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->atomIndex.data(), newMolLookup->atomIndex.size() * sizeof(decltype(newMolLookup->atomIndex.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->atomKind.data(), newMolLookup->atomKind.size() * sizeof(decltype(newMolLookup->atomKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->boxAndKindStart.data(), newMolLookup->boxAndKindStart.size() * sizeof(decltype(newMolLookup->boxAndKindStart.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&newMolLookup->boxAndKindStartCount, 1 * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->boxAndKindSwappableCounts.data(), newMolLookup->boxAndKindSwappableCounts.size() * sizeof(decltype(newMolLookup->boxAndKindSwappableCounts.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->canMoveKind.data(), newMolLookup->canMoveKind.size() * sizeof(decltype(newMolLookup->canMoveKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->canSwapKind.data(), newMolLookup->canSwapKind.size() * sizeof(decltype(newMolLookup->canSwapKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->fixedMolecule.data(), newMolLookup->fixedMolecule.size() * sizeof(decltype(newMolLookup->fixedMolecule.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->molIndex.data(), newMolLookup->molIndex.size() * sizeof(decltype(newMolLookup->molIndex.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->molKind.data(), newMolLookup->molKind.size() * sizeof(decltype(newMolLookup->molKind.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(newMolLookup->molLookup.data(), newMolLookup->molLookup.size() * sizeof(decltype(newMolLookup->molLookup.back())), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&newMolLookup->molLookupCount, 1 * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&newMolLookup->numKinds, 1 * sizeof(uint), MPI_BYTE, exchangePartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Wait(&mpi_req0, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req1, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req2, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req3, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req4, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req5, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req6, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req7, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req8, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req9, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req10, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req11, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req12, MPI_STATUS_IGNORE);
    MPI_Wait(&mpi_req13, MPI_STATUS_IGNORE);

}

#endif