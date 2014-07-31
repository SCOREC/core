#include "IPComMan.h"
#include <stdio.h>
#include <set>

using std::set;


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if ( 0 == rank ) {
        printf("|V| = %d\n", commSize);
    }

    set<int> adjRanks;
    adjRanks.insert(rank);
    adjRanks.insert((rank+1)%commSize);

    IPComMan *ipcm = new IPComMan(MPI_COMM_WORLD, adjRanks, 19880, 64, IPComMan::Neighbors);

    const int iMsgSize = sizeof (int);
    int* piMsgSend = (int*) ipcm->alloc_msg(iMsgSize);
    ipcm->set_fixed_msg_size(iMsgSize);
    *piMsgSend = rank;
    ipcm->send((rank+1)%commSize, (void*) piMsgSend);

    *piMsgSend = -1;
    ipcm->send(rank, (void*) piMsgSend);

    ipcm->finalize_send();
    ipcm->free_msg(piMsgSend);

    void *pvMsgRecv;
    int iPidFrom;
    while (int rc = ipcm->receive(pvMsgRecv, &iPidFrom)) {
        int iMsg = *((int*) pvMsgRecv);
        printf("%d received %d from rank %d\n", rank, iMsg, iPidFrom);
        ipcm->free_msg(pvMsgRecv);
    }

    delete ipcm;
    MPI_Finalize();
    return 0;
}
