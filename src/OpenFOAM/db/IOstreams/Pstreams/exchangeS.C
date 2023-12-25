#include "Pstream.H"
#include "clockTime.H"
#include "labelList.H"
#include "DynamicList.H"
#include <cassert>
#include <mpi.h>

namespace Foam{

// template<>
// void Pstream::exchange<labelList, label>
// (
//     const UList<labelList>& sendBufs,
//     List<labelList>& recvBufs,
//     const int tag,
//     const label comm,
//     const bool block
// )
// {
//     syncClockTime clock;
//     assert(sendBufs.size() == UPstream::nProcs(comm));

//     recvBufs.setSize(sendBufs.size());

//     int mpisize = UPstream::nProcs(comm);
//     int mpirank = UPstream::myProcNo(comm);

//     int* sendcounts = new int[mpisize];
//     int* sdispls = new int[mpisize];
//     int* recvcounts = new int[mpisize];
//     int* rdispls = new int[mpisize];
//     int total_send_count = 0;

//     for(int p = 0; p < mpisize; ++p){
//         sendcounts[p] = sendBufs[p].size();
//         total_send_count += sendcounts[p];
//     }
    
//     MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, getPstreamCommunicator(comm));

//     label* sendbuf = new label[total_send_count];

//     sdispls[0] = 0;
//     for(int p = 1; p < mpisize; ++p){
//         sdispls[p] = sdispls[p - 1] + sendcounts[p - 1];
//     }

//     for(int p = 0; p < mpisize; ++p){
//         std::copy(sendBufs[p].begin(), sendBufs[p].end(), sendbuf + sdispls[p]);
//     }

//     int total_recv_count = 0;
//     for(int p = 0; p < mpisize; ++p){
//         total_recv_count += recvcounts[p];
//     }

//     label* recvbuf = new label[total_recv_count];

//     rdispls[0] = 0;
//     for(int p = 1; p < mpisize; ++p){
//         rdispls[p] = rdispls[p - 1] + recvcounts[p - 1];
//     }

//     MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_LABEL, recvbuf, recvcounts, rdispls, MPI_LABEL, getPstreamCommunicator(comm));

//     for(int p = 0; p < mpisize; ++p){
//         recvBufs[p].resize(recvcounts[p]);
//         std::copy(recvbuf + rdispls[p], recvbuf + rdispls[p] + recvcounts[p], recvBufs[p].data());
//     }

//     delete [] sendbuf;
//     delete [] recvbuf;
//     delete [] sendcounts;
//     delete [] recvcounts;
//     delete [] sdispls;
//     delete [] rdispls;
//     Info << "Foam::Pstream::exchange template labelList time : " << clock.elapsedTime() << endl;
// }

// template<>
// void Pstream::exchange<DynamicList<char>, char>
// (
//     const UList<DynamicList<char>>& sendBufs,
//     List<DynamicList<char>>& recvBufs,
//     const int tag,
//     const label comm,
//     const bool block
// )
// {
//     syncClockTime clock;
//     assert(sendBufs.size() == UPstream::nProcs(comm));

//     recvBufs.setSize(sendBufs.size());

//     int mpisize = UPstream::nProcs(comm);
//     int mpirank = UPstream::myProcNo(comm);

//     int* sendcounts = new int[mpisize];
//     int* recvcounts = new int[mpisize];
//     int* sdispls = new int[mpisize];
//     int* rdispls = new int[mpisize];

//     int total_send_count = 0;

//     int count = 0;
//     for(int p = 0; p < mpisize; ++p){
//         sendcounts[p] = sendBufs[p].size();
//         total_send_count += sendcounts[p];
//         if(sendcounts[p] > 0){
//             count += 1;
//         }
//     }
    
//     char* sendbuf = new char[total_send_count];

//     sdispls[0] = 0;
//     for(int p = 1; p < mpisize; ++p){
//         sdispls[p] = sdispls[p - 1] + sendcounts[p - 1];
//     }

//     for(int p = 0; p < mpisize; ++p){
//         std::copy(sendBufs[p].begin(), sendBufs[p].end(), sendbuf + sdispls[p]);
//     }
    
//     MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, getPstreamCommunicator(comm));

//     int total_recv_count = 0;
//     for(int p = 0; p < mpisize; ++p){
//         total_recv_count += recvcounts[p];
//     }

//     char* recvbuf = new char[total_recv_count];

//     rdispls[0] = 0;
//     for(int p = 1; p < mpisize; ++p){
//         rdispls[p] = rdispls[p - 1] + recvcounts[p - 1];
//     }

//     MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_CHAR, recvbuf, recvcounts, rdispls, MPI_CHAR, getPstreamCommunicator(comm));

//     for(int p = 0; p < mpisize; ++p){
//         recvBufs[p].resize(recvcounts[p]);
//         std::copy(recvbuf + rdispls[p], recvbuf + rdispls[p] + recvcounts[p], recvBufs[p].data());
//     }

//     delete [] sendbuf;
//     delete [] recvbuf;
//     delete [] sendcounts;
//     delete [] recvcounts;
//     delete [] sdispls;
//     delete [] rdispls;
//     Info << "Foam::Pstream::exchange template DynamicList<char> time : " << clock.elapsedTime() << endl;
// }


void Pstream::exchange
(
    const UList<labelList>& sendBufs,
    List<labelList>& recvBufs,
    const int tag,
    const label comm,
    const bool block
)
{
    syncClockTime clock;
    assert(sendBufs.size() == UPstream::nProcs(comm));

    recvBufs.setSize(sendBufs.size());

    int mpisize = UPstream::nProcs(comm);
    int mpirank = UPstream::myProcNo(comm);

    int* sendcounts = new int[mpisize];
    int* sdispls = new int[mpisize];
    int* recvcounts = new int[mpisize];
    int* rdispls = new int[mpisize];
    int total_send_count = 0;

    for(int p = 0; p < mpisize; ++p){
        sendcounts[p] = sendBufs[p].size();
        total_send_count += sendcounts[p];
    }
    
    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, getPstreamCommunicator(comm));

    label* sendbuf = new label[total_send_count];

    sdispls[0] = 0;
    for(int p = 1; p < mpisize; ++p){
        sdispls[p] = sdispls[p - 1] + sendcounts[p - 1];
    }

    for(int p = 0; p < mpisize; ++p){
        std::copy(sendBufs[p].begin(), sendBufs[p].end(), sendbuf + sdispls[p]);
    }

    int total_recv_count = 0;
    for(int p = 0; p < mpisize; ++p){
        total_recv_count += recvcounts[p];
    }

    label* recvbuf = new label[total_recv_count];

    rdispls[0] = 0;
    for(int p = 1; p < mpisize; ++p){
        rdispls[p] = rdispls[p - 1] + recvcounts[p - 1];
    }

    MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_LABEL, recvbuf, recvcounts, rdispls, MPI_LABEL, getPstreamCommunicator(comm));

    for(int p = 0; p < mpisize; ++p){
        recvBufs[p].resize(recvcounts[p]);
        std::copy(recvbuf + rdispls[p], recvbuf + rdispls[p] + recvcounts[p], recvBufs[p].data());
    }

    delete [] sendbuf;
    delete [] recvbuf;
    delete [] sendcounts;
    delete [] recvcounts;
    delete [] sdispls;
    delete [] rdispls;
    Info << "Foam::Pstream::exchange labelList time : " << clock.elapsedTime() << endl;
}

void Pstream::exchange
(
    const UList<DynamicList<char>>& sendBufs,
    List<DynamicList<char>>& recvBufs,
    const int tag,
    const label comm,
    const bool block
)
{
    syncClockTime clock;
    assert(sendBufs.size() == UPstream::nProcs(comm));

    recvBufs.setSize(sendBufs.size());

    int mpisize = UPstream::nProcs(comm);
    int mpirank = UPstream::myProcNo(comm);

    int* sendcounts = new int[mpisize];
    int* recvcounts = new int[mpisize];
    int* sdispls = new int[mpisize];
    int* rdispls = new int[mpisize];

    int total_send_count = 0;

    int count = 0;
    for(int p = 0; p < mpisize; ++p){
        sendcounts[p] = sendBufs[p].size();
        total_send_count += sendcounts[p];
        if(sendcounts[p] > 0){
            count += 1;
        }
    }
    
    char* sendbuf = new char[total_send_count];

    sdispls[0] = 0;
    for(int p = 1; p < mpisize; ++p){
        sdispls[p] = sdispls[p - 1] + sendcounts[p - 1];
    }

    for(int p = 0; p < mpisize; ++p){
        std::copy(sendBufs[p].begin(), sendBufs[p].end(), sendbuf + sdispls[p]);
    }
    
    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, getPstreamCommunicator(comm));

    int total_recv_count = 0;
    for(int p = 0; p < mpisize; ++p){
        total_recv_count += recvcounts[p];
    }

    char* recvbuf = new char[total_recv_count];

    rdispls[0] = 0;
    for(int p = 1; p < mpisize; ++p){
        rdispls[p] = rdispls[p - 1] + recvcounts[p - 1];
    }

    MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_CHAR, recvbuf, recvcounts, rdispls, MPI_CHAR, getPstreamCommunicator(comm));

    for(int p = 0; p < mpisize; ++p){
        recvBufs[p].resize(recvcounts[p]);
        std::copy(recvbuf + rdispls[p], recvbuf + rdispls[p] + recvcounts[p], recvBufs[p].data());
    }

    delete [] sendbuf;
    delete [] recvbuf;
    delete [] sendcounts;
    delete [] recvcounts;
    delete [] sdispls;
    delete [] rdispls;
    Info << "Foam::Pstream::exchange DynamicList<char> time : " << clock.elapsedTime() << endl;
}

}
