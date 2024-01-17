#include "Pstream.H"
#include "clockTime.H"
#include <typeinfo>

namespace Foam{


void Pstream::gatherList(List<label>& Values, const int tag, const label comm)
{
    if(UPstream::parRun()){
        // Info << "Pstream::gatherList template label start" << endl;
        syncClockTime clock;
        int mpirank = UPstream::myProcNo(comm);
        int masterNo = UPstream::masterNo();
        MPI_Gather(&Values[mpirank], 1, MPI_LABEL, Values.data(), 1, MPI_LABEL, masterNo, getPstreamCommunicator(comm));
        // gatherListTree(Values, tag, comm);
        Info << "Pstream::gatherList template label time : " << clock.elapsedTime() << endl;
    }
}

void Pstream::scatterList(List<label>& Values, const int tag, const label comm)
{
    if(UPstream::parRun()){
        // Info << "Pstream::scatterList template label start" << endl;
        syncClockTime clock;
        int mpisize = UPstream::nProcs(comm);
        int masterNo = UPstream::masterNo();
        MPI_Bcast(Values.data(), mpisize, MPI_LABEL, masterNo, getPstreamCommunicator(comm));
        // scatterListTree(Values, tag, comm);
        Info << "Pstream::scatterList template label time : " << clock.elapsedTime() << endl;
    }
}

void Pstream::gatherList(List<labelList>& Values, const int tag, const label comm)
{
    if(UPstream::parRun()){
        // Info << "Pstream::gatherList template labelList start" << endl;
        syncClockTime clock;

        // int mpirank = UPstream::myProcNo(comm);
        // int mpisize = UPstream::nProcs(comm);
        // int masterNo = UPstream::masterNo();

        // int sr_count = Values[mpirank].size();

        // label* buffer; 
        
        // if(UPstream::master(comm)){
        //     buffer = new label[sr_count * mpisize];
        // }

        // MPI_Gather(Values[mpirank].data(), sr_count, MPI_LABEL, buffer, sr_count, MPI_LABEL, masterNo, getPstreamCommunicator(comm));

        // if(UPstream::master(comm)){
        //     for(int p = 0; p < mpisize; ++p){
        //         Values[p].resize(sr_count);
        //         std::copy(buffer + p * sr_count, buffer + (p + 1) * sr_count, Values[p].data());
        //     }
        //     delete [] buffer;
        // }


        gatherListTree(Values, tag, comm);

        Info << "Pstream::gatherList template labelList time : " << clock.elapsedTime() << endl;
    }
}

void Pstream::scatterList(List<labelList>& Values, const int tag, const label comm)
{
    if(UPstream::parRun()){
        // Info << "Pstream::scatterList template labelList start" << endl;
        syncClockTime clock;
        // int mpirank = UPstream::myProcNo(comm);
        // int mpisize = UPstream::nProcs(comm);
        // int masterNo = UPstream::masterNo();
        // int sr_count = Values[mpirank].size();

        // label* buffer = new label[sr_count * mpisize];

        // if(UPstream::master(comm)){
        //     for(int p = 0; p < mpisize; ++p){
        //         Values[p].resize(sr_count);
        //         std::copy(Values[p].begin(), Values[p].end(), buffer + p * sr_count);
        //     }
        // }

        // MPI_Bcast(buffer, sr_count * mpisize, MPI_LABEL, masterNo, getPstreamCommunicator(comm));

        // if(!UPstream::master(comm)){
        //     for(int p = 0; p < mpisize; ++p){
        //         Values[p].resize(sr_count);
        //         std::copy(buffer + p * sr_count, buffer + (p + 1) * sr_count, Values[p].data());
        //     }
        // }

        // delete [] buffer;
        scatterListTree(Values, tag, comm);
        Info << "Pstream::scatterList template labelList time : " << clock.elapsedTime() << endl;
    }
}

}
