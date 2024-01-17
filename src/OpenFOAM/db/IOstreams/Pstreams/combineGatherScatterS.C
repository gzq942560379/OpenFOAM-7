#include "Pstream.H"
#include "clockTime.H"
#include <typeinfo>

namespace Foam{


void Pstream::listCombineGather
(
    List<scalar>& Value,
    const maxEqOp<scalar>& cop,
    const int tag,
    const label comm
){
    // Pout << "Pstream::listCombineGather maxEqOp<scalar> start" << endl;
    if(UPstream::parRun()){
        syncClockTime clock;
        List<scalar> global(Value);
        MPI_Reduce(Value.data(), global.data(), Value.size(), MPI_SCALAR, MPI_MAX, 0, getPstreamCommunicator(comm));
        Info << "Foam::Pstream::listCombineGather template scalar maxEqOp time : " << clock.elapsedTime() << endl;
        Value = global;
    }
    // Pout << "Pstream::listCombineGather maxEqOp<scalar> end" << endl;
}

void Pstream::listCombineGather
(
    List<scalar>& Value,
    const minEqOp<scalar>& cop,
    const int tag,
    const label comm
){
    // Pout << "Pstream::listCombineGather minEqOp<scalar> start" << endl;
    if(UPstream::parRun()){
        syncClockTime clock;
        List<scalar> global(Value);
        MPI_Reduce(Value.data(), global.data(), Value.size(), MPI_SCALAR, MPI_MIN, 0, getPstreamCommunicator(comm));
        Info << "Foam::Pstream::listCombineGather template scalar minEqOp time : " << clock.elapsedTime() << endl;
        Value = global;
    }
    // Pout << "Pstream::listCombineGather minEqOp<scalar> end" << endl;
}

void Pstream::listCombineScatter
(
    List<scalar>& Value,
    const int tag,
    const label comm
){
    // Pout << "Pstream::listCombineScatter scalar start" << endl;
    if(UPstream::parRun()){
        syncClockTime clock;
        MPI_Bcast(Value.data(), Value.size(), MPI_SCALAR, 0, getPstreamCommunicator(comm));
        Info << "Foam::Pstream::listCombineScatter template scalar time : " << clock.elapsedTime() << endl;
    }
    // Pout << "Pstream::listCombineScatter scalar start" << endl;
}

}
