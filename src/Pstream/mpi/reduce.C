#include "clockTime.H"
#include "vector.H"
#include "vector2D.H"
#include "PstreamReduceOps.H"
#include "PstreamGlobals.H"
#include "UPstream.H"

void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    // syncClockTime clock;
    scalar ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_SCALAR, MPI_SUM, PstreamGlobals::MPICommunicators_[communicator]);
    Value = ret;
    // Info << "specialisations for reduce sumOp<scalar> : " << clock.elapsedTime() << endl;
}


void Foam::reduce
(
    scalar& Value,
    const maxOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    // syncClockTime clock;
    scalar ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_SCALAR, MPI_MAX, PstreamGlobals::MPICommunicators_[communicator]);
    Value = ret;
    // Info << "specialisations for reduce maxOp<scalar> : " << clock.elapsedTime() << endl;
}

void Foam::reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label communicator
)
{
    // syncClockTime clock;
    scalar ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_SCALAR, MPI_MIN, PstreamGlobals::MPICommunicators_[communicator]);
    Value = ret;
    // Info << "specialisations for reduce minOp<scalar> : " << clock.elapsedTime() << endl;
}

void Foam::reduce
(
    label& Value,
    const sumOp<label>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    label ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_LABEL, MPI_SUM, PstreamGlobals::MPICommunicators_[comm]);
    Value = ret;
    // Info << "specialisations for reduce sumOp<label> : " << clock.elapsedTime() << endl;
}


void Foam::reduce
(
    label& Value,
    const minOp<label>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    label ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_LABEL, MPI_MIN, PstreamGlobals::MPICommunicators_[comm]);
    Value = ret;
    // Info << "specialisations for reduce minOp<label> : " << clock.elapsedTime() << endl;
}


void Foam::reduce
(
    label& Value,
    const maxOp<label>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    label ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_LABEL, MPI_MAX, PstreamGlobals::MPICommunicators_[comm]);
    Value = ret;
    // Info << "specialisations for reduce maxOp<label> : " << clock.elapsedTime() << endl;
}


void Foam::reduce
(
    bool& Value,
    const orOp<bool>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    bool ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_C_BOOL, MPI_LOR, PstreamGlobals::MPICommunicators_[comm]);
    Value = ret;
    // Info << "specialisations for reduce orOp<bool> : " << clock.timeIncrement() << endl;
}


void Foam::reduce
(
    bool& Value,
    const andOp<bool>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    bool ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_C_BOOL, MPI_LAND, PstreamGlobals::MPICommunicators_[comm]);
    Value = ret;
    // Info << "specialisations for reduce andOp<bool> : " << clock.timeIncrement() << endl;
}


void Foam::reduce
(
    vector& Value,
    const sumOp<vector>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[3];
    scalar global[3];
    local[0] = Value.x();
    local[1] = Value.y();
    local[2] = Value.z();
    
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_SUM, PstreamGlobals::MPICommunicators_[comm]);

    Value.x() = global[0];
    Value.y() = global[1];
    Value.z() = global[2];
    // Info << "specialisations for reduce sumOp<vector> : " << clock.elapsedTime() << endl;
}

void Foam::reduce
(
    vector& Value,
    const maxOp<vector>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[3];
    scalar global[3];
    local[0] = Value.x();
    local[1] = Value.y();
    local[2] = Value.z();
    
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_MAX, PstreamGlobals::MPICommunicators_[comm]);

    Value.x() = global[0];
    Value.y() = global[1];
    Value.z() = global[2];
    // Info << "specialisations for reduce maxOp<vector> : " << clock.timeIncrement() << endl;
}

void Foam::reduce
(
    vector& Value,
    const minOp<vector>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[3];
    scalar global[3];
    local[0] = Value.x();
    local[1] = Value.y();
    local[2] = Value.z();
    
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_MIN, PstreamGlobals::MPICommunicators_[comm]);

    Value.x() = global[0];
    Value.y() = global[1];
    Value.z() = global[2];
    // Info << "specialisations for reduce minOp<vector> : " << clock.timeIncrement() << endl;
}

void Foam::reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[2];
    scalar global[2];
    local[0] = Value.x();
    local[1] = Value.y();
    
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, PstreamGlobals::MPICommunicators_[comm]);

    Value.x() = global[0];
    Value.y() = global[1];
    // Info << "specialisations for reduce sumOp<vector> : " << clock.timeIncrement() << endl;
}

void Foam::reduce
(
    vector2D& Value,
    const maxOp<vector2D>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[2];
    scalar global[2];
    local[0] = Value.x();
    local[1] = Value.y();
    
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, PstreamGlobals::MPICommunicators_[comm]);

    Value.x() = global[0];
    Value.y() = global[1];
    // Info << "specialisations for reduce maxOp<vector> : " << clock.elapsedTime() << endl;
}

void Foam::reduce
(
    vector2D& Value,
    const minOp<vector2D>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[2];
    scalar global[2];
    local[0] = Value.x();
    local[1] = Value.y();
    
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MIN, PstreamGlobals::MPICommunicators_[comm]);

    Value.x() = global[0];
    Value.y() = global[1];
    // Info << "specialisations for reduce minOp<vector> : " << clock.elapsedTime() << endl;
}

Foam::scalar Foam::returnReduce
(
    const scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_DOUBLE, MPI_SUM, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce sumOp<scalar> : " << clock.elapsedTime() << endl;
    return ret;
}


Foam::scalar Foam::returnReduce
(
    const scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_DOUBLE, MPI_MIN, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce minOp<scalar> : " << clock.elapsedTime() << endl;
    return ret;
}


Foam::scalar Foam::returnReduce
(
    const scalar& Value,
    const maxOp<scalar>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_DOUBLE, MPI_MAX, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce maxOp<scalar> : " << clock.elapsedTime() << endl;
    return ret;
}


Foam::label Foam::returnReduce
(
    const label& Value,
    const sumOp<label>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    label ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_LABEL, MPI_SUM, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce sumOp<label> : " << clock.elapsedTime() << endl;
    return ret;
}


Foam::label Foam::returnReduce
(
    const label& Value,
    const minOp<label>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    label ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_LABEL, MPI_MIN, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce minOp<label> : " << clock.elapsedTime() << endl;
    return ret;
}


Foam::label Foam::returnReduce
(
    const label& Value,
    const maxOp<label>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    label ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_LABEL, MPI_MAX, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce maxOp<label> : " << clock.elapsedTime() << endl;
    return ret;
}


bool Foam::returnReduce
(
    const bool& Value,
    const orOp<bool>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    bool ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_C_BOOL, MPI_LOR, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce orOp<bool> : " << clock.elapsedTime() << endl;
    return ret;
}


bool Foam::returnReduce
(
    const bool& Value,
    const andOp<bool>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    bool ret;
    MPI_Allreduce(&Value, &ret, 1, MPI_C_BOOL, MPI_LAND, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce andOp<bool> : " << clock.elapsedTime() << endl;
    return ret;
}


Foam::vector Foam::returnReduce
(
    const vector& Value,
    const sumOp<vector>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[3];
    scalar global[3];
    local[0] = Value.x();
    local[1] = Value.y();
    local[2] = Value.z();
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_SUM, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce sumOp<vector> : " << clock.elapsedTime() << endl;
    return vector(global[0], global[1], global[2]);
}

Foam::vector Foam::returnReduce
(
    const vector& Value,
    const maxOp<vector>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[3];
    scalar global[3];
    local[0] = Value.x();
    local[1] = Value.y();
    local[2] = Value.z();
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_MAX, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce maxOp<vector> : " << clock.elapsedTime() << endl;
    return vector(global[0], global[1], global[2]);
}

Foam::vector Foam::returnReduce
(
    const vector& Value,
    const minOp<vector>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[3];
    scalar global[3];
    local[0] = Value.x();
    local[1] = Value.y();
    local[2] = Value.z();
    MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_MIN, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce minOp<vector> : " << clock.elapsedTime() << endl;
    return vector(global[0], global[1], global[2]);
}

Foam::vector2D Foam::returnReduce
(
    const vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[2];
    scalar global[2];
    local[0] = Value.x();
    local[1] = Value.y();
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce sumOp<vector> : " << clock.elapsedTime() << endl;
    return vector2D(global[0], global[1]);
}

Foam::vector2D Foam::returnReduce
(
    const vector2D& Value,
    const maxOp<vector2D>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[2];
    scalar global[2];
    local[0] = Value.x();
    local[1] = Value.y();
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MAX, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce maxOp<vector> : " << clock.elapsedTime() << endl;
    return vector2D(global[0], global[1]);
}

Foam::vector2D Foam::returnReduce
(
    const vector2D& Value,
    const minOp<vector2D>& bop,
    const int tag,
    const label comm
){
    // syncClockTime clock;
    scalar local[2];
    scalar global[2];
    local[0] = Value.x();
    local[1] = Value.y();
    MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_MIN, PstreamGlobals::MPICommunicators_[comm]);
    // Info << "specialisations for returnReduce minOp<vector> : " << clock.elapsedTime() << endl;
    return vector2D(global[0], global[1]);
}

void Foam::sumReduce
(
    scalar& Value,
    label& Count,
    const int tag,
    const label communicator
)
{
    if (UPstream::warnComm != -1 && communicator != UPstream::warnComm)
    {
        Pout<< "** reducing:" << Value << " with comm:" << communicator
            << " warnComm:" << UPstream::warnComm
            << endl;
        error::printStack(Pout);
    }
    vector2D twoScalars(Value, scalar(Count));
    reduce(twoScalars, sumOp<vector2D>(), tag, communicator);

    Value = twoScalars.x();
    Count = twoScalars.y();
}


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator,
    label& requestID
)
{
#ifdef MPIX_COMM_TYPE_SHARED
    // Assume mpich2 with non-blocking collectives extensions. Once mpi3
    // is available this will change.
    MPI_Request request;
    scalar v = Value;
    MPIX_Ireduce
    (
        &v,
        &Value,
        1,
        MPI_SCALAR,
        MPI_SUM,
        0,              // root
        PstreamGlobals::MPICommunicators_[communicator],
        &request
    );

    requestID = PstreamGlobals::outstandingRequests_.size();
    PstreamGlobals::outstandingRequests_.append(request);

    if (UPstream::debug)
    {
        Pout<< "UPstream::allocateRequest for non-blocking reduce"
            << " : request:" << requestID
            << endl;
    }
#else
    // Non-blocking not yet implemented in mpi
    reduce(Value, bop, tag, communicator);
    requestID = -1;
#endif
}