#include "clockTime.H"
#include "vector.H"
#include "vector2D.H"
#include "PstreamReduceOps.H"
#include "UPstream.H"


void Foam::reduce
(
    scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label communicator
)
{}


void Foam::reduce
(
    scalar& Value,
    const maxOp<scalar>& bop,
    const int tag,
    const label communicator
)
{}

void Foam::reduce
(
    scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label communicator
)
{}

void Foam::reduce
(
    label& Value,
    const sumOp<label>& bop,
    const int tag,
    const label comm
){}


void Foam::reduce
(
    label& Value,
    const minOp<label>& bop,
    const int tag,
    const label comm
){}


void Foam::reduce
(
    label& Value,
    const maxOp<label>& bop,
    const int tag,
    const label comm
){}


void Foam::reduce
(
    bool& Value,
    const orOp<bool>& bop,
    const int tag,
    const label comm
){}


void Foam::reduce
(
    bool& Value,
    const andOp<bool>& bop,
    const int tag,
    const label comm
){}


void Foam::reduce
(
    vector& Value,
    const sumOp<vector>& bop,
    const int tag,
    const label comm
){}

void Foam::reduce
(
    vector& Value,
    const maxOp<vector>& bop,
    const int tag,
    const label comm
){}

void Foam::reduce
(
    vector& Value,
    const minOp<vector>& bop,
    const int tag,
    const label comm
){}

void Foam::reduce
(
    vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label comm
){}

void Foam::reduce
(
    vector2D& Value,
    const maxOp<vector2D>& bop,
    const int tag,
    const label comm
){}

void Foam::reduce
(
    vector2D& Value,
    const minOp<vector2D>& bop,
    const int tag,
    const label comm
){}

Foam::scalar Foam::returnReduce
(
    const scalar& Value,
    const sumOp<scalar>& bop,
    const int tag,
    const label comm
){
    return Value;
}


Foam::scalar Foam::returnReduce
(
    const scalar& Value,
    const minOp<scalar>& bop,
    const int tag,
    const label comm
){
    return Value;
}


Foam::scalar Foam::returnReduce
(
    const scalar& Value,
    const maxOp<scalar>& bop,
    const int tag,
    const label comm
){
    return Value;
}


Foam::label Foam::returnReduce
(
    const label& Value,
    const sumOp<label>& bop,
    const int tag,
    const label comm
){
    return Value;
}


Foam::label Foam::returnReduce
(
    const label& Value,
    const minOp<label>& bop,
    const int tag,
    const label comm
){
    return Value;
}


Foam::label Foam::returnReduce
(
    const label& Value,
    const maxOp<label>& bop,
    const int tag,
    const label comm
){
    return Value;
}


bool Foam::returnReduce
(
    const bool& Value,
    const orOp<bool>& bop,
    const int tag,
    const label comm
){
    return Value;
}


bool Foam::returnReduce
(
    const bool& Value,
    const andOp<bool>& bop,
    const int tag,
    const label comm
){
    return Value;
}


Foam::vector Foam::returnReduce
(
    const vector& Value,
    const sumOp<vector>& bop,
    const int tag,
    const label comm
){
    return Value;
}

Foam::vector Foam::returnReduce
(
    const vector& Value,
    const maxOp<vector>& bop,
    const int tag,
    const label comm
){
    return Value;
}

Foam::vector Foam::returnReduce
(
    const vector& Value,
    const minOp<vector>& bop,
    const int tag,
    const label comm
){
    return Value;
}

Foam::vector2D Foam::returnReduce
(
    const vector2D& Value,
    const sumOp<vector2D>& bop,
    const int tag,
    const label comm
){
    return Value;
}

Foam::vector2D Foam::returnReduce
(
    const vector2D& Value,
    const maxOp<vector2D>& bop,
    const int tag,
    const label comm
){
    return Value;
}

Foam::vector2D Foam::returnReduce
(
    const vector2D& Value,
    const minOp<vector2D>& bop,
    const int tag,
    const label comm
){
    return Value;
}