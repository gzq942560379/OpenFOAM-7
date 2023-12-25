#include "Pstream.H"

namespace Foam{


void Pstream::scatter(label& Value, const int tag, const label comm)
{
    // error::printStack(Pout());
    // Pout << "Pstream::scatter label start" << endl; 
    // if(UPstream::parRun() && UPstream::nProcs(comm) > 1){
    //     // syncClockTime clock;
    //     MPI_Bcast(&Value, 1, MPI_LABEL, 0, getPstreamCommunicator(comm));
    //     // Info << "Pstream::scatter label time : " << clock.elapsedTime() << endl;
    // }
    // Pout << "Pstream::scatter label end" << endl; 
    // Pout << "Pstream::scatter label start" << endl; 
    if (UPstream::parRun() && UPstream::nProcs(comm) > 1){
        clockTime clock;
        if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
        {
            scatter(UPstream::linearCommunication(comm), Value, tag, comm);
        }
        else
        {
            scatter(UPstream::treeCommunication(comm), Value, tag, comm);
        }
        // Info << "Pstream::scatter label time : " << clock.elapsedTime() << endl;
    }
    // Pout << "Pstream::scatter label end" << endl; 
}

void Pstream::scatter(bool& Value, const int tag, const label comm)
{
    // Pout << "Pstream::scatter bool start" << endl; 
    if(UPstream::parRun() && UPstream::nProcs(comm) > 1){
        // syncClockTime clock;
        MPI_Bcast(&Value, 1, MPI_C_BOOL, 0, getPstreamCommunicator(comm));
        // Info << "Pstream::scatter bool time : " << clock.elapsedTime() << endl;
    }
    // Pout << "Pstream::scatter bool end" << endl; 
}

void Pstream::scatter(scalar& Value, const int tag, const label comm)
{
    // Pout << "Pstream::scatter scalar start" << endl; 
    if(UPstream::parRun() && UPstream::nProcs(comm) > 1){
        // syncClockTime clock;
        MPI_Bcast(&Value, 1, MPI_SCALAR, 0, getPstreamCommunicator(comm));
        // Info << "Pstream::scatter scalar time : " << clock.elapsedTime() << endl;
    }
    // Pout << "Pstream::scatter scalar end" << endl; 
}

void Pstream::scatter(string& Value, const int tag, const label comm)
{
    // Pout << "Pstream::scatter string start" << endl; 
    // if(UPstream::parRun() && UPstream::nProcs(comm) > 1){
    //     syncClockTime clock;
    //     int size;
    //     if(UPstream::master(comm)){
    //         size = Value.size();
    //     }
    //     MPI_Bcast(&size, 1, MPI_INT, 0, getPstreamCommunicator(comm));
    //     char buffer[size];
    //     if(UPstream::master(comm)){
    //         std::copy(Value.begin(), Value.end(), buffer);
    //     }
    //     MPI_Bcast(buffer, size, MPI_CHAR, 0, getPstreamCommunicator(comm));
    //     Info << "Pstream::scatter string time : " << clock.elapsedTime() << endl;
    // }

    if (UPstream::parRun() && UPstream::nProcs(comm) > 1){
        clockTime clock;
        if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
        {
            scatter(UPstream::linearCommunication(comm), Value, tag, comm);
        }
        else
        {
            scatter(UPstream::treeCommunication(comm), Value, tag, comm);
        }
        // Info << "Pstream::scatter string time : " << clock.elapsedTime() << endl;
    }
    // Pout << "Pstream::scatter string end" << endl; 
}

void Pstream::scatter(word& Value, const int tag, const label comm)
{
    // Pout << "Pstream::scatter word start" << endl; 
    // if(UPstream::parRun() && UPstream::nProcs(comm) > 1){
    //     syncClockTime clock;
    //     int size;
    //     if(UPstream::master(comm)){
    //         size = Value.size();
    //     }
    //     MPI_Bcast(&size, 1, MPI_INT, 0, getPstreamCommunicator(comm));
    //     char buffer[size];
    //     if(UPstream::master(comm)){
    //         std::copy(Value.begin(), Value.end(), buffer);
    //     }
    //     MPI_Bcast(buffer, size, MPI_CHAR, 0, getPstreamCommunicator(comm));
    //     Info << "Pstream::scatter word time : " << clock.elapsedTime() << endl;
    // }
    if (UPstream::parRun() && UPstream::nProcs(comm) > 1){
        clockTime clock;
        if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
        {
            scatter(UPstream::linearCommunication(comm), Value, tag, comm);
        }
        else
        {
            scatter(UPstream::treeCommunication(comm), Value, tag, comm);
        }
        // Info << "Pstream::scatter word time : " << clock.elapsedTime() << endl;
    }
    // Pout << "Pstream::scatter word end" << endl; 
}

void Pstream::scatter(fileName& Value, const int tag, const label comm)
{
    // Pout << "Pstream::scatter fileName start" << endl; 
    // if(UPstream::parRun() && UPstream::nProcs(comm) > 1){
    //     syncClockTime clock;
    //     int size;
    //     if(UPstream::master(comm)){
    //         size = Value.size();
    //     }
    //     MPI_Bcast(&size, 1, MPI_INT, 0, getPstreamCommunicator(comm));
    //     char buffer[size];
    //     if(UPstream::master(comm)){
    //         std::copy(Value.begin(), Value.end(), buffer);
    //     }
    //     MPI_Bcast(buffer, size, MPI_CHAR, 0, getPstreamCommunicator(comm));
    //     Info << "Pstream::scatter fileName time : " << clock.elapsedTime() << endl;
    // }
    if (UPstream::parRun() && UPstream::nProcs(comm) > 1){
        clockTime clock;
        if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
        {
            scatter(UPstream::linearCommunication(comm), Value, tag, comm);
        }
        else
        {
            scatter(UPstream::treeCommunication(comm), Value, tag, comm);
        }
        // Info << "Pstream::scatter fileName time : " << clock.elapsedTime() << endl;
    }
    // Pout << "Pstream::scatter fileName end" << endl; 
}

}

