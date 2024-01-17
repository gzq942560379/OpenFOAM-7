/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    Values[UPstream::myProcNo(comm)].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use scatter to distribute the data.

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "contiguous.H"
#include <mpi.h>
#include "clockTime.H"
#include <typeinfo>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// template<class T>
// void Pstream::gatherList
// (
//     const List<UPstream::commsStruct>& comms,
//     List<T>& Values,
//     const int tag,
//     const label comm
// )
// {
//     if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
//     {
//         if (Values.size() != UPstream::nProcs(comm))
//         {
//             FatalErrorInFunction
//                 << "Size of list:" << Values.size()
//                 << " does not equal the number of processors:"
//                 << UPstream::nProcs(comm)
//                 << Foam::abort(FatalError);
//         }

//         // Get my communication order
//         const commsStruct& = comms[UPstream::myProcNo(comm)];

//         // Receive from my downstairs neighbours
//         forAll(myComm.below(), belowI)
//         {
//             label belowID = myComm.below()[belowI];
//             const labelList& belowLeaves = comms[belowID].allBelow();

//             if (contiguous<T>())
//             {
//                 List<T> receivedValues(belowLeaves.size() + 1);

//                 UIPstream::read
//                 (
//                     UPstream::commsTypes::scheduled,
//                     belowID,
//                     reinterpret_cast<char*>(receivedValues.begin()),
//                     receivedValues.byteSize(),
//                     tag,
//                     comm
//                 );

//                 Values[belowID] = receivedValues[0];

//                 forAll(belowLeaves, leafI)
//                 {
//                     Values[belowLeaves[leafI]] = receivedValues[leafI + 1];
//                 }
//             }
//             else
//             {
//                 IPstream fromBelow
//                 (
//                     UPstream::commsTypes::scheduled,
//                     belowID,
//                     0,
//                     tag,
//                     comm
//                 );
//                 fromBelow >> Values[belowID];

//                 if (debug & 2)
//                 {
//                     Pout<< " received through "
//                         << belowID << " data from:" << belowID
//                         << " data:" << Values[belowID] << endl;
//                 }

//                 // Receive from all other processors below belowID
//                 forAll(belowLeaves, leafI)
//                 {
//                     label leafID = belowLeaves[leafI];
//                     fromBelow >> Values[leafID];

//                     if (debug & 2)
//                     {
//                         Pout<< " received through "
//                             << belowID << " data from:" << leafID
//                             << " data:" << Values[leafID] << endl;
//                     }
//                 }
//             }
//         }

//         // Send up from Values:
//         // - my own value first
//         // - all belowLeaves next
//         if (myComm.above() != -1)
//         {
//             const labelList& belowLeaves = myComm.allBelow();

//             if (debug & 2)
//             {
//                 Pout<< " sending to " << myComm.above()
//                     << " data from me:" << UPstream::myProcNo(comm)
//                     << " data:" << Values[UPstream::myProcNo(comm)] << endl;
//             }

//             if (contiguous<T>())
//             {
//                 List<T> sendingValues(belowLeaves.size() + 1);
//                 sendingValues[0] = Values[UPstream::myProcNo(comm)];

//                 forAll(belowLeaves, leafI)
//                 {
//                     sendingValues[leafI + 1] = Values[belowLeaves[leafI]];
//                 }

//                 OPstream::write
//                 (
//                     UPstream::commsTypes::scheduled,
//                     myComm.above(),
//                     reinterpret_cast<const char*>(sendingValues.begin()),
//                     sendingValues.byteSize(),
//                     tag,
//                     comm
//                 );
//             }
//             else
//             {
//                 OPstream toAbove
//                 (
//                     UPstream::commsTypes::scheduled,
//                     myComm.above(),
//                     0,
//                     tag,
//                     comm
//                 );
//                 toAbove << Values[UPstream::myProcNo(comm)];

//                 forAll(belowLeaves, leafI)
//                 {
//                     label leafID = belowLeaves[leafI];

//                     if (debug & 2)
//                     {
//                         Pout<< " sending to "
//                             << myComm.above() << " data from:" << leafID
//                             << " data:" << Values[leafID] << endl;
//                     }
//                     toAbove << Values[leafID];
//                 }
//             }
//         }
//     }
// }

template<class T>
void Pstream::gatherListLinear
(
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if(comm != 0){
        FatalErrorInFunction
            << "gatherListLinear only support comm == 0, but comm = " << comm << endl << flush
            << Foam::abort(FatalError);   
    }

    if (UPstream::parRun() && UPstream::nProcs() > 1)
    {
        if (Values.size() != UPstream::nProcs())
        {
            FatalErrorInFunction
                << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs()
                << Foam::abort(FatalError);
        }

        if(UPstream::master()){
            for(int belowID = firstSlave(); belowID <= lastSlave(); ++belowID){
                IPstream fromBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,
                    tag,
                    comm
                );
                fromBelow >> Values[belowID];
            }
        }else{
            OPstream toAbove
            (
                UPstream::commsTypes::scheduled,
                masterNo(),
                0,
                tag,
                comm
            );
            toAbove << Values[UPstream::myProcNo()];
        }
    } 
}

template<class T>
void Pstream::gatherListTree
(
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if(comm != 0){
        FatalErrorInFunction
            << "gatherListLinear only support comm == 0, but comm = " << comm << endl << flush
            << Foam::abort(FatalError);   
    }

    if (UPstream::parRun() && UPstream::nProcs() > 1)
    {
        if (Values.size() != UPstream::nProcs())
        {
            FatalErrorInFunction
                << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs()
                << Foam::abort(FatalError);
        }

        int nLevels = 1;
        while ((1 << nLevels) < static_cast<int>(UPstream::nProcs()))
        {
            nLevels++;
        }

        for (int level = 0; level < nLevels; ++level){
            int mask = (1 << level);
            if((UPstream::myProcNo() & (mask - 1)) != 0){
                break;
            }
            int peer = UPstream::myProcNo() ^ mask;
            if(peer >= UPstream::nProcs()){
                continue;
            }
            int count = mask;
            if(UPstream::myProcNo() & mask){ // sender
                // std::cout << UPstream::myProcNo() <<  " send (" << UPstream::myProcNo() << ", " << UPstream::myProcNo() + count << ") to " << peer << std::endl;
                if (contiguous<T>()){
                    List<T> sendingValues(count);
                    for(int i = 0; i < count; ++i){
                        sendingValues[i] = Values[UPstream::myProcNo() + i];
                    }
                    OPstream::write
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        reinterpret_cast<const char*>(sendingValues.begin()),
                        sendingValues.byteSize(),
                        tag,
                        comm
                    );
                }else{
                    OPstream toPeer
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        0,
                        tag,
                        comm
                    );
                    for(int i = 0; i < count; ++i){
                        toPeer << Values[UPstream::myProcNo() + i];
                    }
                }
            }else{ // reciver
                // std::cout << UPstream::myProcNo() <<  " recv (" << peer << ", " << peer + count << ") from " << peer << std::endl;
                if (contiguous<T>()){
                    List<T> receivedValues(count);
                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        reinterpret_cast<char*>(receivedValues.begin()),
                        receivedValues.byteSize(),
                        tag,
                        comm
                    );
                    for(int i = 0; i < count; ++i){
                        Values[peer + i] = receivedValues[i];
                    }
                }else{
                    IPstream fromPeer
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        0,
                        tag,
                        comm
                    );
                    for(int i = 0; i < count; ++i){
                        fromPeer >> Values[peer + i];
                    }
                }
            }
        }
    } 
}

template<class T>
void Pstream::gatherList(List<T>& Values, const int tag, const label comm)
{
    syncClockTime clock;
    // Info << "Pstream::gatherList start" << endl;
    // if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    // {
    //     gatherList(UPstream::linearCommunication(comm), Values, tag, comm);
    // }
    // else
    // {
    //     gatherList(UPstream::treeCommunication(comm), Values, tag, comm);
    // }
    // gatherListLinear(Values, tag, comm);
    gatherListTree(Values, tag, comm);
    // Info << "Pstream::gatherList end" << endl;
    Info << "Pstream::gatherList T : " << typeid(T).name() << endl;
    Info << "Pstream::gatherList time : " << clock.elapsedTime() << endl;
}


// template<class T>
// void Pstream::scatterList
// (
//     const List<UPstream::commsStruct>& comms,
//     List<T>& Values,
//     const int tag,
//     const label comm
// )
// {
//     if (UPstream::parRun() && UPstream::nProcs(comm) > 1)
//     {
//         if (Values.size() != UPstream::nProcs(comm))
//         {
//             FatalErrorInFunction
//                 << "Size of list:" << Values.size()
//                 << " does not equal the number of processors:"
//                 << UPstream::nProcs(comm)
//                 << Foam::abort(FatalError);
//         }

//         // Get my communication order
//         const commsStruct& myComm = comms[UPstream::myProcNo(comm)];

//         // Receive from up
//         if (myComm.above() != -1)
//         {
//             const labelList& notBelowLeaves = myComm.allNotBelow();

//             if (contiguous<T>())
//             {
//                 List<T> receivedValues(notBelowLeaves.size());

//                 UIPstream::read
//                 (
//                     UPstream::commsTypes::scheduled,
//                     myComm.above(),
//                     reinterpret_cast<char*>(receivedValues.begin()),
//                     receivedValues.byteSize(),
//                     tag,
//                     comm
//                 );

//                 forAll(notBelowLeaves, leafI)
//                 {
//                     Values[notBelowLeaves[leafI]] = receivedValues[leafI];
//                 }
//             }
//             else
//             {
//                 IPstream fromAbove
//                 (
//                     UPstream::commsTypes::scheduled,
//                     myComm.above(),
//                     0,
//                     tag,
//                     comm
//                 );

//                 forAll(notBelowLeaves, leafI)
//                 {
//                     label leafID = notBelowLeaves[leafI];
//                     fromAbove >> Values[leafID];

//                     if (debug)
//                     {
//                         Pout<< " received through "
//                             << myComm.above() << " data for:" << leafID
//                             << " data:" << Values[leafID] << endl;
//                     }
//                 }
//             }
//         }

//         // Send to my downstairs neighbours
//         forAllReverse(myComm.below(), belowI)
//         {
//             label belowID = myComm.below()[belowI];
//             const labelList& notBelowLeaves = comms[belowID].allNotBelow();

//             if (contiguous<T>())
//             {
//                 List<T> sendingValues(notBelowLeaves.size());

//                 forAll(notBelowLeaves, leafI)
//                 {
//                     sendingValues[leafI] = Values[notBelowLeaves[leafI]];
//                 }

//                 OPstream::write
//                 (
//                     UPstream::commsTypes::scheduled,
//                     belowID,
//                     reinterpret_cast<const char*>(sendingValues.begin()),
//                     sendingValues.byteSize(),
//                     tag,
//                     comm
//                 );
//             }
//             else
//             {
//                 OPstream toBelow
//                 (
//                     UPstream::commsTypes::scheduled,
//                     belowID,
//                     0,
//                     tag,
//                     comm
//                 );

//                 // Send data destined for all other processors below belowID
//                 forAll(notBelowLeaves, leafI)
//                 {
//                     label leafID = notBelowLeaves[leafI];
//                     toBelow << Values[leafID];

//                     if (debug)
//                     {
//                         Pout<< " sent through "
//                             << belowID << " data for:" << leafID
//                             << " data:" << Values[leafID] << endl;
//                     }
//                 }
//             }
//         }
//     }
// }

template<class T>
void Pstream::scatterListLinear
(
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if(comm != 0){
        FatalErrorInFunction
            << "gatherListLinear only support comm == 0, but comm = " << comm << endl << flush
            << Foam::abort(FatalError);   
    }
    if (UPstream::parRun() && UPstream::nProcs() > 1)
    {
        if (Values.size() != UPstream::nProcs())
        {
            FatalErrorInFunction
                << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs()
                << Foam::abort(FatalError);
        }
        if(UPstream::master()){
            for(int belowID = firstSlave(); belowID <= lastSlave(); ++belowID){
                OPstream toBelow
                (
                    UPstream::commsTypes::scheduled,
                    belowID,
                    0,
                    tag,
                    comm
                );
                for(int proc = 0; proc < UPstream::nProcs(); ++proc){
                    if(proc != belowID){
                        toBelow << Values[proc];
                    }
                }
            }
        }else{
            IPstream fromAbove
            (
                UPstream::commsTypes::scheduled,
                masterNo(),
                0,
                tag,
                comm
            );
            for(int proc = 0; proc < UPstream::nProcs(); ++proc){
                if(proc != UPstream::myProcNo()){
                    fromAbove >> Values[proc];
                }
            }
        }
    }
}

template<class T>
void Pstream::scatterListTree
(
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if(comm != 0){
        FatalErrorInFunction
            << "gatherListLinear only support comm == 0, but comm = " << comm << endl << flush
            << Foam::abort(FatalError);   
    }
    if (UPstream::parRun() && UPstream::nProcs() > 1)
    {
        if (Values.size() != UPstream::nProcs())
        {
            FatalErrorInFunction
                << "Size of list:" << Values.size()
                << " does not equal the number of processors:"
                << UPstream::nProcs()
                << Foam::abort(FatalError);
        }

        int nLevels = 1;
        while ((1 << nLevels) < static_cast<int>(UPstream::nProcs()))
        {
            nLevels++;
        }

        for (int level = nLevels - 1; level >= 0; --level){
            int mask = (1 << level);
            if((UPstream::myProcNo() & (mask - 1)) != 0){
                continue;
            }
            int peer = UPstream::myProcNo() ^ mask;
            if(peer >= UPstream::nProcs()){
                continue;
            }
            int ignore_count = mask;
            int count = UPstream::nProcs() - ignore_count;
            if(UPstream::myProcNo() & mask){ // reciver
                // std::cout << UPstream::myProcNo() <<  " recv (" << 0 << ", " <<  UPstream::myProcNo() << "), (" << UPstream::myProcNo() + ignore_count << ", " << UPstream::nProcs() << ") from " << peer << std::endl;
                if (contiguous<T>())
                {
                    List<T> receivedValues(count);
                    UIPstream::read
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        reinterpret_cast<char*>(receivedValues.begin()),
                        receivedValues.byteSize(),
                        tag,
                        comm
                    );
                    for(int i = 0; i < UPstream::myProcNo(); ++i){
                        Values[i] = receivedValues[i];
                    }
                    for(int i = UPstream::myProcNo() + ignore_count; i < UPstream::nProcs(); ++i){
                        Values[i] = receivedValues[i - ignore_count];
                    }
                }
                else
                {
                    IPstream fromPeer
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        0,
                        tag,
                        comm
                    );
                    for(int i = 0; i < UPstream::myProcNo(); ++i){
                        fromPeer >> Values[i];
                    }
                    for(int i = UPstream::myProcNo() + ignore_count; i < UPstream::nProcs(); ++i){
                        fromPeer >> Values[i];
                    }
                }
            }else{ //sender
                // std::cout << UPstream::myProcNo() <<  " send (" << 0 << ", " <<  peer << "), (" << peer + ignore_count << ", " << UPstream::nProcs() << ") to " << peer << std::endl;
                if (contiguous<T>())
                {
                    List<T> sendingValues(count);
                    for(int i = 0; i < peer; ++i){
                        sendingValues[i] = Values[i];
                    }
                    for(int i = peer + ignore_count; i < UPstream::nProcs(); ++i){
                        sendingValues[i - ignore_count] = Values[i];
                    }
                    OPstream::write
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        reinterpret_cast<const char*>(sendingValues.begin()),
                        sendingValues.byteSize(),
                        tag,
                        comm
                    );
                }
                else
                {
                    OPstream toPeer
                    (
                        UPstream::commsTypes::scheduled,
                        peer,
                        0,
                        tag,
                        comm
                    );
                    for(int i = 0; i < peer; ++i){
                        toPeer << Values[i];
                    }
                    for(int i = peer + ignore_count; i < UPstream::nProcs(); ++i){
                        toPeer << Values[i];
                    }
                }
            }
        }
    }
}

template<class T>
void Pstream::scatterList(List<T>& Values, const int tag, const label comm)
{
    syncClockTime clock;
    // Info << "Pstream::scatterList start" << endl;
    // if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    // {
    //     scatterList(UPstream::linearCommunication(comm), Values, tag, comm);
    // }
    // else
    // {
    //     scatterList(UPstream::treeCommunication(comm), Values, tag, comm);
    // }
    scatterListTree(Values, tag, comm);
    // Info << "scatterList time : " << scatterList_end - scatterList_start << endl;
    // Info << "Pstream::scatterList end" << endl;
    Info << "Pstream::scatterList T : " << typeid(T).name() << endl;
    Info << "Pstream::scatterList time : " << clock.elapsedTime() << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //