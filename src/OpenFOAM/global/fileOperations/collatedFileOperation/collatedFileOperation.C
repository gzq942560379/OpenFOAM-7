/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "collatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "threadedCollatedOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "masterOFstream.H"
#include "OFstream.H"
#include "clockTime.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(collatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        collatedFileOperation,
        word
    );

    float collatedFileOperation::maxThreadFileBufferSize
    (
        debug::floatOptimisationSwitch("maxThreadFileBufferSize", 1e9)
    );
    registerOptSwitch
    (
        "maxThreadFileBufferSize",
        float,
        collatedFileOperation::maxThreadFileBufferSize
    );

    // Mark as needing threaded mpi
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        collatedFileOperationInitialise,
        word,
        collated
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fileOperations::collatedFileOperation::ioRanks()
{
    labelList ioRanks;

    string ioRanksString(getEnv("FOAM_IORANKS"));
    if (!ioRanksString.empty())
    {
        IStringStream is(ioRanksString);
        is >> ioRanks;
    }

    return ioRanks;
}


bool Foam::fileOperations::collatedFileOperation::isMasterRank
(
    const label proci
)
const
{
    if (Pstream::parRun())
    {
        return Pstream::master(comm_);
    }
    else
    {
        // Use any IO ranks
        if (ioRanks_.size())
        {
            // Find myself in IO rank
            return findIndex(ioRanks_, proci) != -1;
        }
        else
        {
            // Assume all in single communicator
            return proci == 0;
        }
    }
}


bool Foam::fileOperations::collatedFileOperation::appendObject
(
    const regIOobject& io,
    const fileName& pathName,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    // Append to processors/ file

    label proci = detectProcessorPath(io.objectPath());

    if (debug)
    {
        Pout<< "collatedFileOperation::writeObject :"
            << " For local object : " << io.name()
            << " appending processor " << proci
            << " data to " << pathName << endl;
    }

    if (proci == -1)
    {
        FatalErrorInFunction
            << "Not a valid processor path " << pathName
            << exit(FatalError);
    }

    const bool isMaster = isMasterRank(proci);

    // Determine the local rank if the pathName is a per-rank one
    label localProci = proci;
    {
        fileName path, procDir, local;
        label groupStart, groupSize, nProcs;
        splitProcessorPath
        (
            pathName,
            path,
            procDir,
            local,
            groupStart,
            groupSize,
            nProcs
        );
        if (groupSize > 0 && groupStart != -1)
        {
            localProci = proci-groupStart;
        }
    }


    // Create string from all data to write
    string buf;
    {
        OStringStream os(fmt, ver);
        if (isMaster)
        {
            if (!io.writeHeader(os))
            {
                return false;
            }
        }

        // Write the data to the Ostream
        if (!io.writeData(os))
        {
            return false;
        }

        if (isMaster)
        {
            IOobject::writeEndDivider(os);
        }

        buf = os.str();
    }


    // Note: cannot do append + compression. This is a limitation
    // of ogzstream (or rather most compressed formats)

    OFstream os
    (
        pathName,
        IOstream::BINARY,
        ver,
        IOstream::UNCOMPRESSED, // no compression
        !isMaster
    );

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Cannot open for appending"
            << exit(FatalIOError);
    }

    if (isMaster)
    {
        IOobject::writeBanner(os)
            << "FoamFile\n{\n"
            << "    version     " << os.version() << ";\n"
            << "    format      " << os.format() << ";\n"
            << "    class       " << decomposedBlockData::typeName
            << ";\n"
            << "    location    " << pathName << ";\n"
            << "    object      " << pathName.name() << ";\n"
            << "}" << nl;
        IOobject::writeDivider(os) << nl;
    }

    // Write data
    UList<char> slice
    (
        const_cast<char*>(buf.data()),
        label(buf.size())
    );
    os << nl << "// Processor" << localProci << nl << slice << nl;

    return os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::collatedFileOperation
(
    const bool verbose
)
:
    masterUncollatedFileOperation
    (
        (
            ioRanks().size()
          ? UPstream::allocateCommunicator
            (
                UPstream::worldComm,
                subRanks(Pstream::nProcs())
            )
          : UPstream::worldComm
        ),
        false
    ),
    myComm_(comm_),
    writer_(maxThreadFileBufferSize, comm_),
    nProcs_(Pstream::nProcs()),
    ioRanks_(ioRanks())
{
    if (verbose)
    {
        InfoHeader
            << "I/O    : " << typeName
            << " (maxThreadFileBufferSize " << maxThreadFileBufferSize
            << ')' << endl;

        if (maxThreadFileBufferSize == 0)
        {
            InfoHeader
                << "         Threading not activated "
                   "since maxThreadFileBufferSize = 0." << nl
                << "         Writing may run slowly for large file sizes."
                << endl;
        }
        else
        {
            InfoHeader
                << "         Threading activated "
                   "since maxThreadFileBufferSize > 0." << nl
                << "         Requires large enough buffer to collect all data"
                    " or thread support " << nl
                << "         enabled in MPI. If thread support cannot be "
                   "enabled, deactivate" << nl
                << "         threading by setting maxThreadFileBufferSize "
                    "to 0 in" << nl
                << "         $FOAM_ETC/controlDict"
                << endl;
        }

        if (ioRanks_.size())
        {
            // Print a bit of information
            stringList ioRanks(Pstream::nProcs());
            if (Pstream::master(comm_))
            {
                ioRanks[Pstream::myProcNo()] = hostName()+"."+name(pid());
            }
            Pstream::gatherList(ioRanks);

            InfoHeader << "         IO nodes:" << endl;
            forAll(ioRanks, proci)
            {
                if (!ioRanks[proci].empty())
                {
                    InfoHeader << "             " << ioRanks[proci] << endl;
                }
            }
        }


        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::inotifyMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify" << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::timeStampMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
    }
}


Foam::fileOperations::collatedFileOperation::collatedFileOperation
(
    const label comm,
    const labelList& ioRanks,
    const word& typeName,
    const bool verbose
)
:
    masterUncollatedFileOperation(comm, false),
    myComm_(-1),
    writer_(maxThreadFileBufferSize, comm),
    nProcs_(Pstream::nProcs()),
    ioRanks_(ioRanks)
{
    if (verbose)
    {
        InfoHeader
            << "I/O    : " << typeName
            << " (maxThreadFileBufferSize " << maxThreadFileBufferSize
            << ')' << endl;

        if (maxThreadFileBufferSize == 0)
        {
            InfoHeader
                << "         Threading not activated "
                   "since maxThreadFileBufferSize = 0." << nl
                << "         Writing may run slowly for large file sizes."
                << endl;
        }
        else
        {
            InfoHeader
                << "         Threading activated "
                   "since maxThreadFileBufferSize > 0." << nl
                << "         Requires large enough buffer to collect all data"
                    " or thread support " << nl
                << "         enabled in MPI. If thread support cannot be "
                   "enabled, deactivate" << nl
                << "         threading by setting maxThreadFileBufferSize "
                    "to 0 in" << nl
                << "         $FOAM_ETC/controlDict"
                << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::inotifyMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to inotify" << endl;
        }

        if
        (
            regIOobject::fileModificationChecking
         == regIOobject::timeStampMaster
        )
        {
            WarningInFunction
                << "Resetting fileModificationChecking to timeStamp" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::~collatedFileOperation()
{
    if (myComm_ != -1 && myComm_ != UPstream::worldComm)
    {
        UPstream::freeCommunicator(myComm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::collatedFileOperation::objectPath
(
    const IOobject& io,
    const word& typeName
) const
{
    // Replacement for objectPath
    if (io.time().processorCase())
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::PROCOBJECT,
            "dummy",        // not used for processorsobject
            io.instance()
        );
    }
    else
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::OBJECT,
            word::null,
            io.instance()
        );
    }
}


bool Foam::fileOperations::collatedFileOperation::writeObject
(
    const regIOobject& io,
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool write
) const
{
    const Time& tm = io.time();
    const fileName& inst = io.instance();

    if (inst.isAbsolute() || !tm.processorCase())
    {
        mkDir(io.path());
        fileName pathName(io.objectPath());

        if (debug)
        {
            Pout<< "collatedFileOperation::writeObject :"
                << " For object : " << io.name()
                << " falling back to master-only output to " << io.path()
                << endl;
        }

        masterOFstream os
        (
            pathName,
            fmt,
            ver,
            cmp,
            false,
            write
        );

        // If any of these fail, return (leave error handling to Ostream class)
        if (!os.good())
        {
            return false;
        }
        if (!io.writeHeader(os))
        {
            return false;
        }
        // Write the data to the Ostream
        if (!io.writeData(os))
        {
            return false;
        }
        IOobject::writeEndDivider(os);

        return true;
    }
    else
    {
        // Construct the equivalent processors/ directory
        fileName path(processorsPath(io, inst, processorsDir(io)));

        mkDir(path);
        fileName pathName(path/io.name());

        if (io.global())
        {
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For global object : " << io.name()
                    << " falling back to master-only output to " << pathName
                    << endl;
            }

            masterOFstream os
            (
                pathName,
                fmt,
                ver,
                cmp,
                false,
                write
            );

            // If any of these fail, return (leave error handling to Ostream
            // class)
            if (!os.good())
            {
                return false;
            }
            if (!io.writeHeader(os))
            {
                return false;
            }
            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }
            IOobject::writeEndDivider(os);

            return true;
        }
        else if (!Pstream::parRun())
        {
            // Special path for e.g. decomposePar. Append to
            // processorsDDD/ file
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " appending to " << pathName << endl;
            }

            return appendObject(io, pathName, fmt, ver, cmp);
        }
        else
        {
            // Re-check static maxThreadFileBufferSize variable to see
            // if needs to use threading
            bool useThread = (maxThreadFileBufferSize > 0);

            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " starting collating output to " << pathName
                    << " useThread:" << useThread << endl;
            }

            if (!useThread)
            {
                writer_.waitAll();
            }

            threadedCollatedOFstream os
            (
                writer_,
                pathName,
                fmt,
                ver,
                cmp,
                useThread
            );

            // If any of these fail, return (leave error handling to Ostream
            // class)
            if (!os.good())
            {
                return false;
            }
            if (Pstream::master(comm_) && !io.writeHeader(os))
            {
                return false;
            }
            // Write the data to the Ostream
            if (!io.writeData(os))
            {
                return false;
            }
            if (Pstream::master(comm_))
            {
                IOobject::writeEndDivider(os);
            }

            return true;
        }
    }
}

void Foam::fileOperations::collatedFileOperation::flush() const
{
    if (debug)
    {
        Pout<< "collatedFileOperation::flush : clearing and waiting for thread"
            << endl;
    }
    masterUncollatedFileOperation::flush();
    // Wait for thread to finish (note: also removes thread)
    writer_.waitAll();
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const fileName& fName
) const
{
    if (Pstream::parRun())
    {
        const List<int>& procs(UPstream::procID(comm_));

        word procDir(processorsBaseDir+Foam::name(Pstream::nProcs()));

        if (procs.size() != Pstream::nProcs())
        {
            procDir +=
              + "_"
              + Foam::name(procs[0])
              + "-"
              + Foam::name(procs.last());
        }
        return procDir;
    }
    else
    {
        word procDir(processorsBaseDir+Foam::name(nProcs_));

        if (ioRanks_.size())
        {
            // Detect current processor number
            label proci = detectProcessorPath(fName);

            if (proci != -1)
            {
                // Find lowest io rank
                label minProc = 0;
                label maxProc = nProcs_-1;
                forAll(ioRanks_, i)
                {
                    if (ioRanks_[i] >= nProcs_)
                    {
                        break;
                    }
                    else if (ioRanks_[i] <= proci)
                    {
                        minProc = ioRanks_[i];
                    }
                    else
                    {
                        maxProc = ioRanks_[i]-1;
                        break;
                    }
                }
                procDir +=
                  + "_"
                  + Foam::name(minProc)
                  + "-"
                  + Foam::name(maxProc);
            }
        }

        return procDir;
    }
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const IOobject& io
) const
{
    return processorsDir(io.objectPath());
}


void Foam::fileOperations::collatedFileOperation::setNProcs(const label nProcs)
{
    nProcs_ = nProcs;

    if (debug)
    {
        Pout<< "collatedFileOperation::setNProcs :"
            << " Setting number of processors to " << nProcs_ << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fileOperations::collatedFileOperation::mkDir
(
    const fileName& dir,
    mode_t mode
) const
{
    return masterOp<mode_t, mkDirOp>
    (
        dir,
        mkDirOp(mode),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::chMod
(
    const fileName& fName,
    mode_t mode
) const
{
    return masterOp<mode_t, chModOp>
    (
        fName,
        chModOp(mode),
        Pstream::msgType(),
        comm_
    );
}


mode_t Foam::fileOperations::collatedFileOperation::mode
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<mode_t, modeOp>
    (
        fName,
        modeOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


Foam::fileType Foam::fileOperations::collatedFileOperation::type
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return fileType
    (
        masterOp<label, typeOp>
        (
            fName,
            typeOp(checkVariants, followLink),
            Pstream::msgType(),
            comm_
        )
    );
}


bool Foam::fileOperations::collatedFileOperation::exists
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<bool, existsOp>
    (
        fName,
        existsOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::isDir
(
    const fileName& fName,
    const bool followLink
) const
{
    return masterOp<bool, isDirOp>
    (
        fName,
        isDirOp(followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::isFile
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<bool, isFileOp>
    (
        fName,
        isFileOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


off_t Foam::fileOperations::collatedFileOperation::fileSize
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<off_t, fileSizeOp>
    (
        fName,
        fileSizeOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


time_t Foam::fileOperations::collatedFileOperation::lastModified
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<time_t, lastModifiedOp>
    (
        fName,
        lastModifiedOp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


double Foam::fileOperations::collatedFileOperation::highResLastModified
(
    const fileName& fName,
    const bool checkVariants,
    const bool followLink
) const
{
    return masterOp<double, lastModifiedHROp>
    (
        fName,
        lastModifiedHROp(checkVariants, followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::mvBak
(
    const fileName& fName,
    const std::string& ext
) const
{
    return masterOp<bool, mvBakOp>
    (
        fName,
        mvBakOp(ext),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::rm
(
    const fileName& fName
) const
{
    return masterOp<bool, rmOp>
    (
        fName,
        rmOp(),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::rmDir
(
    const fileName& dir
) const
{
    return masterOp<bool, rmDirOp>
    (
        dir,
        rmDirOp(),
        Pstream::msgType(),
        comm_
    );
}


Foam::fileNameList Foam::fileOperations::collatedFileOperation::readDir
(
    const fileName& dir,
    const fileType type,
    const bool filtergz,
    const bool followLink
) const
{
    return masterOp<fileNameList, readDirOp>
    (
        dir,
        readDirOp(type, filtergz, followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::cp
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool, cpOp>
    (
        src,
        dst,
        cpOp(followLink),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::ln
(
    const fileName& src,
    const fileName& dst
) const
{
    return masterOp<bool, lnOp>
    (
        src,
        dst,
        lnOp(),
        Pstream::msgType(),
        comm_
    );
}


bool Foam::fileOperations::collatedFileOperation::mv
(
    const fileName& src,
    const fileName& dst,
    const bool followLink
) const
{
    return masterOp<bool, mvOp>
    (
        src,
        dst,
        mvOp(followLink),
        Pstream::msgType(),
        comm_
    );
}


Foam::fileName Foam::fileOperations::collatedFileOperation::filePath
(
    const bool checkGlobal,
    const IOobject& io,
    const word& typeName
) const
{
    if (debug)
    {
        Pout<< "collatedFileOperation::filePath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    // Now that we have an IOobject path use it to detect & cache
    // processor directory naming
    (void)lookupProcessorsPath(io.objectPath());

    // Trigger caching of times
    (void)findTimes(io.time().path(), io.time().constant());

    // Determine master filePath and scatter

    fileName objPath;
    pathType searchType = NOTFOUND;
    word procsDir;
    word newInstancePath;

    if (Pstream::master(comm_))
    {
        // All masters search locally. Note that global objects might
        // fail (except on master). This gets handled later on (in PARENTOBJECT)
        objPath = filePathInfo
        (
            checkGlobal,
            true,
            io,
            searchType,
            procsDir,
            newInstancePath
        );

        if (debug)
        {
            Pout<< "collatedFileOperation::filePath :"
                << " master objPath:" << objPath
                << " searchType:" << fileOperation::pathTypeNames_[searchType]
                << " procsDir:" << procsDir << " instance:" << newInstancePath
                << endl;
        }

    }

    // Scatter the information about where the master found the object
    // Note: use the worldComm to make sure all processors decide
    //       the same type. Only procsDir is allowed to differ; searchType
    //       and instance have to be same
    {
        label masterType(searchType);
        Pstream::scatter(masterType);
        searchType = pathType(masterType);
    }
    Pstream::scatter(newInstancePath);

    if
    (
        checkGlobal
     || searchType == fileOperation::PARENTOBJECT
     || searchType == fileOperation::PROCBASEOBJECT
     || searchType == fileOperation::PROCBASEINSTANCE
     || io.local() == "uniform"
    )
    {
            // Distribute master path. This makes sure it is seen as uniform
            // and only gets read from the master.
            Pstream::scatter(objPath);
            Pstream::scatter(procsDir);
    }
    else
    {
        Pstream::scatter(procsDir, Pstream::msgType(), comm_);

        // Use the master type to determine if additional information is
        // needed to construct the local equivalent
        switch (searchType)
        {
            case fileOperation::PARENTOBJECT:
            case fileOperation::PROCBASEOBJECT:
            case fileOperation::PROCBASEINSTANCE:
            {
                // Already handled above
            }
            break;

            case fileOperation::ABSOLUTE:
            case fileOperation::WRITEOBJECT:
            case fileOperation::PROCUNCOLLATED:
            case fileOperation::PROCOBJECT:
            case fileOperation::FINDINSTANCE:
            case fileOperation::PROCUNCOLLATEDINSTANCE:
            case fileOperation::PROCINSTANCE:
            {
                // Construct equivalent local path
                objPath = localObjectPath
                (
                    io,
                    searchType,
                    procsDir,
                    newInstancePath
                );
            }
            break;

            case fileOperation::OBJECT:
            case fileOperation::NOTFOUND:
            {
                // Retest all processors separately since some processors might
                // have the file and some not (e.g. lagrangian data)

                objPath = masterOp<fileName, fileOrNullOp>
                (
                    io.objectPath(),
                    fileOrNullOp(true),
                    Pstream::msgType(),
                    comm_
                );
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< "collatedFileOperation::filePath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


Foam::fileName Foam::fileOperations::collatedFileOperation::dirPath
(
    const bool checkGlobal,
    const IOobject& io
) const
{
    if (debug)
    {
        Pout<< "collatedFileOperation::dirPath :"
            << " objectPath:" << io.objectPath()
            << " checkGlobal:" << checkGlobal << endl;
    }

    // Now that we have an IOobject path use it to detect & cache
    // processor directory naming
    (void)lookupProcessorsPath(io.objectPath());

    // Determine master dirPath and scatter

    fileName objPath;
    pathType searchType = NOTFOUND;
    word procsDir;
    word newInstancePath;

    if (Pstream::master(comm_))
    {
        objPath = filePathInfo
        (
            checkGlobal,
            false,
            io,
            searchType,
            procsDir,
            newInstancePath
        );
    }

    {
        label masterType(searchType);
        Pstream::scatter(masterType);   //, Pstream::msgType(), comm_);
        searchType = pathType(masterType);
    }
    Pstream::scatter(newInstancePath);  //, Pstream::msgType(), comm_);

    if
    (
        checkGlobal
     || searchType == fileOperation::PARENTOBJECT
     || searchType == fileOperation::PROCBASEOBJECT
     || searchType == fileOperation::PROCBASEINSTANCE
     || io.local() == "uniform"
    )
    {
            // Distribute master path. This makes sure it is seen as uniform
            // and only gets read from the master.
            Pstream::scatter(objPath);
            Pstream::scatter(procsDir);
    }
    else
    {
        Pstream::scatter(procsDir, Pstream::msgType(), comm_);

        // Use the master type to determine if additional information is
        // needed to construct the local equivalent
        switch (searchType)
        {
            case fileOperation::PARENTOBJECT:
            case fileOperation::PROCBASEOBJECT:
            case fileOperation::PROCBASEINSTANCE:
            {
                // Already handled above
            }
            break;

            case fileOperation::ABSOLUTE:
            case fileOperation::WRITEOBJECT:
            case fileOperation::PROCUNCOLLATED:
            case fileOperation::PROCOBJECT:
            case fileOperation::FINDINSTANCE:
            case fileOperation::PROCUNCOLLATEDINSTANCE:
            case fileOperation::PROCINSTANCE:
            {
                // Construct equivalent local path
                objPath = localObjectPath
                (
                    io,
                    searchType,
                    procsDir,
                    newInstancePath
                );
            }
            break;

            case fileOperation::OBJECT:
            case fileOperation::NOTFOUND:
            {
                // Retest all processors separately since some processors might
                // have the file and some not (e.g. lagrangian data)
                objPath = masterOp<fileName, fileOrNullOp>
                (
                    io.objectPath(),
                    fileOrNullOp(false),
                    Pstream::msgType(),
                    comm_
                );
            }
            break;
        }
    }

    if (debug)
    {
        Pout<< "collatedFileOperation::dirPath :"
            << " Returning from file searching:" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    filePath  :" << objPath << endl << endl;
    }
    return objPath;
}


bool Foam::fileOperations::collatedFileOperation::readHeader
(
    IOobject& io,
    const fileName& fName,
    const word& typeName
) const
{
    bool ok = false;

    if (debug)
    {
        Pout<< "collatedFileOperation::readHeader :" << endl
            << "    objectPath:" << io.objectPath() << endl
            << "    fName     :" << fName << endl;
    }

    // Get filePaths on world master
    // fileNameList filePaths(Pstream::nProcs(Pstream::worldComm));
    // filePaths[Pstream::myProcNo(Pstream::worldComm)] = fName;
    // Pstream::gatherList(filePaths, Pstream::msgType(), Pstream::worldComm);
    // bool uniform = uniformFile(filePaths);
    // Pstream::scatter(uniform, Pstream::msgType(), Pstream::worldComm);

    // fileName fNameMaster;
    // if(Pstream::master(Pstream::worldComm)){
    //     fNameMaster = fName;
    // }
    // Pstream::scatter(fNameMaster);
    // bool uniform = (fNameMaster == fName);
    // reduce(uniform, andOp<bool>());

    bool uniform = true;


    if (uniform)
    {
        if (Pstream::master(Pstream::worldComm))
        {
            if (!fName.empty())
            {
                IFstream is(fName);

                if (is.good())
                {
                    ok = io.readHeader(is);
                    if (io.headerClassName() == decomposedBlockData::typeName)
                    {
                        // Read the header inside the container (master data)
                        ok = decomposedBlockData::readMasterHeader(io, is);
                    }
                }
            }
        }
        Pstream::scatter(ok, Pstream::msgType(), Pstream::worldComm);
        Pstream::scatter
        (
            io.headerClassName(),
            Pstream::msgType(),
            Pstream::worldComm
        );
        Pstream::scatter(io.note(), Pstream::msgType(), Pstream::worldComm);
    }
    else
    {
        fileNameList filePaths;
        if (Pstream::nProcs(comm_) != Pstream::nProcs(Pstream::worldComm))
        {
            // Re-gather file paths on local master
            filePaths.setSize(Pstream::nProcs(comm_));
            filePaths[Pstream::myProcNo(comm_)] = fName;
            Pstream::gatherList(filePaths, Pstream::msgType(), comm_);
        }

        boolList result(Pstream::nProcs(comm_), false);
        wordList headerClassName(Pstream::nProcs(comm_));
        stringList note(Pstream::nProcs(comm_));
        if (Pstream::master(comm_))
        {
            forAll(filePaths, proci)
            {
                if (!filePaths[proci].empty())
                {
                    if (proci > 0 && filePaths[proci] == filePaths[proci-1])
                    {
                        result[proci] = result[proci-1];
                        headerClassName[proci] = headerClassName[proci-1];
                        note[proci] = note[proci-1];
                    }
                    else
                    {
                        IFstream is(filePaths[proci]);

                        if (is.good())
                        {
                            result[proci] = io.readHeader(is);
                            if
                            (
                                io.headerClassName()
                             == decomposedBlockData::typeName
                            )
                            {
                                // Read the header inside the container (master
                                // data)
                                result[proci] = decomposedBlockData::
                                readMasterHeader
                                (
                                    io,
                                    is
                                );
                            }
                            headerClassName[proci] = io.headerClassName();
                            note[proci] = io.note();
                        }
                    }
                }
            }
        }
        ok = scatterList(result, Pstream::msgType(), comm_);
        io.headerClassName() = scatterList
        (
            headerClassName,
            Pstream::msgType(),
            comm_
        );
        io.note() = scatterList(note, Pstream::msgType(), comm_);
    }

    if (debug)
    {
        Pout<< "collatedFileOperation::readHeader :" << " ok:" << ok
            << " class:" << io.headerClassName() << endl;
    }
    return ok;
}
// ************************************************************************* //
