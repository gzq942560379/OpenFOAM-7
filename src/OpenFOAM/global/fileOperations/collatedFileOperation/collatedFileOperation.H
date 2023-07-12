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

Class
    Foam::fileOperations::collatedFileOperation

Description
    Version of masterUncollatedFileOperation that collates regIOobjects
    into a container in the processors/ subdirectory.

    Uses threading if maxThreadFileBufferSize > 0.

See also
    masterUncollatedFileOperation

SourceFiles
    collatedFileOperation.C

\*---------------------------------------------------------------------------*/

#ifndef fileOperations_collatedFileOperation_H
#define fileOperations_collatedFileOperation_H

#include "masterUncollatedFileOperation.H"
#include "OFstreamCollator.H"
#include "fileOperationInitialise.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileOperations
{

/*---------------------------------------------------------------------------*\
                    Class collatedFileOperation Declaration
\*---------------------------------------------------------------------------*/

class collatedFileOperation
:
    public masterUncollatedFileOperation
{
protected:

    // Protected data

        //- Any communicator allocated by me
        const label myComm_;

        //- Threaded writer
        mutable OFstreamCollator writer_;

        // For non-parallel operation

            //- Number of processors (overall)
            label nProcs_;

            //- Ranks of IO handlers
            const labelList ioRanks_;


   // Private Member Functions

        static labelList ioRanks();

        //- Is proci master of communicator (in parallel) or master of
        //  the io ranks (non-parallel)
        bool isMasterRank(const label proci) const;

        //- Append to processors/ file
        bool appendObject
        (
            const regIOobject& io,
            const fileName& pathName,
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;


public:

        //- Runtime type information
        TypeName("collated");


    // Static data

        //- Max size of thread buffer size. This is the overall size of
        //  all files to be written. Starts blocking if not enough size.
        //  Read as float to enable easy specification of large sizes.
        static float maxThreadFileBufferSize;


    // Constructors

        //- Construct null
        collatedFileOperation(const bool verbose);

        //- Construct from user communicator
        collatedFileOperation
        (
            const label comm,
            const labelList& ioRanks,
            const word& typeName,
            const bool verbose
        );


    //- Destructor
    virtual ~collatedFileOperation();


    // Member Functions

        // (reg)IOobject functionality

            //- Generate disk file name for object. Opposite of filePath.
            virtual fileName objectPath
            (
                const IOobject& io,
                const word& typeName
            ) const;

            //- Writes a regIOobject (so header, contents and divider).
            //  Returns success state.
            virtual bool writeObject
            (
                const regIOobject&,
                IOstream::streamFormat format=IOstream::ASCII,
                IOstream::versionNumber version=IOstream::currentVersion,
                IOstream::compressionType compression=IOstream::UNCOMPRESSED,
                const bool write = true
            ) const;

        // Other

            //- Forcibly wait until all output done. Flush any cached data
            virtual void flush() const;

            //- Actual name of processors dir
            virtual word processorsDir(const IOobject&) const;

            //- Actual name of processors dir
            virtual word processorsDir(const fileName&) const;

            //- Set number of processor directories/results. Only used in
            //  decomposePar
            virtual void setNProcs(const label nProcs);

            template<class Type, class fileOp>
            Type masterOp
            (
                const fileName&,
                const fileOp& fop,
                const int tag,
                const label comm
            ) const;

            template<class Type, class fileOp>
            Type masterOp
            (
                const fileName&,
                const fileName&,
                const fileOp& fop,
                const int tag,
                const label comm
            ) const;

       // OSSpecific equivalents

            //- Make directory
            virtual bool mkDir(const fileName&, mode_t=0777) const;

            //- Set the file mode
            virtual bool chMod(const fileName&, const mode_t) const;

            //- Return the file mode
            virtual mode_t mode
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return the file type: directory, file or link
            virtual fileType type
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Does the name exist (as directory or file) in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool exists
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Does the name exist as a directory in the file system?
            virtual bool isDir
            (
                const fileName&,
                const bool followLink = true
            ) const;

            //- Does the name exist as a file in the file system?
            //  Optionally enable/disable check for gzip file.
            virtual bool isFile
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return size of file
            virtual off_t fileSize
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual time_t lastModified
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Return time of last file modification
            virtual double highResLastModified
            (
                const fileName&,
                const bool checkVariants = true,
                const bool followLink = true
            ) const;

            //- Read a directory and return the entries as a string list
            virtual fileNameList readDir
            (
                const fileName&,
                const fileType = fileType::file,
                const bool filterVariants = true,
                const bool followLink = true
            ) const;

            //- Copy, recursively if necessary, the source to the destination
            virtual bool cp
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = true
            ) const;

            //- Create a softlink. dst should not exist. Returns true if
            //  successful.
            virtual bool ln(const fileName& src, const fileName& dst) const;

            //- Rename src to dst
            virtual bool mv
            (
                const fileName& src,
                const fileName& dst,
                const bool followLink = false
            ) const;

            //- Rename to a corresponding backup file
            //  If the backup file already exists, attempt with
            //  "01" .. "99" suffix
            virtual bool mvBak
            (
                const fileName&,
                const std::string& ext = "bak"
            ) const;

            //- Remove a file, returning true if successful otherwise false
            virtual bool rm(const fileName&) const;

            //- Remove a directory and its contents
            virtual bool rmDir(const fileName&) const;

            //- Search for an object. checkGlobal : also check undecomposed case
            virtual fileName filePath
            (
                const bool checkGlobal,
                const IOobject&,
                const word& typeName
            ) const;

            //- Search for a directory. checkGlobal : also check undecomposed
            //  case
            virtual fileName dirPath
            (
                const bool checkGlobal,
                const IOobject&
            ) const;

            //- Read object header from supplied file
            virtual bool readHeader
            (
                IOobject&,
                const fileName&,
                const word& typeName
            ) const;
};


/*---------------------------------------------------------------------------*\
               Class collatedFileOperationInitialise Declaration
\*---------------------------------------------------------------------------*/

class collatedFileOperationInitialise
:
    public masterUncollatedFileOperationInitialise
{
public:

    // Constructors

        //- Construct from components
        collatedFileOperationInitialise(int& argc, char**& argv)
        :
            masterUncollatedFileOperationInitialise(argc, argv)
        {}


    //- Destructor
    virtual ~collatedFileOperationInitialise()
    {}


    // Member Functions

        //- Needs threading
        virtual bool needsThreading() const
        {
            return
                collatedFileOperation::maxThreadFileBufferSize
              > 0;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileOperations
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "collatedFileOperationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
