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

Application
    decomposePar

Description
    Automatically decomposes a mesh and fields of a case for parallel
    execution of OpenFOAM.

Usage
    \b decomposePar [OPTION]

    Options:
      - \par -copyUniform \n
        Copy any \a uniform directories too.

      - \par -constant

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "fvCFD.H"
#include "IOobjectList.H"
#include "domainDecomposition.H"
#include "labelIOField.H"
#include "labelFieldIOField.H"
#include "scalarIOField.H"
#include "scalarFieldIOField.H"
#include "vectorIOField.H"
#include "vectorFieldIOField.H"
#include "sphericalTensorIOField.H"
#include "sphericalTensorFieldIOField.H"
#include "symmTensorIOField.H"
#include "symmTensorFieldIOField.H"
#include "tensorIOField.H"
#include "tensorFieldIOField.H"
#include "pointFields.H"
#include "regionProperties.H"

#include "readFields.H"
#include "dimFieldDecomposer.H"
#include "fvFieldDecomposer.H"
#include "pointFieldDecomposer.H"
#include "lagrangianFieldDecomposer.H"
#include "decompositionModel.H"
#include <cassert>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    Info << "nProc : " << UPstream::nProcs() << endl;
    Info << "rank : " << UPstream::myProcNo() << endl;

    #include "createTime.H"

    Time globalTime
    (
        Time::controlDictName,
        args.rootPath(),
        args.globalCaseName()
    );

    const word& regionName = polyMesh::defaultRegion;
    const word& regionDir = Foam::regionDir(regionName);

    const IOobject dictIO
    (
        IOobject
        (
            "decomposeParDict",
            runTime.time().system(),
            regionDir, // use region if non-standard
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    label nProcs = UPstream::nProcs();
    const label nDomains = readLabel(IOdictionary(dictIO).lookup("numberOfSubdomains"));

    const_cast<fileOperation&>(fileHandler()).setNProcs(nDomains);

    if (nProcs != nDomains)
    {
        FatalErrorInFunction
            << "nProcs != nDomains" << nl
            << "nProcs : " << nProcs << nl
            << "nDomains : " << nDomains << nl
            << exit(FatalError);
    }

    Info<< "Create complete mesh" << endl;

    fvMesh mesh(
        IOobject
        (
            regionName,
            globalTime.timeName(),
            globalTime
        )
    );

    IOobjectList objects(mesh, runTime.timeName());
    PtrList<volScalarField> volScalarFields;
    readFields(mesh, objects, volScalarFields);
    PtrList<volVectorField> volVectorFields;
    readFields(mesh, objects, volVectorFields);

    fvMesh procMesh(
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime
        )
    );

    labelIOList faceProcAddressing(
        IOobject
        (
            "faceProcAddressing",
            procMesh.facesInstance(),
            procMesh.meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelIOList cellProcAddressing(
        IOobject
        (
            "cellProcAddressing",
            procMesh.facesInstance(),
            procMesh.meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    labelIOList boundaryProcAddressing(
        IOobject
        (
            "boundaryProcAddressing",
            procMesh.facesInstance(),
            procMesh.meshSubDir,
            procMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvFieldDecomposer fieldDecomposer
    (
        mesh,
        procMesh,
        faceProcAddressing,
        cellProcAddressing,
        boundaryProcAddressing
    );

    fieldDecomposer.decomposeFields(volScalarFields);
    fieldDecomposer.decomposeFields(volVectorFields);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
