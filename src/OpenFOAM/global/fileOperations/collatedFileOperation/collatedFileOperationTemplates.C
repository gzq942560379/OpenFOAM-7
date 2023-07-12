/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "Pstream.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IFstream.H"
#include <mpi.h>
#include <typeinfo>
#include <cassert>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class fileOp>
Type Foam::fileOperations::collatedFileOperation::masterOp
(
    const fileName& fName,
    const fileOp& fop,
    const int tag,
    const label comm
) const
{
    if (IFstream::debug)
    {
        Pout<< "collatedFileOperation::masterOp : Operation "
            << typeid(fileOp).name()
            << " on " << fName << endl;
    }
    if (Pstream::parRun())
    {
        Type result;
        if (Pstream::master(comm))
        {
            result = fop(fName);
        }
        Pstream::scatter(result);
        return result;
    }
    else
    {
        return fop(fName);
    }
}


template<class Type, class fileOp>
Type Foam::fileOperations::collatedFileOperation::masterOp
(
    const fileName& src,
    const fileName& dest,
    const fileOp& fop,
    const int tag,
    const label comm
) const
{
    if (IFstream::debug)
    {
        Pout<< "collatedFileOperation : Operation on src:" << src
            << " dest:" << dest << endl;
    }
    if (Pstream::parRun())
    {
        FatalErrorInFunction
            << "Not support collatedFileOperation " << typeid(fileOp).name()
            << exit(FatalError);
        Type ret;
        return ret;
    }
    else
    {
        return fop(src, dest);
    }
}


// ************************************************************************* //
