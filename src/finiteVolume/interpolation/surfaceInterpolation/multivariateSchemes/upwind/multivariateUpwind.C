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
    Upwind differencing scheme class.

\*---------------------------------------------------------------------------*/

#include "multivariateUpwind.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeMultivariateSurfaceInterpolationScheme(multivariateUpwind)
}

template<class Type>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Type>>
Foam::multivariateSurfaceInterpolationScheme<Type>::NewUpwind
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Type>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<multivariateSurfaceInterpolationScheme<Type>>(
        new multivariateUpwind<Type>(mesh, vtfs, faceFlux, schemeData)
    );
}

template<>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Foam::scalar>>
Foam::multivariateSurfaceInterpolationScheme<Foam::scalar>::NewUpwind
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Foam::scalar>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<multivariateSurfaceInterpolationScheme<Foam::scalar>>(
        new multivariateUpwind<Foam::scalar>(mesh, vtfs, faceFlux, schemeData)
    );
}

template<>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Foam::Vector<Foam::scalar>>>
Foam::multivariateSurfaceInterpolationScheme<Foam::Vector<Foam::scalar>>::NewUpwind
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Foam::Vector<Foam::scalar>>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<multivariateSurfaceInterpolationScheme<Foam::Vector<Foam::scalar>>>(
        new multivariateUpwind<Foam::Vector<Foam::scalar>>(mesh, vtfs, faceFlux, schemeData)
    );
}

template<>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Foam::Tensor<Foam::scalar>>>
Foam::multivariateSurfaceInterpolationScheme<Foam::Tensor<Foam::scalar>>::NewUpwind
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Foam::Tensor<Foam::scalar>>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<multivariateSurfaceInterpolationScheme<Foam::Tensor<Foam::scalar>>>(
        new multivariateUpwind<Foam::Tensor<Foam::scalar>>(mesh, vtfs, faceFlux, schemeData)
    );
}

template<>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Foam::SymmTensor<Foam::scalar>>>
Foam::multivariateSurfaceInterpolationScheme<Foam::SymmTensor<Foam::scalar>>::NewUpwind
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Foam::SymmTensor<Foam::scalar>>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<multivariateSurfaceInterpolationScheme<Foam::SymmTensor<Foam::scalar>>>(
        new multivariateUpwind<Foam::SymmTensor<Foam::scalar>>(mesh, vtfs, faceFlux, schemeData)
    );
}

template<>
Foam::tmp<Foam::multivariateSurfaceInterpolationScheme<Foam::SphericalTensor<Foam::scalar>>>
Foam::multivariateSurfaceInterpolationScheme<Foam::SphericalTensor<Foam::scalar>>::NewUpwind
(
    const fvMesh& mesh,
    const multivariateSurfaceInterpolationScheme<Foam::SphericalTensor<Foam::scalar>>::fieldTable& vtfs,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<multivariateSurfaceInterpolationScheme<Foam::SphericalTensor<Foam::scalar>>>(
        new multivariateUpwind<Foam::SphericalTensor<Foam::scalar>>(mesh, vtfs, faceFlux, schemeData)
    );
}

// ************************************************************************* //
