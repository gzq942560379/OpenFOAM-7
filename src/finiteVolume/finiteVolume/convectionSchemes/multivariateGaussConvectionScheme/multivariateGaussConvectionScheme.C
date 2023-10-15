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

\*---------------------------------------------------------------------------*/

#include "multivariateGaussConvectionScheme.H"
#include "gaussConvectionScheme.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
multivariateGaussConvectionScheme<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return gaussConvectionScheme<Type>
    (
        this->mesh(),
        faceFlux,
        tinterpScheme_()(vf)
    ).interpolate(faceFlux, vf);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
multivariateGaussConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return gaussConvectionScheme<Type>
    (
        this->mesh(),
        faceFlux,
        tinterpScheme_()(vf)
    ).flux(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type>>
multivariateGaussConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return gaussConvectionScheme<Type>
    (
        this->mesh(),
        faceFlux,
        tinterpScheme_()(vf)
    ).fvmDiv(faceFlux, vf);
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
multivariateGaussConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return gaussConvectionScheme<Type>
    (
        this->mesh(),
        faceFlux,
        tinterpScheme_()(vf)
    ).fvcDiv(faceFlux, vf);
}

template<class Type>
tmp<convectionScheme<Type>> convectionScheme<Type>::NewMultivariateGauss
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<Type>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<convectionScheme<Type>>(new multivariateGaussConvectionScheme<Type>(mesh, fields, faceFlux, schemeData));
}

template<>
tmp<convectionScheme<scalar>> convectionScheme<scalar>::NewMultivariateGauss
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<scalar>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<convectionScheme<scalar>>(new multivariateGaussConvectionScheme<scalar>(mesh, fields, faceFlux, schemeData));
}

template<>
tmp<convectionScheme<Vector<scalar>>> convectionScheme<Vector<scalar>>::NewMultivariateGauss
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<Vector<scalar>>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<convectionScheme<Vector<scalar>>>(new multivariateGaussConvectionScheme<Vector<scalar>>(mesh, fields, faceFlux, schemeData));
}

template<>
tmp<convectionScheme<Tensor<scalar>>> convectionScheme<Tensor<scalar>>::NewMultivariateGauss
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<Tensor<scalar>>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<convectionScheme<Tensor<scalar>>>(new multivariateGaussConvectionScheme<Tensor<scalar>>(mesh, fields, faceFlux, schemeData));
}

template<>
tmp<convectionScheme<SymmTensor<scalar>>> convectionScheme<SymmTensor<scalar>>::NewMultivariateGauss
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<SymmTensor<scalar>>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<convectionScheme<SymmTensor<scalar>>>(new multivariateGaussConvectionScheme<SymmTensor<scalar>>(mesh, fields, faceFlux, schemeData));
}

template<>
tmp<convectionScheme<SphericalTensor<scalar>>> convectionScheme<SphericalTensor<scalar>>::NewMultivariateGauss
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<SphericalTensor<scalar>>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
){
    return tmp<convectionScheme<SphericalTensor<scalar>>>(new multivariateGaussConvectionScheme<SphericalTensor<scalar>>(mesh, fields, faceFlux, schemeData));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
