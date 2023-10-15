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

#include "processorFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    procPatch_(refCast<const processorFvPatch>(p))
{}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& f
)
:
    coupledFvsPatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFvPatch>(p))
{}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isType<processorFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not processor type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const processorFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFvPatch>(p))
{
    if (!isType<processorFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const processorFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{}


template<class Type>
Foam::processorFvsPatchField<Type>::processorFvsPatchField
(
    const processorFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorFvsPatchField<Type>::~processorFvsPatchField()
{}

template<class Type>
Foam::tmp<Foam::fvsPatchField<Type>>
Foam::fvsPatchField<Type>::NewProcessorType
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
){
    return tmp<fvsPatchField<Type>>(new processorFvsPatchField<Type>(p, iF));
}

template<>
Foam::tmp<Foam::fvsPatchField<Foam::scalar>>
Foam::fvsPatchField<Foam::scalar>::NewProcessorType
(
    const fvPatch& p,
    const DimensionedField<Foam::scalar, surfaceMesh>& iF
){
    return tmp<fvsPatchField<Foam::scalar>>(new processorFvsPatchField<Foam::scalar>(p, iF));
}

template<>
Foam::tmp<Foam::fvsPatchField<Foam::Vector<Foam::scalar>>>
Foam::fvsPatchField<Foam::Vector<Foam::scalar>>::NewProcessorType
(
    const fvPatch& p,
    const DimensionedField<Foam::Vector<Foam::scalar>, surfaceMesh>& iF
){
    return tmp<fvsPatchField<Foam::Vector<Foam::scalar>>>(new processorFvsPatchField<Foam::Vector<Foam::scalar>>(p, iF));
}


template<>
Foam::tmp<Foam::fvsPatchField<Foam::SymmTensor<Foam::scalar>>>
Foam::fvsPatchField<Foam::SymmTensor<Foam::scalar>>::NewProcessorType
(
    const fvPatch& p,
    const DimensionedField<Foam::SymmTensor<Foam::scalar>, surfaceMesh>& iF
){
    return tmp<fvsPatchField<Foam::SymmTensor<Foam::scalar>>>(new processorFvsPatchField<Foam::SymmTensor<Foam::scalar>>(p, iF));
}

template<>
Foam::tmp<Foam::fvsPatchField<Foam::Tensor<Foam::scalar>>>
Foam::fvsPatchField<Foam::Tensor<Foam::scalar>>::NewProcessorType
(
    const fvPatch& p,
    const DimensionedField<Foam::Tensor<Foam::scalar>, surfaceMesh>& iF
){
    return tmp<fvsPatchField<Foam::Tensor<Foam::scalar>>>(new processorFvsPatchField<Foam::Tensor<Foam::scalar>>(p, iF));
}

template<>
Foam::tmp<Foam::fvsPatchField<Foam::SphericalTensor<Foam::scalar>>>
Foam::fvsPatchField<Foam::SphericalTensor<Foam::scalar>>::NewProcessorType
(
    const fvPatch& p,
    const DimensionedField<Foam::SphericalTensor<Foam::scalar>, surfaceMesh>& iF
){
    return tmp<fvsPatchField<Foam::SphericalTensor<Foam::scalar>>>(new processorFvsPatchField<Foam::SphericalTensor<Foam::scalar>>(p, iF));
}
// ************************************************************************* //
