/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "capillaryTransportModelBase.H"
#include "fvMesh.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    defineTypeNameAndDebug(capillaryTransportModelBase, 0);
    defineRunTimeSelectionTable(capillaryTransportModelBase, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::capillaryTransportModelBase::capillaryTransportModelBase
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    mesh_(mesh),
    alpha1_(alpha1),
    alpha2_(alpha2),
    coeffDict_(dict),
    pCapillary_
    (
        IOobject
        (
            "pCapillary",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("pCapillary", dimPressure, 0.0)
    ),
    enabled_(dict.getOrDefault<Switch>("enabled", true)),
    wettingPhase_(dict.getOrDefault<label>("wettingPhase", 1)),
    maxPc_
    (
        "maxPc", dimPressure,
        dict.getOrDefault<scalar>("maxPc", 1e5)
    ),
    minPc_
    (
        "minPc", dimPressure,
        dict.getOrDefault<scalar>("minPc", 0.0)
    ),
    minKr1_
    (
        "minKr1", dimless,
        dict.getOrDefault<scalar>("minKr1", 1e-6)
    ),
    minKr2_
    (
        "minKr2", dimless,
        dict.getOrDefault<scalar>("minKr2", 1e-6)
    ),
    D_
    (
        "D", dimensionSet(0, -2, 0, 0, 0, 0, 0), dict.getOrDefault<scalar>("D", 1e12)
    ),
    Kr1_
    (
        IOobject
        (
            "Kr1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Kr1", dimless, 1.0)
    ),
    Kr2_
    (
        IOobject
        (
            "Kr2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Kr2", dimless, 1.0)
    )
{
    if (enabled_)
    {
        Info<< "Capillary transport model base enabled." << endl;
        Info<< "    Darcy Coefficient: " << D_.value() << " m^-2" << endl;
        Info<< "    Wetting Phase: " << wettingPhase_ << endl;
    }
    else
    {
        Info<< "Capillary transport model base disabled." << endl;
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::capillaryTransportModelBase> Foam::capillaryTransportModelBase::New
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting capillary transport model " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "capillaryTransportModelBase",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<capillaryTransportModelBase>(ctorPtr(mesh, alpha1, alpha2, dict));
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::capillaryTransportModelBase::update()
{
    if (!enabled_) return;

    calculateCapillaryPressure();
    calculateRelativePermeability();
}

Foam::tmp<Foam::volVectorField>
Foam::capillaryTransportModelBase::capillaryMomentumSource() const
{
    return tmp<volVectorField>
    (
        new volVectorField
        (
            IOobject
            (
                "capillaryPressureGradient",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fvc::grad(pCapillary_)
        )
    );
}

Foam::tmp<Foam::volVectorField>
Foam::capillaryTransportModelBase::capillaryPressureGradient() const
{
    return this->capillaryMomentumSource();
}

bool Foam::capillaryTransportModelBase::writeData(Ostream& os) const
{
    return true;
}

bool Foam::capillaryTransportModelBase::read(const dictionary& dict)
{
    coeffDict_ = dict;
    enabled_ = dict.getOrDefault<Switch>("enabled", true);
    wettingPhase_ = dict.getOrDefault<label>("wettingPhase", wettingPhase_);
    D_.value() = dict.getOrDefault<scalar>("D", D_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::capillaryTransportModelBase& model)
{
    model.writeData(os);
    return os;
}