/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "leverett.H"
#include "fvMesh.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    // Define the run-time selection entry
    defineTypeNameAndDebug(leverett, 0);
    addToRunTimeSelectionTable(capillaryTransportModelBase, leverett, dictionary);
}

// capillary pressure calculation - Leverett J-function
void Foam::leverett::calculateCapillaryPressure()
{
    forAll(pCapillary_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar sigma = surfaceTension_.value();
        scalar thetaDeg = contactAngle_.value();
        scalar thetaRad = thetaDeg * Foam::constant::mathematical::pi / 180.0;

        scalar D = D_.value();
        scalar phi = porosity_.value();

        scalar maxPcLimit = maxPc_.value();
        scalar minPcLimit = minPc_.value();

        scalar J;

        if (thetaDeg < 90.0)
        {
            J = (1.417 * alphaW) - (2.120 * pow(alphaW, 2)) + (1.263 * pow(alphaW, 3));
        }
        else
        {
            scalar oneMinusAlpha = 1.0 - alphaW;
            J = (1.417 * oneMinusAlpha) - (2.120 * pow(oneMinusAlpha, 2)) + (1.263 * pow(oneMinusAlpha, 3));
        }

        scalar Pc_calculated = (sigma * cos(thetaRad)) * sqrt(D * phi) * J;

        pCapillary_[cellI] = max(min(Pc_calculated, maxPcLimit), minPcLimit);
    }
    pCapillary_.correctBoundaryConditions();
}

// Relative Permeability Calculation (typically simple power laws)
void Foam::leverett::calculateRelativePermeability()
{
    forAll(Kr1_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];
        scalar alphanW = (wettingPhase_ == 1) ? alpha2_[cellI] : alpha1_[cellI];

        scalar Krw_val; // Wetting phase relative permeability
        scalar Krnw_val; // Non-wetting phase relative permeability

        if (alphaW <= 0.0)
        {
            Krw_val = 0.0;
            Krnw_val = 0.0;
        }
        else if (alphaW >= 1.0)
        {
            Krw_val = 1.0;
            Krnw_val = 0.0;
        }
        else
        {
            Krw_val = pow(alphaW, 3);
            Krnw_val = pow(alphanW, 3);
        }

        Krw_val = max(Krw_val, minKr1_.value());
        Krnw_val = max(Krnw_val, minKr2_.value());

        if (wettingPhase_ == 1)
        {
            Kr1_[cellI] = Krw_val;
            Kr2_[cellI] = Krnw_val;
        }
        else
        {
            Kr1_[cellI] = Krnw_val;
            Kr2_[cellI] = Krw_val;
        }
    }
    Kr1_.correctBoundaryConditions();
    Kr2_.correctBoundaryConditions();
}

Foam::leverett::leverett
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    capillaryTransportModelBase(mesh, alpha1, alpha2, dict), // constructor
    contactAngle_
    (
        "contactAngle",
        dimless,
        dict.getOrDefault<scalar>("contactAngle", 60.0)
    ),
    surfaceTension_
    (
        "surfaceTension",
        dimForce/dimLength,
        dict.getOrDefault<scalar>("surfaceTension", 0.072)
    ),
    porosity_
    (
        "porosity",
        dimless,
        dict.getOrDefault<scalar>("porosity", 0.3)
    )
{
    if (enabled_)
    {
        Info<< "Capillary transport model enabled (Leverett):" << endl;
        Info<< "    Surface tension: " << surfaceTension_.value() << " N/m" << endl;
        Info<< "    Contact Angle: " << contactAngle_.value() << " deg" << endl;
        Info<< "    Porosity: " << porosity_.value() << endl;
    }

    calculateCapillaryPressure(); // Initial calculation
    calculateRelativePermeability(); // Initial calculation
}

bool Foam::leverett::writeData(Ostream& os) const
{
    capillaryTransportModelBase::writeData(os);
    return os.good();
}

bool Foam::leverett::read(const dictionary& dict)
{
    capillaryTransportModelBase::read(dict);
    surfaceTension_.value() = dict.getOrDefault<scalar>("surfaceTension", surfaceTension_.value());
    porosity_.value() = dict.getOrDefault<scalar>("porosity", porosity_.value());
    contactAngle_.value() = dict.getOrDefault<scalar>("contactAngle", contactAngle_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::leverett& model)
{
    model.writeData(os);
    return os;
}
