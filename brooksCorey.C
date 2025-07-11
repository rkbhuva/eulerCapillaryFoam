/*---------------------------------------------------------------------------*\
  =========                 |\
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |\
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |\
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
\*---------------------------------------------------------------------------*/

#include "brooksCorey.H"
#include "fvMesh.H"
#include "fvc.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"

namespace Foam
{
    // Define the run-time selection entry
    defineTypeNameAndDebug(brooksCorey, 0);
    addToRunTimeSelectionTable(capillaryTransportModelBase, brooksCorey, dictionary);
}

// capillary pressure calculation - Brooks Corey
void Foam::brooksCorey::calculateCapillaryPressure()
{
    forAll(pCapillary_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar Srw = residualSaturationWetting_.value();
        scalar Srnw = residualSaturationNonWetting_.value();
        scalar lambda = max(lambda_.value(), SMALL);
        scalar Pd = entryPressure_.value();

        scalar Se = (alphaW - Srw) / (1.0 - Srw - Srnw + SMALL);
        Se = max(1e-12, min(1.0 - 1e-12, Se));

        scalar maxPcLimit = maxPc_.value();
        scalar minPcLimit = minPc_.value();

        scalar Pc_calculated;

        if (Se > 1.0 - 1e-9)
        {
            Pc_calculated = 0.0;
        }
        else if (Se < 1e-9)
        {
            Pc_calculated = maxPcLimit;
        }
        else
        {
            Pc_calculated = Pd * pow(Se, -1.0/lambda);
        }

        pCapillary_[cellI] = max(min(Pc_calculated, maxPcLimit), minPcLimit);
    }
    pCapillary_.correctBoundaryConditions();
}

// Relative Permeability Calculation - Brooks Corey Model
void Foam::brooksCorey::calculateRelativePermeability()
{
    forAll(Kr1_, cellI)
    {
        scalar alphaW = (wettingPhase_ == 1) ? alpha1_[cellI] : alpha2_[cellI];

        scalar Srw = residualSaturationWetting_.value();
        scalar Srnw = residualSaturationNonWetting_.value();
        scalar lambda = max(lambda_.value(), SMALL);

        scalar Se = (alphaW - Srw) / (1.0 - Srw - Srnw + SMALL);
        Se = max(1e-12, min(1.0 - 1e-12, Se));

        scalar Krw_val; // Wetting phase relative permeability
        scalar Krnw_val; // Non-wetting phase relative permeability

        if (Se <= 0.0)
        {
            Krw_val = 0.0;
            Krnw_val = 0.0;
        }
        else if (Se >= 1.0)
        {
            Krw_val = 1.0;
            Krnw_val = 0.0;
        }
        else
        {
            Krw_val = pow(Se, (2.0 + 3.0 * lambda) / lambda);
            Krnw_val = pow(1.0 - Se, 2.0) * (1.0 - pow(Se, (2.0 + lambda) / lambda));
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

Foam::brooksCorey::brooksCorey
(
    const fvMesh& mesh,
    const volScalarField& alpha1,
    const volScalarField& alpha2,
    const dictionary& dict
)
:
    capillaryTransportModelBase(mesh, alpha1, alpha2, dict), // constructor
    residualSaturationWetting_
    (
        "residualSaturationWetting", dimless, dict.getOrDefault<scalar>("residualSaturationWetting", 0.05)
    ),
    residualSaturationNonWetting_
    (
        "residualSaturationNonWetting", dimless, dict.getOrDefault<scalar>("residualSaturationNonWetting", 0.1)
    ),
    lambda_
    (
        "lambda", dimless, dict.getOrDefault<scalar>("lambda", 2.0)
    ),
    entryPressure_
    (
        "entryPressure", dimPressure, dict.getOrDefault<scalar>("entryPressure", 1.0)
    )
{
    if (enabled_)
    {
        Info<< "Capillary transport model enabled (Brooks-Corey):" << endl;
        Info<< "    Residual Saturation (Wetting): " << residualSaturationWetting_.value() << endl;
        Info<< "    Residual Saturation (Non-Wetting): " << residualSaturationNonWetting_.value() << endl;
        Info<< "    Lambda (Pore Size Distribution Index): " << lambda_.value() << endl;
        Info<< "    Entry (Displacement) Pressure: " << entryPressure_.value() << " Pa" << endl;
    }

    calculateCapillaryPressure(); // Initial calculation
    calculateRelativePermeability(); // Initial calculation
}

bool Foam::brooksCorey::writeData(Ostream& os) const
{
    capillaryTransportModelBase::writeData(os);
    return os.good();
}

bool Foam::brooksCorey::read(const dictionary& dict)
{
    capillaryTransportModelBase::read(dict);
    residualSaturationWetting_.value() = dict.getOrDefault<scalar>("residualSaturationWetting", residualSaturationWetting_.value());
    residualSaturationNonWetting_.value() = dict.getOrDefault<scalar>("residualSaturationNonWetting", residualSaturationNonWetting_.value());
    lambda_.value() = dict.getOrDefault<scalar>("lambda", lambda_.value());
    entryPressure_.value() = dict.getOrDefault<scalar>("entryPressure", entryPressure_.value());
    return true;
}

Foam::Ostream& operator<<(Foam::Ostream& os, const Foam::brooksCorey& model)
{
    model.writeData(os);
    return os;
}