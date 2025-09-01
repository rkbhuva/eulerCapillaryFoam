/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Implementation
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

\*---------------------------------------------------------------------------*/

#ifndef brooksCorey_H
#define brooksCorey_H

#include "capillaryTransportModelBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class brooksCorey Declaration
\*---------------------------------------------------------------------------*/

class brooksCorey
:
    public capillaryTransportModelBase
{
    // Private Data (specific to Brooks-Corey)
    dimensionedScalar residualSaturationWetting_;
    dimensionedScalar residualSaturationNonWetting_;
    dimensionedScalar lambda_;
    dimensionedScalar entryPressure_;

    // Private Member Functions
    void calculateCapillaryPressure() override;
    void calculateRelativePermeability() override;

    //- No copy construct
    brooksCorey(const brooksCorey&) = delete;

    //- No copy assignment
    void operator=(const brooksCorey&) = delete;

public:

    TypeName("brooksCorey");

    // Constructor
    brooksCorey
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~brooksCorey() = default;

    // Accessors specific to Brooks-Corey
    const dimensionedScalar& residualSaturationWetting() const { return residualSaturationWetting_; }
    const dimensionedScalar& residualSaturationNonWetting() const { return residualSaturationNonWetting_; }
    const dimensionedScalar& lambda() const { return lambda_; }
    const dimensionedScalar& entryPressure() const { return entryPressure_; }

    // I-O
    bool writeData(Ostream& os) const override;
    bool read(const dictionary& dict) override;
};

// Ostream operator (for Info<< model)
Foam::Ostream& operator<<(Foam::Ostream&, const brooksCorey&);
}

#endif
