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

#ifndef leverett_H
#define leverett_H

#include "capillaryTransportModelBase.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class leverett Declaration
\*---------------------------------------------------------------------------*/
class leverett
:
    public capillaryTransportModelBase
{
    // Private Data (specific to Leverett)
    dimensionedScalar contactAngle_;
    dimensionedScalar surfaceTension_;
    dimensionedScalar porosity_;

    // Private Member Functions
    void calculateCapillaryPressure() override;
    void calculateRelativePermeability() override;

    //- No copy construct
    leverett(const leverett&) = delete;

    //- No copy assignment
    void operator=(const leverett&) = delete;

public:

    TypeName("leverett");

    // Constructor
    leverett
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~leverett() = default;

    // Accessors specific to Leverett
    const dimensionedScalar& contactAngle() const { return contactAngle_; }
    const dimensionedScalar& surfaceTension() const { return surfaceTension_; }
    const dimensionedScalar& porosity() const { return porosity_; }

    // I-O
    bool writeData(Ostream& os) const override;
    bool read(const dictionary& dict) override;
};

// Ostream operator (for Info<< model)
Foam::Ostream& operator<<(Foam::Ostream&, const leverett&);

}

#endif
