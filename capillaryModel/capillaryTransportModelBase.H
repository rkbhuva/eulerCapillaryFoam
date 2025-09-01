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

#ifndef capillaryTransportModelBase_H
#define capillaryTransportModelBase_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "dimensionedScalar.H"
#include "dimensionedVector.H"
#include "Switch.H"
#include "IOobject.H"
#include "autoPtr.H"
#include "className.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class capillaryTransportModelBase
{
protected:
    // Protected Data (common to all models)
    const fvMesh& mesh_;
    const volScalarField& alpha1_;
    const volScalarField& alpha2_;
    dictionary coeffDict_;

    volScalarField pCapillary_;
    Switch enabled_;
    label wettingPhase_;

    dimensionedScalar maxPc_;
    dimensionedScalar minPc_;
    dimensionedScalar minKr1_;
    dimensionedScalar minKr2_;

    dimensionedScalar D_;
    volScalarField Kr1_;
    volScalarField Kr2_;

private:
    //- No copy construct
    capillaryTransportModelBase(const capillaryTransportModelBase&) = delete;

    //- No copy assignment
    void operator=(const capillaryTransportModelBase&) = delete;

public:

    TypeName("capillaryTransportModelBase");

    declareRunTimeSelectionTable
    (
        autoPtr,
        capillaryTransportModelBase,
        dictionary,
        (
            const fvMesh& mesh,
            const volScalarField& alpha1,
            const volScalarField& alpha2,
            const dictionary& dict
        ),
        (mesh, alpha1, alpha2, dict)
    );

    // Constructor
    capillaryTransportModelBase
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    virtual ~capillaryTransportModelBase() = default;

    // Selector
    static autoPtr<capillaryTransportModelBase> New
    (
        const fvMesh& mesh,
        const volScalarField& alpha1,
        const volScalarField& alpha2,
        const dictionary& dict
    );

    // Accessors (common to all models)
    const volScalarField& pCapillary() const { return pCapillary_; }
    bool enabled() const { return enabled_; }
    label wettingPhase() const { return wettingPhase_; }
    const dimensionedScalar& D() const { return D_; }
    const volScalarField& Kr1() const { return Kr1_; }
    const volScalarField& Kr2() const { return Kr2_; }
    const dimensionedScalar& maxPc() const { return maxPc_; }
    const dimensionedScalar& minPc() const { return minPc_; }
    const dimensionedScalar& minKr1() const { return minKr1_; }
    const dimensionedScalar& minKr2() const { return minKr2_; }

    virtual void calculateCapillaryPressure() = 0;
    virtual void calculateRelativePermeability() = 0;
    virtual void update();

    virtual tmp<volVectorField> capillaryMomentumSource() const;
    virtual tmp<volVectorField> capillaryPressureGradient() const;
    virtual bool writeData(Ostream& os) const;
    virtual bool read(const dictionary& dict);
};

// Ostream operator
Foam::Ostream& operator<<(Foam::Ostream&, const capillaryTransportModelBase&);

}

#endif
