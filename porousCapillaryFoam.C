/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Application
    twoPhaseEulerFoam

Group
    grpMultiphaseSolvers

Description
    Solver for a system of two compressible fluid phases with one dispersed
    phase. Eg, gas bubbles in a liquid including heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fixedValueFvsPatchFields.H"
#include "capillaryTransportModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for a system of two compressible fluid phases with one"
        " dispersed phase.\n"
        "Eg, gas bubbles in a liquid including heat-transfer."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    Info<< "Reading capillaryTransportModel dictionary\n" << endl;
    IOdictionary capillaryTransportModelDict
    (
        IOobject
        (
            "capillaryTransportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    Info<< "Creating capillaryTransportModel\n" << endl;
    autoPtr<Foam::capillaryTransportModelBase> capillaryModel
    (
        Foam::capillaryTransportModelBase::New
        (
            mesh,
            alpha1,
            alpha2,
            capillaryTransportModelDict.subDict("capillaryTransportModel")
        )
    );

    #include "createTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"

    bool faceMomentum
    (
        pimple.dict().getOrDefault("faceMomentum", false)
    );

    bool implicitPhasePressure
    (
        mesh.solverDict(alpha1.name()).getOrDefault
        (
            "implicitPhasePressure", false
        )
    );

    #include "pUf/createDDtU.H"
    #include "pU/createDDtU.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        ++runTime;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            #include "contErrs.H"

            if (capillaryModel.valid())
            {
                capillaryModel->update();

                // Print pCapillary min, max, and mean
                Info << "Capillary Pressure - min: "
                     << min(capillaryModel->pCapillary()).value() << ", max: "
                     << max(capillaryModel->pCapillary()).value() << ", mean: "
                     << capillaryModel->pCapillary().average().value() << endl;

                // Print Kr1 and Kr2 min, max, and mean
                Info << "Relative Permeability Phase1 (Kr1) - min: "
                     << min(capillaryModel->Kr1()).value() << ", max: "
                     << max(capillaryModel->Kr1()).value() << ", mean: "
                     << capillaryModel->Kr1().average().value() << endl;
                Info << "Relative Permeability Phase2 (Kr2) - min: "
                     << min(capillaryModel->Kr2()).value() << ", max: "
                     << max(capillaryModel->Kr2()).value() << ", mean: "
                     << capillaryModel->Kr2().average().value() << endl;
            }
            else
            {
                FatalErrorInFunction
                    << "capillaryModel pointer is not initialized or valid.\n"
                    << "Ensure capillaryTransportModel is properly constructed in createFields.H"
                    << abort(FatalError);
            }

            if (faceMomentum)
            {
                #include "pUf/UEqns.H"
                //#include "EEqns.H"
                #include "pUf/pEqn.H"
                #include "pUf/DDtU.H"
            }
            else
            {
                #include "pU/UEqns.H"
                //#include "EEqns.H"
                #include "pU/pEqn.H"
                #include "pU/DDtU.H"
            }

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }

        #include "write.H"

        Info << "min/max U1: " << min(U1).value() << ", "<< max(U1).value() << endl;
        Info << "min/max U2: " << min(U2).value() << ", "<< max(U2).value() << endl;
        Info << "average alpha1/alpha2: " << alpha1.average().value() << ", "
        << alpha2.average().value() << endl;

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
