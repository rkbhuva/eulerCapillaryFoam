Info<< "Constructing momentum equations" << endl;

MRF.correctBoundaryVelocity(U1);
MRF.correctBoundaryVelocity(U2);
MRF.correctBoundaryVelocity(U);

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVolume/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVolume/dimTime);

volScalarField Kd(fluid.Kd());
const dimensionedScalar D = capillaryModel->D();
const volScalarField& Kr1 = capillaryModel->Kr1();
const volScalarField& Kr2 = capillaryModel->Kr2();

{
    volScalarField Vm(fluid.Vm());

    volVectorField gradPCap = capillaryModel->capillaryMomentumSource();

    {
        U1Eqn =
        (
            fvm::ddt(alpha1, rho1, U1) + fvm::div(alphaRhoPhi1, U1)
          - fvm::Sp(contErr1, U1)
          + MRF.DDt(alpha1*rho1 + Vm, U1)
          + phase1.turbulence().divDevRhoReff(U1)
         ==
          - Vm
           *(
                fvm::ddt(U1)
              + fvm::div(phi1, U1)
              - fvm::Sp(fvc::div(phi1), U1)
              - DDtU2
            )
          + fvOptions(alpha1, rho1, U1)
        );

        U1Eqn += alpha1 * fvm::Sp(alpha1 * mu1 * D / max(Kr1, SMALL), U1);

        if (capillaryModel->wettingPhase() == 1)
        {
            U1Eqn += -alpha1 * gradPCap;
        }
        else
        {
            U1Eqn += alpha1 * gradPCap;
        }
        U1Eqn.relax();
        U1Eqn += fvm::Sp(Kd, U1);
        fvOptions.constrain(U1Eqn);
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            fvm::ddt(alpha2, rho2, U2) + fvm::div(alphaRhoPhi2, U2)
          - fvm::Sp(contErr2, U2)
          + MRF.DDt(alpha2*rho2 + Vm, U2)
          + phase2.turbulence().divDevRhoReff(U2)
         ==
          - Vm
           *(
                fvm::ddt(U2)
              + fvm::div(phi2, U2)
              - fvm::Sp(fvc::div(phi2), U2)
              - DDtU1
            )
          + fvOptions(alpha2, rho2, U2)
        );

        U2Eqn += alpha2 * fvm::Sp(alpha2 * mu2 * D / max(Kr2, SMALL), U2);
        
        /*if (capillaryModel->wettingPhase() == 2)
        {
            U2Eqn += -alpha2 * gradPCap;
        }
        else
        {
            U2Eqn += alpha2 * gradPCap;
        }*/
        U2Eqn.relax();
        U2Eqn += fvm::Sp(Kd, U2);
        fvOptions.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        fvOptions.correct(U2);
    }
}
