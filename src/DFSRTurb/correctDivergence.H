Ut = vector::zero;

forAll(facesIndx, i)
{
    Ut[facesIndx[i]] = allU[runTime.timeIndex()][i];
}

reduce(Ut, sumOp<vectorField>());

//Free-up previous time-step 
//incase more memory is needed  
allU[runTime.timeIndex()-1].clear();

if(correctDivergence)
{ 
    Ut -= flowDir*(sum(Ut*faceAreas).x()/totArea);

    if(runTime.timeIndex() == 1)
    {
        Ut0 = Ut;
    }

    forAll(mergedFaces, celli)
    {    
        U[celli] = Ut[celli];
    }

    U.boundaryFieldRef()[0] == U.boundaryField()[0].patchInternalField();
    U.boundaryFieldRef()[1] == inletMapper.interpolate((0.5*(Ut0 + Ut))());

    //BC for Phi on the side faces. 
    ////Used to set the gradient to specific values so that wall-normal flow is suppressed
    fixedGradientFvPatchScalarField& PhiSidesFvPatch = refCast<fixedGradientFvPatchScalarField>(Phi.boundaryFieldRef()[0]);
    
    PhiSidesFvPatch.gradient() = U.boundaryField()[0] & sidesPatch.faceNormals();

    phi = fvc::flux(U);

    adjustPhi(phi, U, p);
	
    Info<< "DFSR: Applying divergence-free corrections..."  << nl;

    // Non-orthogonal velocity potential corrector loop
    while (potentialFlow.correctNonOrthogonal())
    {
        fvScalarMatrix PhiEqn
        (
            fvm::laplacian(dimensionedScalar(dimless, 1), Phi)
            ==
            fvc::div(phi)
        );

        PhiEqn.setReference(PhiRefCell, PhiRefValue);
	
        PhiEqn.solve();    

        if(potentialFlow.finalNonOrthogonalIter())
        {
            phi -= PhiEqn.flux();
        }
    }

    Info<< "Continuity error = "
    << mag(fvc::div(phi))().weightedAverage(inletMesh.V()).value()
    << endl;

    //Reconstruct velocity back from
    //corrected face fluxes 
    U = fvc::reconstruct(phi);
    
    
    forAll(mergedFaces, celli)
    {    
        Ut[celli] = U[celli];
    }
    
    Ut0 = Ut;  
}

//Add mean to the fluctuations
Ut += allUav;

