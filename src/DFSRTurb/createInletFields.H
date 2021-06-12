//Create fields at the inlet needed for div-free operation 
//Boundary conditions for Phi and U
wordList PhiBCTypes(3); 
wordList UBCTypes(3); 

//Boundary condition for Phi and 
PhiBCTypes[0] = "fixedGradient"; //sides dl/dx = 0, dl/dy = -Uy, dl/dz = -Uz
PhiBCTypes[1] = "zeroGradient";  //inlet 
PhiBCTypes[2] = "fixedValue";  //outlet

UBCTypes[0] = "slip"; //sides dl/dx = 0, dl/dy = -Uy, dl/dz = -Uz
UBCTypes[1] = "fixedValue";  //inlet 
UBCTypes[2] = "zeroGradient";  //outlet

//Constructing velocity potential field Phi
volScalarField Phi
(
    IOobject
    (
        "Phi",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    inletMesh,
    dimensionedScalar("Phi", dimLength*dimVelocity, 0.0), 
    PhiBCTypes
);

//Instantaneous velocity on all the cell/face centers 
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    inletMesh,
    dimensionedVector("U", dimVelocity, vector::zero),
    UBCTypes
); 

//Creating face fluxes.
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(U)
); 

//Create pressure field 
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    inletMesh,
    dimensionedScalar("p", sqr(dimVelocity), 0.0), 
    PhiBCTypes
);

label PhiRefCell = -Pstream::myProcNo();
scalar PhiRefValue = 0;


//Instantaneous velocity at time t and t_(n-1)
vectorField Ut(n, vector::zero);
vectorField Ut0(n, vector::zero);
