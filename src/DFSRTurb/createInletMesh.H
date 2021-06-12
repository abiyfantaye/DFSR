//Reverse face normals so that the extruded mesh 
//extends in the positive x-direction 
forAll(mergedFaces, facei)
{
    mergedFaces[facei] = reverseList(mergedFaces[facei]);
}

//Create a primitive patch for the inlet 
const PrimitivePatch<faceList, pointField> inletPatch
(
    mergedFaces, 
    mergedPoints
); 

const IOobject io
(
    Foam::fvMesh::defaultRegion,
    runTime.timeName(),
    runTime,
    Foam::IOobject::MUST_READ
);

//Mean velocity at all face centers
const vectorField allUav = prof.getUav(faceCentres);

dimensionedVector Ubar("Ubar", dimVelocity, sum(faceAreas*allUav)/totArea);

//Approximate mesh size in x-direction 
const scalar dx = dt*mag(Ubar).value(); 

//Number of layers for the temp mesh
const label nLayers = 1;

//Created extruded mesh geometry
extrudedPatch infExtrude(io, inletPatch, nLayers, dx);

fvMesh inletMesh
(
    io,
    pointField(infExtrude.points()),
    faceList(infExtrude.faces()),
    cellList(infExtrude.cells()),
    false
);

//Add fvPatchs 
inletMesh.addFvPatches(infExtrude.patches(inletPatch, nLayers));

//Poly mesh at each face 
const polyPatch& sidesPatch = inletMesh.boundaryMesh()["sides"];
const polyPatch& inflowPatch = inletMesh.boundaryMesh()["inflowPatch"];

//Interplation method 
const bool nearestOnly = true;

const scalar perturbation = 1.0e-6;

// Allocate the interpolator at the inlet
pointToPointPlanarInterpolation inletMapper 
(
    faceCentres,
    inflowPatch.faceCentres(),
    perturbation,
    nearestOnly
);

nonOrthogonalSolutionControl potentialFlow(mesh, "potentialFlow");

