//Flow direction 
const vector flowDir(1.0, 0.0, 0.0);

//Get reference of the polyPatch at the inlet 
const polyPatch& pPatch =  mesh.boundaryMesh()[patchName];

const PrimitivePatch<faceList, pointField> priPatch
(
    pPatch.localFaces(), 
    pPatch.localPoints()
); 

//The big mesh file is not needed anymore, 
//clear it out to save memory
mesh.clearOut();

const scalar mergeDist = 1.0e-10;
pointField mergedPoints; 
faceList mergedFaces; 
labelList pointMap;

//Merge faces and points in the master
Foam::PatchTools::gatherAndMerge
(
    mergeDist,
    priPatch,
    mergedPoints,
    mergedFaces,
    pointMap
);

pointField faceCentres; 
scalarField faceAreas; 

if(Pstream::master())
{
   //Calculate face centers and areas
   faceCentres = pointField(mergedFaces.size());
   faceAreas = scalarField(mergedFaces.size());

   forAll(mergedFaces, facei)
   {
      faceCentres[facei] = mergedFaces[facei].centre(mergedPoints);
      faceAreas[facei] = mergedFaces[facei].mag(mergedPoints);
   }

   //Sort points by distance from the first point 
   //Helps to make the CPSD matrix diagonally dominant 
   SortableList<scalar> distance(mag(faceCentres-faceCentres[0]));
   
   faceList facesCopy = mergedFaces; 
   pointField centresCopy = faceCentres; 
   scalarField areasCopy = faceAreas; 


   //Fill sorted cells;
   forAll(faceCentres, i)
   {
      faceCentres[i] = centresCopy[distance.indices()[i]];
      faceAreas[i] = areasCopy[distance.indices()[i]];
      mergedFaces[i] = facesCopy[distance.indices()[i]];
   }   
}

Pstream::scatter(faceCentres);
Pstream::scatter(faceAreas);

//Total face area
const scalar totArea = sum(faceAreas);

const label n = faceCentres.size();
Info << "\nNumber of faces on the patch = " << n << nl << nl;

//Local points on the processor
vectorField locCenters(0);

//Local face centers on the processor
labelList facesIndx(0);

//Fill the local face centers
forAll(faceCentres, i)
{
   if(i%Pstream::nProcs() == Pstream::myProcNo())
   {
      locCenters.append(faceCentres[i]);
      facesIndx.append(i);
   }
}

//Number of points on local processor 
const label nLoc = locCenters.size();

//Get mean U, I, L from the ABL profile.
vectorField Uav = prof.getUav(locCenters);
vectorField I = prof.getI(locCenters);
vectorField L = prof.getL(locCenters);

//Scale the turbulence intensity and length scales 
//by the desired factors to account for the 
//decay of turbulence downstream.
for(label i=0; i<nLoc; i++)
{
   I[i] = cmptMultiply(I[i], scaleI);
   L[i] = cmptMultiply(L[i], scaleI);
}

//Frequency step
scalar df = fMax/nFreq;	

//The smallest frequency in the simulation
scalar f0 = df/n;

//Frequencies used for interpolating the CPSD matrix.
scalarField fInterp(nInterp, 0.0);

//The interpolation method is taken from Tao, Tianyou, et al. 
// "Efficacy of Interpolation-Enhanced Schemes in Random Wind 
//  Field Simulation over Long-Span Bridges." 
//  Journal of Bridge Engineering 23.3 (2017): 0401714.
forAll(fInterp,i)
{
     fInterp[i] = f0*::pow(n*nFreq, i/(nInterp - 1.0));
}

//The range will be (f0-SMALL) to (fMax - f0 + SMALL)
//It will make all the frequencies lie between the interplation range 
fInterp[0] -= SMALL;
fInterp[nInterp - 1] += SMALL;

Info << "\nSelected interpolation frequencies = " << fInterp << nl << nl ;

//Random number generator with a seed
Random rand(seed);

//Number of components 
const label nCmpt = 3;
label cmpt = 0;

wordList cmptName(nCmpt);

cmptName[0] = "Ux";
cmptName[1] = "Uy";
cmptName[2] = "Uz";

//Double indexed phase angles sign and consing computed.
//Stored as a one dimensional array 
scalarField sinPhi(nFreq*n, 0.0);
scalarField cosPhi(nFreq*n, 0.0);

scalar twopi = 2.0*Foam::constant::mathematical::pi;

//The coherency matrix is symmetric and only the lower triangular part is stored.
//The matrix is stored in row-major packed Lapack format with n*(n+1)/2 entries 
//including the diagonal elements
scalarField Lc(n*(n+1)/2, 0.0);

//List to store the decomposed CPSD matrix for each interpolation frequency
List<scalarField>H(n*nLoc, scalarField(nInterp, 0.0));

//List to store the generated time series
List<vectorField> allU(Nt, vectorField(nLoc, vector::zero));

//Number of frequencies in the FFT 
label M(2*nFreq);

//Arrays to hold inverse FFT variables.
fftw_complex *B, *D;

B = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
D = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

//FFTW plans 
fftw_plan pln;

//Other params needed in the FFT
scalar sqrt2df = ::sqrt(2.0*df);
scalar Hjm, freq;

