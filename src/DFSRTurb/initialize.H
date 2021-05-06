//Read the face centers at the inlet and collect the to the master processor.
//Use the patch name to extract from the boundary mesh
// Collect points from all processors and combine them
List<vectorField> gatheredPoints(Pstream::nProcs());
gatheredPoints[Pstream::myProcNo()] = mesh.boundary()[patchName].Cf();
Pstream::gatherList(gatheredPoints);
Pstream::scatterList(gatheredPoints);

//Combine all the points to one list.
vectorField points = ListListOps::combine<vectorField>(gatheredPoints, accessOp<vectorField>());
// vectorField points(combinedPoints.size(), vector::zero);

// SortableList<scalar> dist(mag(combinedPoints-combinedPoints[0]));
// //dist.stor();
// forAll(points, i)
// {
//   points[i] = combinedPoints[dist.indices()[i]];
// }

// Info << dist.indices() << nl;

// // points = combinedPoints;

const label n = points.size();
Info << "\nNumber of points on the patch = " << n << nl << nl;

//Number of point on the processor. 
label nP = n/Pstream::nProcs();

//If not evenly distributed add the remaining 
//points to the last processor. 
if(Pstream::myProcNo()==Pstream::nProcs()-1)
{
  nP = nP + n - nP*Pstream::nProcs();
} 

// Pout << "\npoints: " << nP << nl << nl;

// Calculate the mean velocity, turbulence 
// intensity, and length scale profiles. 
scalarField Uav(n, 0.0);
vectorField I(n, vector::zero);
vectorField L(n, vector::zero);

//Get Uav, I, L from the ABL profile.
prof.calculateWindProfile(points, Uav, I, L);

//Scale the turbulence intensity and length scales 
//by the required factors by the scale factors.
for(label i=0; i<n; i++)
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

Info << "\nSelected interpolation frequencies = " << fInterp << nl << nl ;

// //Initialize the spectrum at any particular frequency
// //It is calculated for each point in x, y and z-direction.
// vectorField S(Np, vector::zero);

//Random number generator
Random rand(seed);

//Double indexed phase angle 
//Stored as a one dimensional array 
vectorField Phi(nFreq*n, vector::zero);

// //Double indexed frequency
// scalarField F(Nf*Np, 0.0);

//Fill the double indexed phase angles and frequencies.  
for(int i = 0; i < nFreq ; i++)
{
  for(int j = 0; j < n ; j++)
  {
    //Pick random phase angels from a uniform distribution 
    //in the interval [0.0, 2*pi]
    Phi[i*n + j]  =  2.0*Foam::constant::mathematical::pi*rand.sample01<vector>(); //Random vector
    //Fill the double indexed frequencies.
    // F[i*n + j] = (i + (j + 1.0)/n)*df;
  }
}


//List to store the decomposed CPSD matrix for each interplation frequency
List<scalarField>H(n*(n+1)/2, scalarField(nInterp, 0.0));

//List to store the generated time-series 
List<vectorField> U(Nt, vectorField(n, vector::zero));

