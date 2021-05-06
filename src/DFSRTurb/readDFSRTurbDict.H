Info<< "Reading DFSRTurbDict\n" << endl;

// IO object to read the DFSRTurbDict
IOdictionary DFSRTurbDict
(
    IOobject
    (
        "DFSRTurbDict",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Maximum frequency of the simulation.
scalar fMax(DFSRTurbDict.lookup<scalar>("fMax"));

//Time step for the inflow.
scalar dt (DFSRTurbDict.lookup<scalar>("timeStep"));

//Start time of the inflow. 
scalar startTime (DFSRTurbDict.lookupOrDefault<scalar>("startTime", 0.0));

//End time for the inflow.
scalar endTime (DFSRTurbDict.lookup<scalar>("endTime"));

//Number of time steps, 2 is added to include begining and end times. 
label Nt = round((endTime - startTime)/dt) + 2;

//Number of frequency segments 
label nFreq(DFSRTurbDict.lookup<label>("nFreq"));

//Number frequencies used for the interpolation. 
label nInterp(DFSRTurbDict.lookup<label>("nInterp"));

//Coherency decay coefficients
tensor C(DFSRTurbDict.lookup("C"));

//A random number to seed the random number generator. 
label seed(DFSRTurbDict.lookupOrDefault<label>("seed", 123456));

//Name of the patch.
word patchName(DFSRTurbDict.lookup("patchName"));

//A factor to scale the turbulence intensity profiles in each direction.
vector scaleI(DFSRTurbDict.lookupOrDefault<vector>("scaleI", vector::one));

//A factor to scale the integral length scale profile in each direction.
vector scaleL(DFSRTurbDict.lookupOrDefault<vector>("scaleL", vector::one));

//A constant mean velocity used to define the coherency function. 
//Can be set as the mean velocity at reference target height.
scalar cohUav(DFSRTurbDict.lookup<scalar>("cohUav"));

//Records sample velocity time history at the specified points if turned on. 
Switch recordSamples (DFSRTurbDict.lookupOrDefault<Switch>("recordSamples", "off"));

//If turned on, corrects the velocity to keep constant flow-rate (continuity at large scale).
Switch corrMassFlux (DFSRTurbDict.lookupOrDefault<Switch>("corrMassFlux", "off"));

//ABL profile, read from "constant/boundaryData/windProfile/" folder.
windProfile prof
(
    fileName     
    (
        runTime.constant()
        /"boundaryData"
        /"windProfile"
    )
);

Info<<"\nTime step(s): "<< dt <<nl;

Info<<"\nSimulation time(s): " <<endTime<<nl;

Info<<"\nNumber of time steps: " <<Nt<<nl;

Info<<"\nMaximum frequency(fmax): " << fMax<<nl;

Info<<"\nNumber of frequency steps (nFreq): " << nFreq<<nl;

Info<<"\nInterpolatingfrequencies (nInterp): " << nInterp<<nl;

Info<<"\nCoherency decay constants(C): " << C<<nl;

Info<<"\nMean velocity used for the coherency function (cohUav): " << cohUav <<nl;

Info<<"\nRandom number generator seed: " << seed <<nl;

Info<<"\nInlet patch name: " << patchName <<nl;

// Info<<"\nNumber of points: " << C<<nl;

// Info<<"\nInlet patch extent: " << C <<nl;

// Info<<"\nGenerating coefficient matrices" << nl <<nl;
