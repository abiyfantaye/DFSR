Info<< "Reading DFSRTurbDict\n" << endl;

// IO object to read the DFSRTurb Dictionary
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
const scalar fMax(DFSRTurbDict.lookup<scalar>("fMax"));

//Time step for the inflow.
scalar dt (DFSRTurbDict.lookup<scalar>("timeStep"));

//Start time of the inflow. 
const scalar startTime (DFSRTurbDict.lookupOrDefault<scalar>("startTime", 0.0));

//End time for the inflow.
const scalar endTime (DFSRTurbDict.lookup<scalar>("endTime"));

if(dt != 1.0/(2*fMax))
{
    dt = 1.0/(2*fMax);

    Info << "Warning: timeStep must be 1/(2*fMax) to use FFT! \n" 
         << "\t timeStep set to: " << dt << nl << endl;
}

//Number of time steps, 2 is added to include begining and end times. 
const label Nt = round((endTime - startTime)/dt) + 1;

//Set run time based on DFSRTurbDict 
runTime.setTime(startTime, 0);
runTime.setDeltaT(dt);
runTime.setEndTime(endTime);

//Number of frequency segments 
const label nFreq(DFSRTurbDict.lookup<label>("nFreq"));

//Number frequencies used for the interpolation. 
const label nInterp(DFSRTurbDict.lookup<label>("nInterp"));

//Coherency decay coefficients
const tensor C(DFSRTurbDict.lookup("C"));

//A random number to seed the random number generator. 
const label seed(DFSRTurbDict.lookupOrDefault<label>("seed", 123456));

//Name of the patch.
const word patchName(DFSRTurbDict.lookup("patchName"));

//A factor to scale the turbulence intensity profiles in each direction.
const vector scaleI(DFSRTurbDict.lookupOrDefault<vector>("scaleI", vector::one));

//A factor to scale the integral length scale profile in each direction.
const vector scaleL(DFSRTurbDict.lookupOrDefault<vector>("scaleL", vector::one));

//A constant mean velocity used to define the coherency function. 
//Can be set as the mean velocity at reference target height.
const scalar cohUav(DFSRTurbDict.lookup<scalar>("cohUav"));

//If turned on, corrects the generated flow to be div-free. 
const Switch correctDivergence (DFSRTurbDict.lookupOrDefault<Switch>("correctDivergence", "on"));

//Writes sample velocity time history at the specified points if turned on. 
const Switch writeSamples (DFSRTurbDict.lookupOrDefault<Switch>("writeSamples", "off"));

//If turned on, writes the inflow data into the 'boundaryData' directory. 
const Switch writeInflow (DFSRTurbDict.lookupOrDefault<Switch>("writeInflow", "on"));

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

Info<<"\nInterpolating frequencies (nInterp): " << nInterp<<nl;

Info<<"\nCoherency decay constants(C): " << C<<nl;

Info<<"\nMean velocity used for the coherency function (cohUav): " << cohUav <<nl;

Info<<"\nRandom number generator seed: " << seed <<nl;

Info<<"\nInlet patch name: " << patchName <<nl;

Info<<"\nDivergence-free correction: " << correctDivergence <<nl;

Info<<"\nWrite sample data: " << writeSamples <<nl;

Info<<"\nWrite inflow data: " << writeInflow <<nl;


