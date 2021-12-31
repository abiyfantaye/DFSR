/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    SRFGTurb

Description
    DFSR is a synthetic inflow turbulence generation method tailored to 
    computational wind engineering applications. The technique uses the 
    spectral representation method to generate ABL-like turbulence with 
    prescribed two-point flow statistics. Once the turbulence is generated, 
    the method performs a posteriori divergence-free operation on it so 
    that the continuity equation is satisfied. 

    Reference:
    \verbatim
        Melaku, A.F. and Bitsuamlak, G.T., 2021. 
        A divergence-free inflow turbulence generator using spectral 
        representation method for large-eddy simulation of ABL flows. 

        Journal of Wind Engineering and Industrial Aerodynamics, 212, p.104580.
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "nonOrthogonalSolutionControl.H"
#include "OFstream.H"
#include "IFstream.H"
#include "windProfile.H"
#include "extrudedPatch.H"
#include "IPstream.H"
#include "OPstream.H"
#include "fftw3.h"
#include "mkl.h"
#include "Pstream.H"
#include "PrimitivePatch.H"
#include "PatchTools.H"
#include "surfaceWriter.H"
#include "mathematicalConstants.H"
#include "fixedGradientFvPatchFields.H"
#include "pointToPointPlanarInterpolation.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

#include "functions.H"

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readDFSRTurbDict.H"
     
    //Initialize the necessary parameters.
    #include "initialize.H"

    //Generate the turbulence for each componet separately
    while(cmpt < nCmpt)
    {
        // Perform the Cholesky decomposition of the CPSD matrix.
        #include "decomposeMatrix.H"

        // Perform the time series generation.
        #include "generateInflow.H"

        cmpt++;
    }
    //Setup files to write the output 
    #include "setupWritePaths.H"    
    
    //Create temporary inlet mesh for div-free correction
    #include "createInletMesh.H"
   
    // Create fields needed for div-free correction 
    #include "createInletFields.H"

    
    while (runTime.loop())
    {    
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Applies the divergence free corrections.
        #include "correctDivergence.H"
        
        //Writes the generated turbulence 
        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<<"\nEnd "<< nl << nl;

  return 0;
}

// ************************************************************************* //


