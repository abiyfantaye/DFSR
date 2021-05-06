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
    Makes a box of turbulence which conforms to a given energy
    spectrum and is divergence free.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "OFstream.H"
#include "IFstream.H"
#include "windProfile.H"
#include "IPstream.H"
#include "OPstream.H"
#include "fftw3.h"
#include "mkl.h"
#include "Pstream.H"
#include "SortableList.H"
#include "LLTMatrix.H"
#include "SquareMatrix.H"
#include "mathematicalConstants.H"

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
    
    // Perform the Cholesky decomposition of the CPSD matrix.
    #include "decomposeMatrix.H"

    // Perform the time series generation.
    #include "simulate.H"
    
    // //Wrtie the generated velocity field to file.
    // #include "writeToFile.H"

    Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    Info<<"\nEnd "<< nl << nl;

  return 0;
}


// ************************************************************************* //


