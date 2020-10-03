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
#include "ABLProfile.H"
#include "IPstream.H"
#include "OPstream.H"
#include "fftw3.h" 
//#include <boost/algorithm/string.hpp> 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "readSRFGTurbDict.H"
     
    //Initialize the necessary parameters.
    #include "initialize.H"
    
    //Simulate the velocity field.
    #include "simulate.H"

    //Wrtie the generated velocity field to file.
    #include "writeToFile.H"

    Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    Info<<"\nEnd "<< nl << nl;

  return 0;
}


// ************************************************************************* //


