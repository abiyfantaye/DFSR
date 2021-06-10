// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::windProfile

Description
    A class for reading, interpolating and calculating ABL profile. 
    
SourceFiles
    windProfile.C
\*---------------------------------------------------------------------------*/

#ifndef WINDPROFILE_H
#define WINDPROFILE_H

#include "vectorField.H"
#include "scalarField.H"
#include "List.H"
#include "interpolateXY.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class windProfile Declaration
\*---------------------------------------------------------------------------*/

class windProfile
{
	
// Private data

	//Path of the the directory where the wind profiles are stored.
    fileName path_;

    //Elevation where turbulence quantities are defined. 
    scalarField z_;

    //Mean velocity profile. 
    vectorField Uav_;

    //Turbulence intensity profile for u,v and w-components. 
    vectorField I_;
    
    //Integral length scale profiles for u,v and w-components. 
    vectorField L_;

public:
    // Constructors

        //- Construct from a path where the wind profile is stored. 
        windProfile
        (
           fileName path
        );

        windProfile
        (
        	const windProfile& profile
        );


    // Member Functions

        //- Read the wind profiles
       	void read();
        
        //- Interpolates mean velocity profiles for a given points.
        vectorField getUav
        (
        	const vectorField& p
        );


        //- Interpolates turbulence profiles for a given points.
        vectorField getI
        (
        	const vectorField& p
        );

        //- Interpolates length scale profiles for a given points.
        vectorField getL
        (
        	const vectorField& p
        );
      
        //- Returns the directory path to windProfiles.
        const fileName& path() const
        {
        	return path_;
        }
};
	

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif
