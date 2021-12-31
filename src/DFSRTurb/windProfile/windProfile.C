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

\*---------------------------------------------------------------------------*/

#include "windProfile.H"
#include "IFstream.H"
#include "IOmanip.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::windProfile::windProfile
(
	fileName path
)
:
 path_(path),
 z_(0),
 Uav_(0),
 I_(0),
 L_(0)
{
    this->read();
}

Foam::windProfile::windProfile
(
	const windProfile& prof
)
:
 path_(prof.path_),
 z_(prof.z_),
 Uav_(prof.Uav_),
 I_(prof.I_),
 L_(prof.L_)
{
}

void Foam::windProfile::read()
{   
    Info << "\nReading the wind propfile from: "<< path_ << nl;

    ifstream profFile(path_/"profile");
    
    scalar z;
    vector Uav(vector::zero);
    vector I(vector::zero);
    vector L(vector::zero);

    if(!profFile.good())
    {
    	FatalErrorInFunction
    		<< "Cannot read/find wind profile file"
    		<< "Please check your file"
    		<< exit(FatalError);
    }

    while(profFile.good())
    {
	     profFile >> z >> Uav.x() >> I.x() >> I.y() >> I.z() >> L.x() >> L.y() >> L.z();
	     z_.append(z);
	     Uav_.append(Uav);
	     I_.append(I);
	     L_.append(L);
    }

    if(Pstream::master())
    {
        //OpenFoam's 'Info' stream is sometimes crazy for keeping 
        //alignments and/or precssion so use a C-style code instead. 
        printf("\n\tz(m)\t\tU(m/s)\t\tIu\t\tIv\t\tIw\t\tLu(m)\t\tLv(m)\t\tLw(m)\n");
    
        for(label i = 0; i < z_.size(); i++)
        {
            printf("\t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f \t%.6f\n",z_[i], mag(Uav_[i]), I_[i].x(), I_[i].y(), I_[i].z(), L_[i].x(),L_[i].y(), L_[i].z()); 
        }
    }
    
    Info << nl <<nl;
}

Foam::vectorField Foam::windProfile::getUav 
(
	const vectorField& p
)
{
    scalarField z = p & vector(0.0, 0.0, 1.0);	

    return interpolateXY<vector>(z, z_, Uav_);
}


Foam::vectorField Foam::windProfile::getI 
(
	const vectorField& p
)
{
    scalarField z = p & vector(0.0, 0.0, 1.0);	

    return interpolateXY<vector>(z, z_, I_);
}

Foam::vectorField Foam::windProfile::getL 
(
	const vectorField& p
)
{
    scalarField z = p & vector(0.0, 0.0, 1.0);	
    return interpolateXY<vector>(z, z_, L_);
}

