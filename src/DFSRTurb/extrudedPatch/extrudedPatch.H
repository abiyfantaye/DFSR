/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::extrudedPatch

Description

SourceFiles
    extrudedPatch.C
    extrudedPatchTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extrudedPatch_H
#define extrudedPatch_H

#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class extrudedPatch Declaration
\*---------------------------------------------------------------------------*/

class extrudedPatch
:
    public polyMesh
{        
    // Private Member Functions

        //- Do edge and face use points in same order?
        static bool sameOrder(const face&, const edge&);

        //- Construct and return the extruded mesh points
        template<class FaceList, class PointField>
        pointField extrudedPoints
        (
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const label nLayers,
            const scalar dx
        );

        //- Construct and return the extruded mesh faces
        template<class FaceList, class PointField>
        faceList extrudedFaces
        (
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const label nLayers
        );

        //- Construct and return the extruded mesh cells
        template<class FaceList, class PointField>
        cellList extrudedCells
        (
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const label nLayers
        );


public:

    // Constructors

        //- Construct from the primitivePatch to extrude
        template<class FaceList, class PointField>
        extrudedPatch
        (
            const IOobject&,
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const label nLayers,
            const scalar dx
        );

        //- Disallow default bitwise copy construction
        extrudedPatch(const extrudedPatch&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const extrudedPatch&) = delete;
        
        template<class FaceList, class PointField>
        List<polyPatch*> patches
        ( 
            const PrimitivePatch<FaceList, PointField>& extrudePatch,
            const label nLayers
        ); 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "extrudedPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
