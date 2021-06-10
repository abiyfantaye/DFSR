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

\*---------------------------------------------------------------------------*/

#include "extrudedPatch.H"
#include "wallPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
template<class FaceList, class PointField>
Foam::pointField Foam::extrudedPatch::extrudedPoints
(
    const PrimitivePatch<FaceList, PointField>& extrudePatch,
    const label nLayers,
    const scalar dx
)
{
    const pointField& surfacePoints = extrudePatch.localPoints();
    const vectorField& surfaceNormals = extrudePatch.pointNormals();

    // const label nLayers = 1;

    pointField ePoints((nLayers + 1)*surfacePoints.size());

    for (label layer=0; layer<=nLayers; layer++)
    {
        label offset = layer*surfacePoints.size();

        forAll(surfacePoints, i)
        {
            ePoints[offset + i] = surfacePoints[i] + dx*layer*surfaceNormals[i];
        }
    }

    return ePoints;
}

template<class FaceList, class PointField>
Foam::faceList Foam::extrudedPatch::extrudedFaces
(
    const PrimitivePatch<FaceList, PointField>& extrudePatch,
    const label nLayers

)
{
    const pointField& surfacePoints = extrudePatch.localPoints();
    const List<face>& surfaceFaces = extrudePatch.localFaces();
    const edgeList& surfaceEdges = extrudePatch.edges();
    const label nInternalEdges = extrudePatch.nInternalEdges();

    // const label nLayers = 1;

    label nFaces =
        (nLayers + 1)*surfaceFaces.size() + nLayers*surfaceEdges.size();

    faceList eFaces(nFaces);

    labelList quad(4);
    label facei = 0;

    // Internal faces
    for (label layer=0; layer<nLayers; layer++)
    {
        label currentLayerOffset = layer*surfacePoints.size();
        label nextLayerOffset = currentLayerOffset + surfacePoints.size();

        // Vertical faces from layer to layer+1
        for (label edgeI=0; edgeI<nInternalEdges; edgeI++)
        {
            const edge& e = surfaceEdges[edgeI];
            const labelList& edgeFaces = extrudePatch.edgeFaces()[edgeI];

            face& f = eFaces[facei++];
            f.setSize(4);

            if
            (
                (edgeFaces[0] < edgeFaces[1])
             == sameOrder(surfaceFaces[edgeFaces[0]], e)
            )
            {
                f[0] = e[0] + currentLayerOffset;
                f[1] = e[1] + currentLayerOffset;
                f[2] = e[1] + nextLayerOffset;
                f[3] = e[0] + nextLayerOffset;
            }
            else
            {
                f[0] = e[1] + currentLayerOffset;
                f[1] = e[0] + currentLayerOffset;
                f[2] = e[0] + nextLayerOffset;
                f[3] = e[1] + nextLayerOffset;
            }
        }

        // Faces between layer and layer+1
        if (layer < nLayers-1)
        {
            forAll(surfaceFaces, i)
            {
                eFaces[facei++] =
                    face
                    (
                        surfaceFaces[i] //.reverseFace()
                      + nextLayerOffset
                    );
            }
        }
    }

    // External side faces
    for (label layer=0; layer<nLayers; layer++)
    {
        label currentLayerOffset = layer*surfacePoints.size();
        label nextLayerOffset = currentLayerOffset + surfacePoints.size();

        // Side faces across layer
        for (label edgeI=nInternalEdges; edgeI<surfaceEdges.size(); edgeI++)
        {
            const edge& e = surfaceEdges[edgeI];
            const labelList& edgeFaces = extrudePatch.edgeFaces()[edgeI];

            face& f = eFaces[facei++];
            f.setSize(4);

            if (sameOrder(surfaceFaces[edgeFaces[0]], e))
            {
                f[0] = e[0] + currentLayerOffset;
                f[1] = e[1] + currentLayerOffset;
                f[2] = e[1] + nextLayerOffset;
                f[3] = e[0] + nextLayerOffset;
            }
            else
            {
                f[0] = e[1] + currentLayerOffset;
                f[1] = e[0] + currentLayerOffset;
                f[2] = e[0] + nextLayerOffset;
                f[3] = e[1] + nextLayerOffset;
            }
        }
    }

    // Bottom faces
    forAll(surfaceFaces, i)
    {
        eFaces[facei++] = face(surfaceFaces[i]).reverseFace();
    }

    // Top faces
    forAll(surfaceFaces, i)
    {
        eFaces[facei++] =
            face
            (
                surfaceFaces[i]
              + nLayers*surfacePoints.size()
            );
    }

    return eFaces;
}

template<class FaceList, class PointField>
Foam::cellList Foam::extrudedPatch::extrudedCells
(
    const PrimitivePatch<FaceList, PointField>& extrudePatch,
    const label nLayers
)
{
    const List<face>& surfaceFaces = extrudePatch.localFaces();
    const edgeList& surfaceEdges = extrudePatch.edges();
    const label nInternalEdges = extrudePatch.nInternalEdges();

    // const label nLayers = 1;

    cellList eCells(nLayers*surfaceFaces.size());


    // Size the cells
    forAll(surfaceFaces, i)
    {
        const face& f = surfaceFaces[i];

        for (label layer=0; layer<nLayers; layer++)
        {
            eCells[i + layer*surfaceFaces.size()].setSize(f.size() + 2);
        }
    }

    // Current face count per cell.
    labelList nCellFaces(eCells.size(), 0);

    label facei = 0;

    for (label layer=0; layer<nLayers; layer++)
    {
        // Side faces from layer to layer+1
        for (label i=0; i<nInternalEdges; i++)
        {
            // Get patch faces using edge
            const labelList& edgeFaces = extrudePatch.edgeFaces()[i];

            // Get cells on this layer
            label cell0 = layer*surfaceFaces.size() + edgeFaces[0];
            label cell1 = layer*surfaceFaces.size() + edgeFaces[1];

            eCells[cell0][nCellFaces[cell0]++] = facei;
            eCells[cell1][nCellFaces[cell1]++] = facei;

            facei++;
        }

        // Faces between layer and layer+1
        if (layer < nLayers-1)
        {
            forAll(surfaceFaces, i)
            {
                label cell0 = layer*surfaceFaces.size() + i;
                label cell1 = (layer+1)*surfaceFaces.size() + i;

                eCells[cell0][nCellFaces[cell0]++] = facei;
                eCells[cell1][nCellFaces[cell1]++] = facei;

                facei++;
            }
        }
    }

    // External side faces
    for (label layer=0; layer<nLayers; layer++)
    {
        // Side faces across layer
        for (label i=nInternalEdges; i<surfaceEdges.size(); i++)
        {
            // Get patch faces using edge
            const labelList& edgeFaces = extrudePatch.edgeFaces()[i];

            // Get cells on this layer
            label cell0 = layer*surfaceFaces.size() + edgeFaces[0];

            eCells[cell0][nCellFaces[cell0]++] = facei;

            facei++;
        }
    }

    // Top faces
    forAll(surfaceFaces, i)
    {
        eCells[i][nCellFaces[i]++] = facei;

        facei++;
    }

    // Bottom faces
    forAll(surfaceFaces, i)
    {
        label cell0 = (nLayers-1)*surfaceFaces.size() + i;

        eCells[cell0][nCellFaces[cell0]++] = facei;

        facei++;
    }

    return eCells;
}

template<class FaceList, class PointField>
Foam::List<Foam::polyPatch*> Foam::extrudedPatch::patches
(
    const PrimitivePatch<FaceList, PointField>& extrudePatch,
    const label nLayers

)
{
    List<polyPatch*> patches(3);

    label facei = nInternalFaces();

    label sz = nLayers*(extrudePatch.nEdges() - extrudePatch.nInternalEdges());

    patches[0] = new polyPatch
    (
        "sides",
        sz,
        facei,
        0,
        boundaryMesh(),
        polyPatch::typeName
    );

    facei += sz;

    patches[1] = new polyPatch
    (
        "inflowPatch",
        extrudePatch.size(),
        facei,
        1,
        boundaryMesh(),
        polyPatch::typeName
    );

    facei += extrudePatch.size();

    patches[2] = new polyPatch
    (
        "outflowPatch",
        extrudePatch.size(),
        facei,
        2,
        boundaryMesh(),
        polyPatch::typeName
    );


    return patches;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class FaceList, class PointField>
Foam::extrudedPatch::extrudedPatch
(
    const IOobject& io,
    const PrimitivePatch<FaceList, PointField>& extrudePatch,
    const label nLayers,
    const scalar dx
)
:
    polyMesh
    (
        io,
        extrudedPoints(extrudePatch, nLayers, dx),
        extrudedFaces(extrudePatch, nLayers),
        extrudedCells(extrudePatch, nLayers)
    )
{
    
    addPatches(patches(extrudePatch, nLayers));
}


// ************************************************************************* //
