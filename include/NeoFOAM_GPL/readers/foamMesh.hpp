// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "fvMesh.H"
#include "foamFields.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/blas/primitives/label.hpp"

NeoFOAM::unstructuredMesh readOpenFOAMMesh(const Foam::fvMesh &mesh)
{
    const int32_t nCells = mesh.nCells();
    const int32_t nInternalFaces = mesh.nInternalFaces();

    Foam::scalarField magFaceAreas = mag(mesh.faceAreas());

    NeoFOAM::unstructuredMesh uMesh(
        read_vectorField("points", mesh.points()),
        read_scalarField("V", mesh.cellVolumes()),
        read_vectorField("cellCentres", mesh.cellCentres()),
        read_vectorField("Sf", mesh.faceAreas() ),
        read_vectorField("faceCentres", mesh.faceCentres()),
        read_scalarField("magFaceAreas", magFaceAreas),
        read_labelField("owner", mesh.faceOwner()),
        read_labelField("neighbour", mesh.faceNeighbour()),
        nCells,
        nInternalFaces
    );


    return uMesh;
}