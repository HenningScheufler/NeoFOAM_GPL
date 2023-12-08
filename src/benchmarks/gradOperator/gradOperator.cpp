/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    gradOperator


Description


\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvm.H"
#include "linear.H"
#include "uniformDimensionedFields.H"
#include "calculatedFvPatchFields.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "adjustPhi.H"
#include "findRefCell.H"
#include "IOMRFZoneList.H"
#include "constants.H"
#include "gravityMeshObject.H"

#include "columnFvMesh.H"

#include "OSspecific.H"
#include "argList.H"
#include "timeSelector.H"

#include "profiling.H"
#include "NeoFOAM/blas/fields.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/grad/gaussGreenGrad.hpp"
#include "NeoFOAM_GPL/readers/foamMesh.hpp"
#include "Kokkos_Core.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    Kokkos::initialize(argc, argv);
    {
        #include "addProfilingOption.H"
        #include "addCheckCaseOptions.H"
        #include "setRootCase.H"
        #include "createTime.H"
        #include "createMesh.H"
        #include "createFields.H"

        runTime++;
        int N = 100;

        for (int i =0; i<N;i++)
        {
            addProfiling(foamGradOperator, "foamGradOperator");
            auto test = Foam::fvc::grad(T);
        }

        NeoFOAM::unstructuredMesh uMesh = readOpenFOAMMesh(mesh);
        NeoFOAM::scalarField Temperature("T", T.internalField().size());
        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator");
            NeoFOAM::vectorField gradPhi = NeoFOAM::gaussGreenGrad(uMesh).grad(Temperature);
            // gradGPU.grad(gradPhi,Temperature);
            Kokkos::fence();
        }

        auto gradGPU = NeoFOAM::gaussGreenGrad(uMesh);
        NeoFOAM::vectorField gradPhi = create_gradField(uMesh.nCells_);
        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator inject");
            // NeoFOAM::vectorField gradPhi = NeoFOAM::gaussGreenGrad(uMesh).grad(Temperature);
            gradGPU.grad(gradPhi,Temperature);
            Kokkos::fence();
        }

        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator inject atomic");
            // NeoFOAM::vectorField gradPhi = NeoFOAM::gaussGreenGrad(uMesh).grad(Temperature);
            gradGPU.grad_atomic(gradPhi,Temperature);
            Kokkos::fence();
        }

        Foam::profiling::print(Foam::Info);
        // runTime.write();
    

        Foam::Info << "End\n"
                   << Foam::endl;
    }
    Kokkos::finalize();

    

    return 0;
}

// ************************************************************************* //
