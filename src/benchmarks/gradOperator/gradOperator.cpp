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
#include "NeoFOAM/blas/field.hpp"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/cellCentredFiniteVolume/grad/gaussGreenGrad.hpp"
#include "NeoFOAM_GPL/readers/foamMesh.hpp"
#include "Kokkos_Core.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template <typename T>
void print_field(T &a)
{
    std::cout << "a has a size of: " << a.size() << std::endl;
    auto tmp_view = Kokkos::create_mirror_view(a.field());
    Kokkos::deep_copy(tmp_view, a.field());
    for (int i = 0; i < a.size(); i++)
    {
        std::cout << "tmp_view: " << tmp_view(i) << " at: " << i << std::endl;
    }
}



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
        Temperature.apply(KOKKOS_LAMBDA(int i) { return i; });
        Foam::volScalarField Temperature_write("Temperature", T*1.0);
        auto tmp_temp = Kokkos::create_mirror_view(Temperature.field());
        Kokkos::deep_copy(tmp_temp, Temperature.field());
        forAll(Temperature_write, celli)
        {
            Temperature_write[celli] = tmp_temp(celli);
            T[celli] = celli;
        }
        T.correctBoundaryConditions();
        T.write();
        Temperature_write.write();
        NeoFOAM::gaussGreenGrad::RegisterFunctions();

        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator");
            NeoFOAM::vectorField test = NeoFOAM::gaussGreenGrad(uMesh).grad(Temperature);
            Kokkos::fence();
        }


        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator runtime selection");
            NeoFOAM::vectorField test = NeoFOAM::gaussGreenGrad(uMesh).grad(Temperature,"atomic");
            Kokkos::fence();
        }
        Foam::volVectorField gradPhiFoam("gradPhiFoam",Foam::fvc::grad(T));
        Foam::volVectorField gradPhiGPU("gradPhiGPU",Foam::fvc::grad(T));
        NeoFOAM::vectorField GPUField = NeoFOAM::gaussGreenGrad(uMesh).grad(Temperature,"atomic");
        Kokkos::fence();
        auto testview = Kokkos::create_mirror_view(GPUField.field());
        Kokkos::deep_copy(testview, GPUField.field());
        forAll(gradPhiGPU, celli)
        {
            gradPhiGPU[celli] = Foam::vector(testview(celli)(0),testview(celli)(1),testview(celli)(2));
        }
        gradPhiGPU.write();
        gradPhiFoam.write();

        auto gradGPU = NeoFOAM::gaussGreenGrad(uMesh);
        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator allocate");
            NeoFOAM::gaussGreenGrad(uMesh).grad_allocate(Temperature);
            Kokkos::fence();
        }

        NeoFOAM::vectorField gradPhi = NeoFOAM::create_gradField(uMesh.nCells());
        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator inject atomic");
            gradGPU.grad_atomic(gradPhi,Temperature);
            Kokkos::fence();
        }

        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator inject allocted vector");
            Kokkos::View<NeoFOAM::vector *> gradPhi_host("test",uMesh.nCells());
            gradGPU.grad_atomic(gradPhi,Temperature);
            Kokkos::fence();
        }
        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator inject allocted scalar[3]");
            Kokkos::View<double *[3]> gradPhi_host("test",uMesh.nCells());
            gradGPU.grad_atomic(gradPhi,Temperature);
            Kokkos::fence();
        }
        for (int i =0; i<N;i++)
        {
            addProfiling(neofoamGradOperator, "neofoamGradOperator inject allocted outside");
            NeoFOAM::vectorField gradPhi2 = NeoFOAM::create_gradField(uMesh.nCells());
            gradGPU.grad_atomic(gradPhi2,Temperature);
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
