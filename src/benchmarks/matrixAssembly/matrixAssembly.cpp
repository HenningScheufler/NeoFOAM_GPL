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
    matrixAssembly


Description


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "profiling.H"
// #include "benchmark/benchmark.h"
// #include <nanobench.h>



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// void BM_TEqnTemporalOperator(volScalarField& T) {
//     fvScalarMatrix TEqnTemporalOperator(fvm::ddt(T));
// }

// void BM_TEqnDiffusionOperator(volScalarField& Gamma,volScalarField& T) {
//     fvScalarMatrix TEqnDiffusionOperator(fvm::laplacian(Gamma, T));
// }

// void BM_TEqnConvectionOperator(surfaceScalarField& phi, volScalarField& T) {
//     fvScalarMatrix TEqnConvectionOperator(fvm::div(phi, T));
// }

// void BM_EnergyEquation(volScalarField& rho,surfaceScalarField& phi, volScalarField& Gamma,volScalarField& T) {
//     fvScalarMatrix EnergyEquation(fvm::ddt(rho, T) + fvm::div(phi, T) - fvm::laplacian(Gamma, T));
// }

// void gen(std::string const& typeName, char const* mustacheTemplate,

//     ankerl::nanobench::Bench const& bench) {


//     std::ofstream templateOut("mustache.template." + typeName);

//     templateOut << mustacheTemplate;


//     std::ofstream renderOut("mustache.render." + typeName);

//     ankerl::nanobench::render(mustacheTemplate, bench, renderOut);

// }

int main(int argc, char *argv[])
{
    #include "addProfilingOption.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    // ankerl::nanobench::Bench bench;
    // bench.title("matrixAssembly");

    // bench.run("BM_TEqnTemporalOperator", [&] {
    //     BM_TEqnTemporalOperator(T);
    // });

    // bench.run("BM_TEqnDiffusionOperator", [&] {
    //     BM_TEqnDiffusionOperator(Gamma,T);
    // });

    // bench.run("BM_TEqnConvectionOperator", [&] {
    //     BM_TEqnConvectionOperator(phi,T);
    // });

    //     bench.run("BM_EnergyEquation", [&] {
    //     BM_EnergyEquation(rho,phi,Gamma,T);
    // });

    // gen("json", ankerl::nanobench::templates::json(), bench);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
