// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/blas/field.hpp"
#include "scalarField.H"
#include "vectorField.H"

// namespace NeoFOAM
// {

    NeoFOAM::labelField read_labelField(std::string field_name,const Foam::labelList &field)
    {
        const int32_t size = field.size();

        // Create host-side views
        Kokkos::View<NeoFOAM::label *, Kokkos::HostSpace> field_host("field_host",size);
        for (Foam::label i = 0; i < size; i++)
        {
            field_host(i) = field[i];
        }

        // Create device-side views
        NeoFOAM::labelField nffield(field_name, size);

        // Copy the data from the host to the device
        Kokkos::deep_copy(nffield.field() , field_host);

        return nffield;
    };


    NeoFOAM::scalarField read_scalarField(std::string field_name,const Foam::scalarField &field)
    {
        const int32_t size = field.size();

        // Create host-side views
        Kokkos::View<NeoFOAM::scalar *, Kokkos::HostSpace> field_host("field_host",size);
        for (Foam::label i = 0; i < size; i++)
        {
            field_host(i) = field[i];
        }

        // Create device-side views
        NeoFOAM::scalarField nffield(field_name, size);

        // Copy the data from the host to the device
        Kokkos::deep_copy(nffield.field() , field_host);

        return nffield;
    };


    NeoFOAM::vectorField read_vectorField(std::string field_name,const Foam::vectorField &field)
    {
        const int32_t size = field.size();

        // Create host-side views
        Kokkos::View<NeoFOAM::vector *, Kokkos::HostSpace> field_host("field_host",size);
        for (Foam::label i = 0; i < size; i++)
        {
            field_host(i) = NeoFOAM::vector(field[i][0], field[i][1], field[i][2]);;
        }

        // Create device-side views
        NeoFOAM::vectorField nffield(field_name, size);

        // Copy the data from the host to the device
        Kokkos::deep_copy(nffield.field() , field_host);

        return nffield;
    };

    

// } // namespace NeoFOAM