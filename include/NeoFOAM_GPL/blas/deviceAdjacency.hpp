// SPDX-License-Identifier: MPL-2.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include <Kokkos_Core.hpp>
#include <iostream>
#include "primitives/scalar.hpp"

namespace NeoFOAM
{
    // Basically a CSR approach to connectivity - avoids 'vector of vectors' memory layout for connectivity.
    template <typename Tlabel>
    class deviceAdjacency
    {
    public:
        KOKKOS_FUNCTION
        deviceAdjacency(const deviceAdjacency<T> &rhs)
            : size_(rhs.size_), adjacency_(rhs.adjacency_), offset_(rhs.offset_)
        {

        }

        KOKKOS_FUNCTION
        deviceAdjacency(const Kokkos::View<T *> &adjacency, const Kokkos::View<T *> &offset)
            : size_(field.size()), adjacency_(adjacency_), offset_(rhs.offset)
        {

        }

        deviceAdjacency(const std::string &name, const int size)
            : size_(size), offset_(Kokkos::View<T *>(name, size)) // note adjacency not sized
        {

        }


        [[nodiscard]] inline auto data()
        {
            return {adjacency_.data(), offset_.data()};
        }

        [[nodiscard]] inline std::string name()
        {
            return offset_.name(); 
        }

        [[nodiscard]] inline int size()
        {
            return offset_.size();
        }
    
        [[nodiscard]] inline std::span<const Tlabel> at(const std::Tlabel& index) const {
            return (*this)[index];
        }

        [[nodiscard]] inline std::span<const std::Tlabel> operator[](const std::Tlabel& index) const {
            return {adjacency_.begin() + static_cast<s_size_t>(offset_[index]),
                    offset_[index + 1] - offset_[index]};
        }

    private:
        Kokkos::View<Tlabel *> adjacency_;
        Kokkos::View<Tlabel *> offset_;    // NOTE used .name here for class name

        int size_;
    };
} // namespace NeoFOAM
