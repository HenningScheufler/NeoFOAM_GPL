#include "fvMesh.H"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"

NeoFOAM::unstructuredMesh readOpenFOAMMesh(const Foam::fvMesh &mesh)
{
    const Foam::label nCells = mesh.nCells();
    const Foam::label nInternalFaces = mesh.nInternalFaces();

    // Create unmanged host-side views
    Kokkos::View<NeoFOAM::vector *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> Sf_host(mesh.Sf().primitiveFieldRef().data(), nInternalFaces);
    Kokkos::View<int32_t *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> owner_host(mesh.owner().data(), nInternalFaces);
    Kokkos::View<int32_t *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> neighbour_host(mesh.neigbour().data(), nInternalFaces);
    Kokkos::View<scalar *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> V_host(mesh.V().field().data(), nCells);

    // Create device-side views
    auto Sf_device = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), Sf_host);
    auto owner_device = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), owner_host);
    auto neighbour_device = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), neighbour_host);
    auto V_device = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), V_host);

    // Copy the data from the host to the device
    Kokkos::deep_copy(Sf_device, Sf_host);
    Kokkos::deep_copy(owner_device, owner_host);
    Kokkos::deep_copy(neighbour_device, neighbour_host);
    Kokkos::deep_copy(V_device, V_host);

    // Create the unstructuredMesh object
    NeoFOAM::unstructuredMesh uMesh(Sf_device, owner_device, neighbour_device, V_device, nCells, nInternalFaces);

    return uMesh;
}