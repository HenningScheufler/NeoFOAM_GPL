#include "fvMesh.H"
#include "NeoFOAM/mesh/unstructuredMesh/unstructuredMesh.hpp"
#include "NeoFOAM/blas/primitives/label.hpp"

NeoFOAM::unstructuredMesh readOpenFOAMMesh(const Foam::fvMesh &mesh)
{
    const int32_t nCells = mesh.nCells();
    const int32_t nInternalFaces = mesh.nInternalFaces();

    // Create host-side views
    Kokkos::View<NeoFOAM::vector *, Kokkos::HostSpace> Sf_host("Sf_host",nInternalFaces);
    for (Foam::label facei = 0; facei < nInternalFaces; facei++)
    {
        Foam::vector sf = mesh.Sf()[facei];
        Sf_host(facei) = NeoFOAM::vector(sf[0], sf[1], sf[2]);
    }
    Kokkos::View<int32_t *, Kokkos::HostSpace> owner_host("owner_host",nInternalFaces);
    for (Foam::label facei = 0; facei < nInternalFaces; facei++)
    {
        owner_host(facei) = mesh.owner()[facei];
    }
    Kokkos::View<int32_t *, Kokkos::HostSpace> neighbour_host("neighbour_host",nInternalFaces);
    for (Foam::label facei = 0; facei < nInternalFaces; facei++)
    {
        neighbour_host(facei) = mesh.neighbour()[facei];
    }
    Kokkos::View<NeoFOAM::scalar *, Kokkos::HostSpace> V_host("V_host",nCells);
    for (Foam::label celli = 0; celli < nCells; celli++)
    {
        V_host(celli) = mesh.V()[celli];
    }

    // Create device-side views
    NeoFOAM::vectorField Sf("Sf", nInternalFaces);
    NeoFOAM::labelField  owner("owner", nInternalFaces);
    NeoFOAM::labelField  neighbour("neighbour", nInternalFaces);
    NeoFOAM::scalarField V("V", nCells);

    // Copy the data from the host to the device
    Kokkos::deep_copy(Sf.field() , Sf_host);
    Kokkos::deep_copy(owner.field() , owner_host);
    Kokkos::deep_copy(neighbour.field() , neighbour_host);
    Kokkos::deep_copy(V.field() , V_host);


    // NeoFOAM::vectorField Sf_("Sf",nInternalFaces);  // area vector
    // NeoFOAM::labelField owner_("owner",nInternalFaces);  // owner cell
    // NeoFOAM::labelField neighbour_("neighbour",nInternalFaces);  // neighbour cell
    // NeoFOAM::scalarField V_("V",nCells);  // cell volume
    // int32_t nCells_ = nCells;  // number of cells
    // int32_t nInternalFaces_ = nInternalFaces;  // number of internal faces
    NeoFOAM::unstructuredMesh uMesh(Sf, owner, neighbour, V, nCells, nInternalFaces);

    return uMesh;
}