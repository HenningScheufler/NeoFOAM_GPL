# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors
add_executable(test_treeAMR "test_treeAMR.cpp")

# find_package(Kokkos REQUIRED)

set_target_properties(test_treeAMR PROPERTIES CXX_STANDARD 17)

# target_compile_definitions(test_treeAMR PUBLIC WM_LABEL_SIZE=32 NoRepository WM_DP OPENFOAM=2306 OMPI_SKIP_MPICXX)

# target_include_directories(
# test_treeAMR
#  SYSTEM
#   PUBLIC $ENV{FOAM_SRC}/finiteVolume/lnInclude
#          $ENV{FOAM_SRC}/meshTools/lnInclude
#          $ENV{FOAM_SRC}/meshTools/lnInclude
#          $ENV{FOAM_SRC}/OpenFOAM/lnInclude
#          $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude
#          ${CMAKE_SOURCE_DIR}/include
#          ${CMAKE_CURRENT_SOURCE_DIR})

# target_link_libraries(
#   test_treeAMR
#   PUBLIC $ENV{FOAM_LIBBIN}/libOpenFOAM.so
#          $ENV{FOAM_LIBBIN}/libfiniteVolume.so
#          $ENV{FOAM_LIBBIN}/libmeshTools.so
#          $ENV{FOAM_LIBBIN}/$ENV{FOAM_MPI}/libPstream.so
#          NeoFOAM
#          Kokkos::kokkos)

install(TARGETS test_treeAMR DESTINATION bin)

