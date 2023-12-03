# SPDX-License-Identifier: Unlicense
# SPDX-FileCopyrightText: 2023 NeoFOAM authors

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(NEOFOAM_GPL_INSTALL_PKGCONFIG_DIR "${CMAKE_INSTALL_FULL_LIBDIR}/pkgconfig")
set(NEOFOAM_GPL_INSTALL_CONFIG_DIR "${CMAKE_INSTALL_FULL_LIBDIR}/cmake/NeoFOAM")
set(NEOFOAM_GPL_INSTALL_MODULE_DIR "${CMAKE_INSTALL_FULL_LIBDIR}/cmake/NeoFOAM/Modules")

function(neofoam_GPL_install)
	install(TARGETS NeoFOAM EXPORT NeoFOAMTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include
)

# install the public header files
install(
DIRECTORY "${NeoFOAM_GPL_SOURCE_DIR}/include/"
DESTINATION "${CMAKE_INSTALL_FULL_INCLUDEDIR}"
FILES_MATCHING
PATTERN "*.hpp")


include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/NeoFOAM_GPL/NeoFOAMConfigVersion.cmake"
  VERSION ${NeoFOAM_VERSION}
  COMPATIBILITY AnyNewerVersion
)

export(EXPORT NeoFOAMTargets
	FILE "${CMAKE_CURRENT_BINARY_DIR}/NeoFOAMTargets/NeoFOAMTargets.cmake"
)
configure_file(cmake/NeoFOAMConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/NeoFOAM/NeoFOAMConfig.cmake"
  COPYONLY
)

set(ConfigPackageLocation lib/cmake/NeoFOAM)
install(EXPORT NeoFOAMTargets
  FILE
    NeoFOAMTargets.cmake
  DESTINATION
    ${ConfigPackageLocation}
)
install(
  FILES
    cmake/NeoFOAMConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/NeoFOAM/NeoFOAMConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
endfunction()
