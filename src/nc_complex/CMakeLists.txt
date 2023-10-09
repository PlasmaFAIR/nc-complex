cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

find_program(CYTHON "cython")

add_custom_command(
  OUTPUT _c_nc_complex.c
  DEPENDS _c_nc_complex.pyx
  VERBATIM
  COMMAND "${CYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/_c_nc_complex.pyx"
          --output-file "${CMAKE_CURRENT_BINARY_DIR}/_c_nc_complex.c"
          --3str
)

python_add_library(_c_nc_complex MODULE "${CMAKE_CURRENT_BINARY_DIR}/_c_nc_complex.c"
  WITH_SOABI)
target_link_libraries(_c_nc_complex PUBLIC nc_complex)

set_target_properties(_c_nc_complex PROPERTIES
  INSTALL_RPATH "$ORIGIN"
  BUILD_RPATH_USE_ORIGIN TRUE
  BUILD_WITH_INSTALL_RPATH TRUE
  )

# We make a copy of the bundled _netCDF4.pyx in order to give it a
# different name and avoid clashes with the original when importing
# (even though this is in a different namespace, because it's a
# compiled module, it'll end up with the same name for the actual .so,
# and we have fun with linking/dynamic loading)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/netCDF4/src/netCDF4/_netCDF4.pyx"
  "${CMAKE_CURRENT_BINARY_DIR}/_pf_netCDF4.pyx"
  COPYONLY
)
add_custom_command(
  OUTPUT _pf_netCDF4.c
  DEPENDS netCDF4/src/netCDF4/_netCDF4.pyx
  VERBATIM
  COMMAND "${CYTHON}" "${CMAKE_CURRENT_BINARY_DIR}/_pf_netCDF4.pyx"
          --output-file "${CMAKE_CURRENT_BINARY_DIR}/_pf_netCDF4.c"
          --include-dir "${CMAKE_CURRENT_BINARY_DIR}"
          --include-dir "${CMAKE_CURRENT_SOURCE_DIR}"
          --include-dir "${CMAKE_CURRENT_SOURCE_DIR}/netCDF4/include"
          --3str
)

python_add_library(_pf_netCDF4 MODULE
  "${CMAKE_CURRENT_BINARY_DIR}/_pf_netCDF4.c"
  WITH_SOABI)
target_include_directories(_pf_netCDF4 PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/netCDF4/include")
target_link_libraries(_pf_netCDF4 PUBLIC netCDF::netcdf PRIVATE Python::NumPy nc_complex)
set_target_properties(_pf_netCDF4 PROPERTIES
  INSTALL_RPATH "$ORIGIN"
  BUILD_RPATH_USE_ORIGIN TRUE
  BUILD_WITH_INSTALL_RPATH TRUE
  )

install(TARGETS _c_nc_complex _pf_netCDF4 DESTINATION ${SKBUILD_PROJECT_NAME})
