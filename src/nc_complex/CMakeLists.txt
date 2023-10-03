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
          --output-file "${CMAKE_CURRENT_BINARY_DIR}/_c_nc_complex.c")

python_add_library(_c_nc_complex MODULE "${CMAKE_CURRENT_BINARY_DIR}/_c_nc_complex.c"
  WITH_SOABI)
target_link_libraries(_c_nc_complex PUBLIC nc_complex)

set_target_properties(_c_nc_complex PROPERTIES
  INSTALL_RPATH "$ORIGIN"
  BUILD_RPATH_USE_ORIGIN TRUE
  BUILD_WITH_INSTALL_RPATH TRUE
  )

add_custom_command(
  OUTPUT _netCDF4.c
  DEPENDS netCDF4/src/netCDF4/_netCDF4.pyx
  VERBATIM
  COMMAND "${CYTHON}" "${CMAKE_CURRENT_SOURCE_DIR}/netCDF4/src/netCDF4/_netCDF4.pyx"
          --output-file "${CMAKE_CURRENT_BINARY_DIR}/_netCDF4.c"
          --include-dir "${CMAKE_CURRENT_BINARY_DIR}"
          --include-dir "${CMAKE_CURRENT_SOURCE_DIR}"
          --include-dir "${CMAKE_CURRENT_SOURCE_DIR}/netCDF4/include"
)

# TODO: work out how to configure this properly
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/constants.pyx.in"
  "${CMAKE_CURRENT_BINARY_DIR}/constants.pyx")

python_add_library(_netCDF4 MODULE "${CMAKE_CURRENT_BINARY_DIR}/_netCDF4.c"
  WITH_SOABI)
target_link_libraries(_netCDF4 PUBLIC netCDF::netcdf PRIVATE Python::NumPy nc_complex)

install(TARGETS _c_nc_complex DESTINATION ${SKBUILD_PROJECT_NAME})
