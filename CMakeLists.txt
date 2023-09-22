cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_EXTENSIONS OFF)

project(nc-complex
  DESCRIPTION "Complex numbers for netCDF"
  VERSION 0.1.0
  LANGUAGES C CXX Fortran)

include(cmake/PreventInSourceBuilds.cmake)
include(ProjectOptions.cmake)

nc_complex_setup_options()
nc_complex_global_options()
nc_complex_local_options()

set(CMAKE_CXX_VISIBILITY_PRESET hidden)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(netCDF REQUIRED)

add_subdirectory(src)
if (SKBUILD)
  add_subdirectory(src/nc_complex)
endif()

if (PROJECT_IS_TOP_LEVEL)
  include(CTest)
  if (BUILD_TESTING)
    add_subdirectory(tests)
  endif()
endif()
