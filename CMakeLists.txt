cmake_minimum_required(VERSION 3.20)

# Only set the cxx_standard if it is not set by someone else
if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
if (NOT DEFINED CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

project(nc-complex
  DESCRIPTION "Complex numbers for netCDF"
  VERSION 0.1.0
  LANGUAGES C CXX)

include(cmake/PreventInSourceBuilds.cmake)
include(ProjectOptions.cmake)

nc_complex_setup_options()
nc_complex_global_options()
include(Dependencies.cmake)
nc_complex_setup_dependencies()
nc_complex_local_options()

set(CMAKE_CXX_VISIBILITY_PRESET hidden)

target_compile_features(nc_complex_options INTERFACE
  cxx_std_${CMAKE_CXX_STANDARD}
  c_std_${CMAKE_C_STANDARD}
)

if (nc_complex_ENABLE_Fortran)
  enable_language(Fortran)
endif()

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(netCDF REQUIRED)

add_subdirectory(src)

if (PROJECT_IS_TOP_LEVEL)
  include(CTest)
  if (BUILD_TESTING)
    add_subdirectory(tests)
  endif()

  # MSVC doesn't really handle complex numbers in C
  if (nc_complex_BUILD_EXAMPLES AND NOT WIN32)
    add_subdirectory(examples/C)
  endif()
endif()
