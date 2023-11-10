include(GenerateExportHeader)

add_library(nc_complex nc_complex.c)
add_library(plasmafair::nc_complex ALIAS nc_complex)

target_include_directories(nc_complex
  ${WARNING_GUARD}
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(nc_complex
  PUBLIC
  netCDF::netcdf
  PRIVATE
  nc_complex_warnings
  nc_complex_options
)

target_compile_features(nc_complex PUBLIC c_std_99)

set_target_properties(nc_complex PROPERTIES
  INSTALL_RPATH "$ORIGIN"
  BUILD_RPATH_USE_ORIGIN TRUE
  BUILD_WITH_INSTALL_RPATH TRUE
  )

generate_export_header(nc_complex
  EXPORT_FILE_NAME ${PROJECT_BINARY_DIR}/include/nc_complex/nc_complex_export.h)

if (NOT BUILD_SHARED_LIBS)
  target_compile_definitions(nc_complex PUBLIC NC_COMPLEX_STATIC_DEFINE)
endif()

# This target gets rerun every build, ensuring that the git version
# information is always up-to-date. The file it generates only gets
# touched if its contents are different, which means rebuilds aren't
# triggered unnecessarily
add_custom_target(
  version_header
  ${CMAKE_COMMAND} -D SRC="${CMAKE_CURRENT_LIST_DIR}/nc_complex_version.h.in"
                   -D DST="${PROJECT_BINARY_DIR}/include/nc_complex_version.h"
                   -D GIT_EXECUTABLE="${GIT_EXECUTABLE}"
                   -D GIT_DIR=".."
                   -D NC_COMPLEX_GIT_DEBUG="${NC_COMPLEX_GIT_DEBUG}"
                   -P "${PROJECT_SOURCE_DIR}/cmake/get_git_version.cmake"
  COMMENT "Generating git version information"
  )
add_dependencies(nc_complex version_header)
