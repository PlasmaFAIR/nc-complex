include(GenerateExportHeader)

add_library(nc_complex nc_complex.c)
add_library(plasmafair::nc_complex ALIAS nc_complex)

target_include_directories(nc_complex
  ${WARNING_GUARD}
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_compile_features(nc_complex PUBLIC c_std_99)

target_link_libraries(nc_complex PUBLIC netCDF::netcdf)

if (PROJECT_IS_TOP_LEVEL AND SKBUILD)
  install(TARGETS nc_complex DESTINATION ${SKBUILD_PROJECT_NAME})
endif()
