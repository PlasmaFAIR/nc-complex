# Call nx-config with an argument, and append the resulting path to a list
# Taken from https://github.com/LiamBindle/geos-chem/blob/feature/CMake/CMakeScripts/FindNetCDF.cmake
function(inspect_netcdf_config VAR NX_CONFIG ARG)
  execute_process(
    COMMAND ${NX_CONFIG} ${ARG}
    OUTPUT_VARIABLE NX_CONFIG_OUTPUT
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if (NX_CONFIG_OUTPUT)
    set(${VAR} ${NX_CONFIG_OUTPUT} PARENT_SCOPE)
  endif()
endfunction()
