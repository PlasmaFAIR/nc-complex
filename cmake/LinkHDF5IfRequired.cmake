include(CheckCSourceCompiles)

# Link `TARGET` against HDF5 if required
#
# This is necessary because the netCDF4 Python API uses
# `H5get_libversion` directly, and on some systems HDF5 is not
# automatically found (either through netCDF or by being in a default
# include directory)
function(link_hdf5_if_required TARGET)
  set(CMAKE_REQUIRED_LIBRARIES netCDF::netcdf)
  check_c_source_compiles("
#include <H5public.h>
int main() {}"
    _hdf5_links)

  if (NOT _hdf5_links)
    message(STATUS "Linking against HDF5 explicitly")
    find_package(HDF5 COMPONENTS C REQUIRED)
    target_link_libraries(${TARGET} PRIVATE hdf5::hdf5)
  endif()
endfunction()
