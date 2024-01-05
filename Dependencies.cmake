include(cmake/CPM.cmake)

# Done as a function so that updates to variables like
# CMAKE_CXX_FLAGS don't propagate out to other
# targets
function(nc_complex_setup_dependencies)
  include(CTest)
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

  if(BUILD_TESTING AND NOT TARGET Catch2::Catch2WithMain)
    cpmaddpackage("gh:catchorg/Catch2@3.3.2")
  endif()

  find_package(netCDF REQUIRED)

  if (nc_complex_BUILD_CXX)
    if (nc_complex_DOWNLOAD_NETCDF_CXX)
      set(NCXX_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
      cpmaddpackage("gh:ZedThree/netcdf-cxx4#722cbd5")
    else()
      find_package(netCDFCxx REQUIRED)
    endif()
  endif()
endfunction()
