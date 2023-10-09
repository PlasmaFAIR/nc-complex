include(cmake/CPM.cmake)

# Done as a function so that updates to variables like
# CMAKE_CXX_FLAGS don't propagate out to other
# targets
function(nc_complex_setup_dependencies)
  if(NOT TARGET Catch2::Catch2WithMain AND NOT SKBUILD)
    cpmaddpackage("gh:catchorg/Catch2@3.3.2")
  endif()
endfunction()
