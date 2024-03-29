cmake_minimum_required(VERSION 3.15...3.23)

project(CmakeConfigPackageTests LANGUAGES CXX)

# ---- Test as standalone project the exported config package ----

if(PROJECT_IS_TOP_LEVEL OR TEST_INSTALLED_VERSION)
  enable_testing()

  find_package(nc_complex CONFIG REQUIRED) # for intro, project_options, ...

  if(NOT TARGET nc_complex_options)
    message(FATAL_ERROR "Required config package not found!")
    return() # be strictly paranoid for Template Janitor github action! CK
  endif()
endif()


include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)

add_library(test_utils
  STATIC
  test_utilities.cpp
  test_utilities.h
)
target_link_libraries(test_utils
  PRIVATE
  nc_complex_warnings
  nc_complex_options
  PUBLIC
  plasmafair::nc_complex
)

add_executable(test_c_api test_c_api.cpp)
target_link_libraries(test_c_api
  PRIVATE
  nc_complex_warnings
  nc_complex_options
  test_utils
  Catch2::Catch2WithMain
)
target_compile_features(test_c_api PUBLIC cxx_std_17)

if(WIN32 AND BUILD_SHARED_LIBS)
  add_custom_command(
    TARGET test_c_api
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:test_c_api> $<TARGET_FILE_DIR:test_c_api>
    COMMAND_EXPAND_LISTS)
endif()

catch_discover_tests(
  test_c_api
  TEST_PREFIX "unittests.c."
)

if (nc_complex_BUILD_CXX)
  add_executable(test_cpp_api test_cpp_api.cpp)
  target_link_libraries(test_cpp_api
    PRIVATE
    nc_complex_warnings
    nc_complex_options
    test_utils
    plasmafair::nc_complex_cpp
    Catch2::Catch2WithMain
  )
  target_compile_features(test_cpp_api PUBLIC cxx_std_17)

  catch_discover_tests(
    test_cpp_api
    TEST_PREFIX "unittests.cpp."
  )
endif()
