cmake_minimum_required(VERSION 3.15...3.23)

project(CmakeConfigPackageTests LANGUAGES CXX)

# ---- Test as standalone project the exported config package ----

if(PROJECT_IS_TOP_LEVEL OR TEST_INSTALLED_VERSION)
  enable_testing()

  find_package(myproject CONFIG REQUIRED) # for intro, project_options, ...

  if(NOT TARGET myproject_options)
    message(FATAL_ERROR "Requiered config package not found!")
    return() # be strictly paranoid for Template Janitor github action! CK
  endif()
endif()

add_executable(basic_read_write basic.cpp)
target_link_libraries(basic_read_write plasmafair::nc_complex)
add_test(NAME basic_read_write COMMAND basic_read_write)

add_executable(create_with_api with_api.cpp)
target_link_libraries(create_with_api plasmafair::nc_complex)
add_test(NAME create_with_api COMMAND create_with_api)
