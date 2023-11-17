include(cmake/SystemLink.cmake)
include(cmake/LibFuzzer.cmake)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)


macro(nc_complex_supports_sanitizers)
  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND NOT WIN32)
    set(SUPPORTS_UBSAN ON)
  else()
    set(SUPPORTS_UBSAN OFF)
  endif()

  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND WIN32)
    set(SUPPORTS_ASAN OFF)
  else()
    set(SUPPORTS_ASAN ON)
  endif()
endmacro()

macro(nc_complex_setup_options)
  option(BUILD_SHARED_LIBS "Build shared libs" ON)
  option(nc_complex_ENABLE_HARDENING "Enable hardening" ON)
  option(nc_complex_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  cmake_dependent_option(
    nc_complex_ENABLE_GLOBAL_HARDENING
    "Attempt to push hardening options to built dependencies"
    ON
    nc_complex_ENABLE_HARDENING
    OFF)

  nc_complex_supports_sanitizers()

  if(NOT PROJECT_IS_TOP_LEVEL OR nc_complex_PACKAGING_MAINTAINER_MODE)
    option(nc_complex_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(nc_complex_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(nc_complex_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(nc_complex_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
    option(nc_complex_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(nc_complex_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" OFF)
    option(nc_complex_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(nc_complex_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(nc_complex_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(nc_complex_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
    option(nc_complex_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
    option(nc_complex_ENABLE_PCH "Enable precompiled headers" OFF)
  else()
    option(nc_complex_ENABLE_IPO "Enable IPO/LTO" ON)
    option(nc_complex_WARNINGS_AS_ERRORS "Treat Warnings As Errors" ON)
    option(nc_complex_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(nc_complex_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" ${SUPPORTS_ASAN})
    option(nc_complex_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(nc_complex_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" ${SUPPORTS_UBSAN})
    option(nc_complex_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(nc_complex_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(nc_complex_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(nc_complex_ENABLE_CLANG_TIDY "Enable clang-tidy" ON)
    option(nc_complex_ENABLE_CPPCHECK "Enable cpp-check analysis" ON)
    option(nc_complex_ENABLE_PCH "Enable precompiled headers" OFF)
    option(nc_complex_BUILD_EXAMPLES "Build nc-compelx examples" ON)
  endif()

  if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
      nc_complex_ENABLE_IPO
      nc_complex_WARNINGS_AS_ERRORS
      nc_complex_ENABLE_USER_LINKER
      nc_complex_ENABLE_SANITIZER_ADDRESS
      nc_complex_ENABLE_SANITIZER_LEAK
      nc_complex_ENABLE_SANITIZER_UNDEFINED
      nc_complex_ENABLE_SANITIZER_THREAD
      nc_complex_ENABLE_SANITIZER_MEMORY
      nc_complex_ENABLE_UNITY_BUILD
      nc_complex_ENABLE_CLANG_TIDY
      nc_complex_ENABLE_CPPCHECK
      nc_complex_ENABLE_COVERAGE
      nc_complex_ENABLE_PCH
      nc_complex_ENABLE_CACHE
      nc_complex_BUILD_EXAMPLES
    )
  endif()

  option(nc_complex_ENABLE_Fortran "Build Fortran API" OFF)
  option(nc_complex_ENABLE_CXX "Build C++ API" OFF)

  if (nc_complex_ENABLE_CXX)
    nc_complex_check_libfuzzer_support(LIBFUZZER_SUPPORTED)
    if(LIBFUZZER_SUPPORTED AND (nc_complex_ENABLE_SANITIZER_ADDRESS OR nc_complex_ENABLE_SANITIZER_THREAD OR nc_complex_ENABLE_SANITIZER_UNDEFINED))
      set(DEFAULT_FUZZER ON)
    else()
      set(DEFAULT_FUZZER OFF)
    endif()

    option(nc_complex_BUILD_FUZZ_TESTS "Enable fuzz testing executable" ${DEFAULT_FUZZER})
  endif()

endmacro()

macro(nc_complex_global_options)
  if(nc_complex_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    nc_complex_enable_ipo()
  endif()

  nc_complex_supports_sanitizers()

  if(nc_complex_ENABLE_HARDENING AND nc_complex_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR nc_complex_ENABLE_SANITIZER_UNDEFINED
       OR nc_complex_ENABLE_SANITIZER_ADDRESS
       OR nc_complex_ENABLE_SANITIZER_THREAD
       OR nc_complex_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    message("${nc_complex_ENABLE_HARDENING} ${ENABLE_UBSAN_MINIMAL_RUNTIME} ${nc_complex_ENABLE_SANITIZER_UNDEFINED}")
    nc_complex_enable_hardening(nc_complex_options ON ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()

macro(nc_complex_local_options)
  if(PROJECT_IS_TOP_LEVEL)
    include(cmake/StandardProjectSettings.cmake)
  endif()

  add_library(nc_complex_warnings INTERFACE)
  add_library(nc_complex_options INTERFACE)

  include(cmake/CompilerWarnings.cmake)
  nc_complex_set_project_warnings(
    nc_complex_warnings
    ${nc_complex_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
    "")

  if(nc_complex_ENABLE_USER_LINKER)
    include(cmake/Linker.cmake)
    configure_linker(nc_complex_options)
  endif()

  include(cmake/Sanitizers.cmake)
  nc_complex_enable_sanitizers(
    nc_complex_options
    ${nc_complex_ENABLE_SANITIZER_ADDRESS}
    ${nc_complex_ENABLE_SANITIZER_LEAK}
    ${nc_complex_ENABLE_SANITIZER_UNDEFINED}
    ${nc_complex_ENABLE_SANITIZER_THREAD}
    ${nc_complex_ENABLE_SANITIZER_MEMORY})

  set_target_properties(nc_complex_options PROPERTIES UNITY_BUILD ${nc_complex_ENABLE_UNITY_BUILD})

  if(nc_complex_ENABLE_PCH)
    target_precompile_headers(
      nc_complex_options
      INTERFACE
      <vector>
      <string>
      <utility>)
  endif()

  include(cmake/StaticAnalyzers.cmake)
  if(nc_complex_ENABLE_CLANG_TIDY)
    nc_complex_enable_clang_tidy(nc_complex_options ${nc_complex_WARNINGS_AS_ERRORS})
  endif()

  if(nc_complex_ENABLE_CPPCHECK)
    nc_complex_enable_cppcheck(${nc_complex_WARNINGS_AS_ERRORS} "" # override cppcheck options
    )
  endif()

  if(nc_complex_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    nc_complex_enable_coverage(nc_complex_options)
  endif()

  if(nc_complex_WARNINGS_AS_ERRORS)
    check_cxx_compiler_flag("-Wl,--fatal-warnings" LINKER_FATAL_WARNINGS)
    if(LINKER_FATAL_WARNINGS)
      # This is not working consistently, so disabling for now
      # target_link_options(nc_complex_options INTERFACE -Wl,--fatal-warnings)
    endif()
  endif()

  if(nc_complex_ENABLE_HARDENING AND NOT nc_complex_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR nc_complex_ENABLE_SANITIZER_UNDEFINED
       OR nc_complex_ENABLE_SANITIZER_ADDRESS
       OR nc_complex_ENABLE_SANITIZER_THREAD
       OR nc_complex_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    nc_complex_enable_hardening(nc_complex_options OFF ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()

endmacro()
