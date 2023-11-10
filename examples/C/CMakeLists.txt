add_executable(example_without_nc_complex example_without_nc_complex.c)
target_link_libraries(example_without_nc_complex
  plasmafair::nc_complex
  nc_complex_warnings
  nc_complex_options
)
target_compile_features(example_without_nc_complex PUBLIC c_std_11)

add_executable(example_with_nc_complex example_with_nc_complex.c)
target_link_libraries(example_with_nc_complex
  plasmafair::nc_complex
  nc_complex_warnings
  nc_complex_options
)
target_compile_features(example_with_nc_complex PUBLIC c_std_11)
