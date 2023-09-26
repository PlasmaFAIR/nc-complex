#include "nc_complex/nc_complex.h"

#include <array>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <netcdf.h>
#include <string>

#define CHECK(func)                                                                    \
  do {                                                                                 \
    if (const auto res = (func)) {                                                     \
      printf("Bailing out in file %s, line %d, error: (%d) %s.\n", __FILE__, __LINE__, \
             res, nc_strerror(res));                                                   \
      return false;                                                                    \
    }                                                                                  \
  } while (0)

constexpr int len_x = 3;
constexpr int len_ri = 2;

constexpr std::array<std::complex<double>, len_x> data = {
    {{0., 1.}, {2., 3.}, {4., 5}}};

constexpr auto filename = "with_api.nc";

/// Create a netCDF file with a variety of complex conventions
bool create_file(const std::string &filename) {
  int ncid = 0;
  int res = 0;
  CHECK(nc_create(filename.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

  int x_dim_id = 0;
  CHECK(nc_def_dim(ncid, "x", len_x, &x_dim_id));

  int type_id{};
  CHECK(pfnc_get_double_complex_typeid(ncid, &type_id));
  int var_struct_id = 0;
  const std::array<int, 1> dim_struct_ids{{x_dim_id}};
  CHECK(nc_def_var(ncid, "data_struct", type_id, 1, dim_struct_ids.data(),
                   &var_struct_id));

  CHECK(nc_put_var(ncid, var_struct_id, data.data()));

  CHECK(nc_close(ncid));

  return true;
}

/// Check that a complex array matches the expected result
bool check_data(std::array<std::complex<double>, len_x> const &data_in) {
  if (data_in == data) {
    printf("  Success!\n");
    return true;
  }

  printf("Failed:\n");
  for (const auto &elem : data_in) {
    printf("  %.1f %+.1f\n", elem.real(), elem.imag());
  }
  return false;
}

/// Read the test file using the nc_complex API
bool read_file_nc_complex(const std::string &filename) {
  int ncid = 0;
  int res = 0;
  CHECK(nc_open(filename.c_str(), NC_NOWRITE, &ncid));

  constexpr size_t zeros[NC_MAX_VAR_DIMS] = {0};

  std::array<std::complex<double>, len_x> data_struct_out;
  int var_struct_id = 0;
  CHECK(nc_inq_varid(ncid, "data_struct", &var_struct_id));
  CHECK(pfnc_get_vara_double_complex(ncid, var_struct_id, zeros, nullptr,
                                     cpp_to_c_complex(data_struct_out.data())));
  return check_data(data_struct_out);
}

int main() {
  printf("Creating file... ");
  if (not create_file(filename)) {
    return EXIT_FAILURE;
  }
  printf("done!\n");

  printf("pfnc_get_vara_double_complex, with type checking:\n");
  if (not read_file_nc_complex(filename)) {
    return EXIT_FAILURE;
  }
}
