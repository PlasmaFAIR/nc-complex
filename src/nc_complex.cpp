#include "nc_complex/nc_complex.h"

#include <cstddef>
#include <netcdf.h>

#include <vector>
#include <complex>
#include <string>


namespace plasmafair {

constexpr auto complex_struct_name = "PF_NC_COMPLEX_TYPE";

bool file_has_complex_struct(int ncid) {
  nc_type typeidp{};
  const auto err = nc_inq_typeid(ncid, complex_struct_name, &typeidp);
  return (err == NC_NOERR) and (typeidp > 0);
}

bool compound_type_is_compatible(int nc_typeid) {
  // get fields of type

  // check names

  // check sizes

  return false;
}

bool variable_has_complex_dimension(int ncid, int nc_varid) {
  int num_dims{};
  nc_inq_varndims(ncid, nc_varid, &num_dims);

  std::vector<int> dim_ids(static_cast<std::size_t>(num_dims));
  nc_inq_vardimid(ncid, nc_varid, dim_ids.data());

  for (const auto& dim_id : dim_ids) {
    std::size_t length{};
    nc_inq_dimlen(ncid, dim_id, &length);

    if (length != 2) {
      continue;
    }

    std::string name;
    name.reserve(NC_MAX_NAME + 1);
    nc_inq_dimname(ncid, dim_id, name.data());
    name.shrink_to_fit();

    if (name == "ri") return true;
  }

  return false;
}

int nc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                               const size_t *countp,
                               const std::complex<double> *op);

int nc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                               const size_t *countp, std::complex<double> *ip) {
  return nc_get_var(ncid, varid, ip);
}

namespace details {
// Vector of ones for get/put_var1 functions
static constexpr size_t coord_one[NC_MAX_VAR_DIMS] = {1};
}

} // namespace plasmafair

int nc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                               const size_t *countp,
                               const double _Complex *op) {
  return plasmafair::nc_put_vara_double_complex(
      ncid, varid, startp, countp,
      reinterpret_cast<const std::complex<double> *>(op));
}

int nc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                               const size_t *countp, double _Complex *ip) {
  return plasmafair::nc_get_vara_double_complex(
      ncid, varid, startp, countp,
      reinterpret_cast<std::complex<double> *>(ip));
}

int nc_put_var1_double_complex(int ncid, int varid, const size_t *indexp,
                               const double_complex *data) {
  return nc_put_vara_double_complex(ncid, varid, indexp,
                                    plasmafair::details::coord_one, data);
}

int nc_get_var1_double_complex(int ncid, int varid, const size_t *indexp,
                               double_complex *data) {
  return nc_get_vara_double_complex(ncid, varid,
                                    plasmafair::details::coord_one,
                                    plasmafair::details::coord_one, data);
}
