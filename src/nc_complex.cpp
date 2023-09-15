#include "nc_complex/nc_complex.h"

#include <cstring>
#include <netcdf.h>

#include <cctype>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#define CHECK(func)                                                                    \
  do {                                                                                 \
    if (const auto res = (func)) {                                                     \
      return res;                                                                      \
    }                                                                                  \
  } while (0)

using namespace std::string_view_literals;

namespace plasmafair {

constexpr auto double_complex_struct_name = "_PFNC_DOUBLE_COMPLEX_TYPE";

/// Return true if file already has our complex type
bool file_has_complex_struct(int ncid, nc_type &typeidp) {
  const auto err = nc_inq_typeid(ncid, double_complex_struct_name, &typeidp);
  return (err == NC_NOERR) and (typeidp > 0);
}

/// Create complex datatype if it doesn't already exist
int create_double_complex_struct(int ncid, nc_type &type_id) {
  if (file_has_complex_struct(ncid, type_id)) {
    return NC_NOERR;
  }

  CHECK(nc_def_compound(ncid, sizeof(std::complex<double>), double_complex_struct_name,
                        &type_id));
  CHECK(nc_insert_compound(ncid, type_id, "r", 0, NC_DOUBLE));
  CHECK(nc_insert_compound(ncid, type_id, "i", sizeof(double), NC_DOUBLE));

  return NC_NOERR;
}

/// Return true if a compound type is compatible with a known convention
bool compound_type_is_compatible(int ncid, int nc_typeid) {

  // Does the name matching a known convention?
  std::string name;
  name.reserve(NC_MAX_NAME + 1);
  nc_inq_compound_name(ncid, nc_typeid, name.data());
  if (name == double_complex_struct_name) {
    return true;
  }

  // Does it have exactly two fields?
  std::size_t num_fields{};
  nc_inq_compound_nfields(ncid, nc_typeid, &num_fields);
  if (num_fields != 2) {
    return false;
  }

  // As far as I can tell, all conventions put the real part first and
  // the imaginary part second. I'm pretty sure all natiev language
  // types are also this way round. That means we don't have to worry
  // about trying both combinations!
  std::string real_name;
  real_name.reserve(NC_MAX_NAME + 1);
  std::size_t real_offset;
  nc_type real_field_type;
  int real_rank;
  nc_inq_compound_field(ncid, nc_typeid, 0, real_name.data(), &real_offset,
                        &real_field_type, &real_rank, nullptr);

  // If it's not a floating type, we're not interested
  if (!(real_field_type == NC_FLOAT || real_field_type == NC_DOUBLE)) {
    return false;
  }
  // Also needs to be scalar
  if (real_rank != 0) {
    return false;
  }

  // Now check names. For now, just check it starts with "r", in any case
  if (std::tolower(real_name.front()) != 'r') {
    return false;
  }

  std::string imag_name;
  imag_name.reserve(NC_MAX_NAME + 1);
  std::size_t imag_offset;
  nc_type imag_field_type;
  int imag_rank;
  nc_inq_compound_field(ncid, nc_typeid, 0, imag_name.data(), &imag_offset,
                        &imag_field_type, &imag_rank, nullptr);

  // Both component types better match
  if (imag_field_type != real_field_type) {
    return false;
  }
  if (imag_rank != 0) {
    return false;
  }
  if (std::tolower(imag_name.front()) != 'i') {
    return false;
  }

  return true;
}

/// Return true if a given dimension matches a known convention
bool dimension_is_complex(int ncid, int dim_id) {
  std::size_t length{};
  nc_inq_dimlen(ncid, dim_id, &length);

  // Definitely can only be exactly two. Note that we can't catch
  // unlimited dimensions that only have two records so far.
  if (length != 2) {
    return false;
  }

  // Not sure if this is the best way, but here we are.
  std::string name;
  name.resize(NC_MAX_NAME + 1);
  nc_inq_dimname(ncid, dim_id, name.data());
  name.resize(std::strlen(name.c_str()));

  // Check against known names of complex dimensions
  if (name == "ri"sv) {
    return true;
  }

  return false;
}

/// Return true if a variable uses the dimension-convention
bool variable_has_complex_dimension(int ncid, int nc_varid) {
  int num_dims{};
  nc_inq_varndims(ncid, nc_varid, &num_dims);

  std::vector<int> dim_ids(static_cast<std::size_t>(num_dims));
  nc_inq_vardimid(ncid, nc_varid, dim_ids.data());

  // Now we check if any of the dimensions match one of our known
  // conventions. Do we need to check all of them, or just the
  // first/last?
  for (const auto &dim_id : dim_ids) {
    if (dimension_is_complex(ncid, dim_id)) {
      return true;
    }
  }

  return false;
}

/// Return true if a netCDF datatype is a compound type
bool is_compound_type(int ncid, int type_id) {
  // There appears to be no API for detecting whether a type ID is a
  // primitive type, so we have to check ourselves
  switch (type_id) {
  case NC_NAT:
  case NC_BYTE:
  case NC_CHAR:
  case NC_SHORT:
  case NC_INT:
  case NC_FLOAT:
  case NC_DOUBLE:
  case NC_UBYTE:
  case NC_USHORT:
  case NC_UINT:
  case NC_INT64:
  case NC_UINT64:
  case NC_STRING:
    return false;
  }

  int class_type{};
  nc_inq_user_type(ncid, type_id, nullptr, nullptr, nullptr, nullptr, &class_type);
  return class_type == NC_COMPOUND;
}

/// Return true if the variable matches a known complex convention
bool check_variable_is_double_complex(int ncid, int varid) {
  nc_type var_type_id{};
  if (nc_inq_vartype(ncid, varid, &var_type_id)) {
    return false;
  }

  if (is_compound_type(ncid, var_type_id)) {
    return compound_type_is_compatible(ncid, var_type_id);
  }

  return variable_has_complex_dimension(ncid, varid);
}

int nc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                               const size_t *countp, const std::complex<double> *op) {
  if (!check_variable_is_double_complex(ncid, varid)) {
    return NC_EBADTYPE;
  }

  // TODO: handle start/count/stride correctly for dimension convention
  // TODO: handle converting different float sizes

  return nc_put_vara(ncid, varid, startp, countp, op);
}

int nc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                               const size_t *countp, std::complex<double> *ip) {
  if (!check_variable_is_double_complex(ncid, varid)) {
    return NC_EBADTYPE;
  }

  // TODO: handle start/count/stride correctly for dimension convention
  // TODO: handle converting different float sizes

  return nc_get_vara(ncid, varid, startp, countp, ip);
}

namespace details {
// Vector of ones for get/put_var1 functions
static constexpr size_t coord_one[NC_MAX_VAR_DIMS] = {1};
} // namespace details

} // namespace plasmafair

int pfnc_get_double_complex_typeid(int ncid, nc_type *complex_typeid) {
  return plasmafair::create_double_complex_struct(ncid, *complex_typeid);
}

int pfnc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, const double _Complex *op) {
  return plasmafair::nc_put_vara_double_complex(
      ncid, varid, startp, countp, reinterpret_cast<const std::complex<double> *>(op));
}

int pfnc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, double_complex *ip) {
  return plasmafair::nc_get_vara_double_complex(
      ncid, varid, startp, countp, reinterpret_cast<std::complex<double> *>(ip));
}

int pfnc_put_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 const double_complex *data) {
  return pfnc_put_vara_double_complex(ncid, varid, indexp,
                                      plasmafair::details::coord_one, data);
}

int pfnc_get_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 double_complex *data) {
  return pfnc_get_vara_double_complex(ncid, varid, plasmafair::details::coord_one,
                                      plasmafair::details::coord_one, data);
}
