#include "nc_complex/nc_complex.h"

#include <netcdf.h>

#include <ctype.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(func)                                                                    \
  do {                                                                                 \
    int res;                                                                           \
    if ((res = (func))) {                                                              \
      return res;                                                                      \
    }                                                                                  \
  } while (0)

// Vector of ones for get/put_var1 functions
static const size_t coord_one[NC_MAX_VAR_DIMS] = {1};

static const char *double_complex_struct_name = "_PFNC_DOUBLE_COMPLEX_TYPE";
static const char *float_complex_struct_name = "_PFNC_FLOAT_COMPLEX_TYPE";

/// Return true if file already has our complex type
bool file_has_double_complex_struct(int ncid, nc_type *typeidp) {
  const int err = nc_inq_typeid(ncid, double_complex_struct_name, typeidp);
  return (err == NC_NOERR) && (*typeidp > 0);
}

bool file_has_float_complex_struct(int ncid, nc_type *typeidp) {
  const int err = nc_inq_typeid(ncid, float_complex_struct_name, typeidp);
  return (err == NC_NOERR) && (*typeidp > 0);
}

bool pfnc_is_complex(int ncid, int varid) {
  return pfnc_is_complex_type(ncid, varid) || pfnc_has_complex_dimension(ncid, varid);
}

int pfnc_complex_base_type(int ncid, int nc_typeid, int *base_type_id) {
  if (nc_typeid < NC_MAX_ATOMIC_TYPE) {
    *base_type_id = nc_typeid;
    return NC_NOERR;
  }

  // TODO: This should probably handle vlens too

  return nc_inq_compound_field(ncid, nc_typeid, 0, NULL, NULL, base_type_id, NULL,
                               NULL);
}

int pfnc_inq_var_complex_base_type(int ncid, int varid, int *nc_typeid) {
  nc_type var_type_id;
  int ierr = nc_inq_vartype(ncid, varid, &var_type_id);
  if (ierr != NC_NOERR) {
    return ierr;
  }
  return pfnc_complex_base_type(ncid, var_type_id, nc_typeid);
}

/// Return true if a compound type is compatible with a known convention
bool compound_type_is_compatible(int ncid, int nc_typeid) {

  // Does the name matching a known convention?
  char name[NC_MAX_NAME + 1];
  nc_inq_compound_name(ncid, nc_typeid, name);
  if (name == double_complex_struct_name) {
    return true;
  }

  // Does it have exactly two fields?
  size_t num_fields;
  nc_inq_compound_nfields(ncid, nc_typeid, &num_fields);
  if (num_fields != 2) {
    return false;
  }

  // As far as I can tell, all conventions put the real part first and
  // the imaginary part second. I'm pretty sure all native language
  // types are also this way round. That means we don't have to worry
  // about trying both combinations!
  char real_name[NC_MAX_NAME + 1];
  size_t real_offset;
  nc_type real_field_type;
  int real_rank;
  nc_inq_compound_field(ncid, nc_typeid, 0, real_name, &real_offset, &real_field_type,
                        &real_rank, NULL);

  // If it's not a floating type, we're not interested
  if (!(real_field_type == NC_FLOAT || real_field_type == NC_DOUBLE)) {
    return false;
  }
  // Also needs to be scalar
  if (real_rank != 0) {
    return false;
  }

  // Now check names. For now, just check it starts with "r", in any case
  if (tolower(real_name[0]) != 'r') {
    return false;
  }

  char imag_name[NC_MAX_NAME + 1];
  size_t imag_offset;
  nc_type imag_field_type;
  int imag_rank;
  nc_inq_compound_field(ncid, nc_typeid, 1, imag_name, &imag_offset, &imag_field_type,
                        &imag_rank, NULL);

  // Both component types better match
  if (imag_field_type != real_field_type) {
    return false;
  }
  if (imag_rank != 0) {
    return false;
  }
  if (tolower(imag_name[0]) != 'i') {
    return false;
  }

  return true;
}

/// Return true if a given dimension matches a known convention
bool dimension_is_complex(int ncid, int dim_id) {
  size_t length;
  nc_inq_dimlen(ncid, dim_id, &length);

  // Definitely can only be exactly two. Note that we can't catch
  // unlimited dimensions that only have two records so far.
  if (length != 2) {
    return false;
  }

  // Not sure if this is the best way, but here we are.
  char name[NC_MAX_NAME + 1];
  nc_inq_dimname(ncid, dim_id, name);

  const size_t name_length = strlen(name);

  // Check against known names of complex dimensions
  // TODO: Generalise
  if (strncmp(name, "ri", name_length) == 0) {
    return true;
  }

  if (strncmp(name, "complex", name_length) == 0) {
    return true;
  }

  return false;
}

/// Return true if a variable uses the dimension-convention
bool pfnc_has_complex_dimension(int ncid, int nc_varid) {
  int num_dims;
  nc_inq_varndims(ncid, nc_varid, &num_dims);

  int *dim_ids = (int *)malloc((size_t)num_dims * sizeof(int));
  nc_inq_vardimid(ncid, nc_varid, dim_ids);

  // Now we check if any of the dimensions match one of our known
  // conventions. Do we need to check all of them, or just the
  // first/last?
  for (int i = 0; i < num_dims; i++) {
    if (dimension_is_complex(ncid, dim_ids[i])) {
      free(dim_ids);
      return true;
    }
  }

  free(dim_ids);
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

  int class_type;
  nc_inq_user_type(ncid, type_id, NULL, NULL, NULL, NULL, &class_type);
  return class_type == NC_COMPOUND;
}

/// Copy an array meant for a complex-dimensioned variable
size_t *copy_complex_dim_size_t_array(const size_t *old_array, int numdims,
                                      size_t complex_dim_value) {
  size_t *new_buffer = NULL;

  if (old_array != NULL) {
    new_buffer = (size_t *)malloc(sizeof(size_t) * (size_t)numdims);

    size_t last_dim = (size_t)(numdims - 1);
    for (size_t i = 0; i < last_dim; i++) {
      new_buffer[i] = old_array[i];
    }

    new_buffer[last_dim] = complex_dim_value;
  }
  return new_buffer;
}

ptrdiff_t *copy_complex_dim_ptrdiff_t_array(const ptrdiff_t *old_array, int numdims,
                                            ptrdiff_t complex_dim_value) {
  ptrdiff_t *new_buffer = NULL;

  if (old_array != NULL) {
    new_buffer = (ptrdiff_t *)malloc(sizeof(ptrdiff_t) * (size_t)numdims);

    size_t last_dim = (size_t)(numdims - 1);
    for (size_t i = 0; i < last_dim; i++) {
      new_buffer[i] = old_array[i];
    }

    new_buffer[last_dim] = complex_dim_value;
  }
  return new_buffer;
}

bool pfnc_is_complex_type(int ncid, int varid) {
  nc_type var_type_id;
  if (nc_inq_vartype(ncid, varid, &var_type_id)) {
    return false;
  }

  if (is_compound_type(ncid, var_type_id)) {
    return compound_type_is_compatible(ncid, var_type_id);
  }
  return false;
}

int pfnc_get_double_complex_typeid(int ncid, nc_type *type_id) {
  // TODO: Error if not netCDF4

  if (file_has_double_complex_struct(ncid, type_id)) {
    return NC_NOERR;
  }

  CHECK(nc_def_compound(ncid, sizeof(double_complex), double_complex_struct_name,
                        type_id));
  CHECK(nc_insert_compound(ncid, *type_id, "r", 0, NC_DOUBLE));
  CHECK(nc_insert_compound(ncid, *type_id, "i", sizeof(double), NC_DOUBLE));

  return NC_NOERR;
}

int pfnc_get_float_complex_typeid(int ncid, nc_type *type_id) {
  // TODO: Error if not netCDF4

  if (file_has_float_complex_struct(ncid, type_id)) {
    return NC_NOERR;
  }

  CHECK(
      nc_def_compound(ncid, sizeof(float_complex), float_complex_struct_name, type_id));
  CHECK(nc_insert_compound(ncid, *type_id, "r", 0, NC_FLOAT));
  CHECK(nc_insert_compound(ncid, *type_id, "i", sizeof(float), NC_FLOAT));

  return NC_NOERR;
}

int pfnc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, const double _Complex *op) {
  if (!pfnc_is_complex(ncid, varid)) {
    return NC_EBADTYPE;
  }

  // TODO: handle start/count/stride correctly for dimension convention
  // TODO: handle converting different float sizes

  return nc_put_vara(ncid, varid, startp, countp, op);
}

int pfnc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, double_complex *ip) {
  if (!pfnc_is_complex(ncid, varid)) {
    return NC_EBADTYPE;
  }

  // TODO: handle converting different float sizes

  // Check if we can get away without fudging count/start sizes
  if (((startp == NULL) && (countp == NULL)) ||
      !pfnc_has_complex_dimension(ncid, varid)) {
    return nc_get_vara(ncid, varid, startp, countp, ip);
  }

  // The real variable has a complex dimension, but we're pretending
  // it doesn't, so now we need start/count arrays of the real size

  int numdims = 0;
  {
    const int ierr = nc_inq_varndims(ncid, varid, &numdims);
    if (ierr != NC_NOERR) {
      return ierr;
    }
  }

  // Copy start/count buffers, appending an extra element for the
  // complex dimension. This dimension starts at 0 and has 2 elements
  size_t *start_buffer = copy_complex_dim_size_t_array(startp, numdims, 0);
  size_t *count_buffer = copy_complex_dim_size_t_array(countp, numdims, 2);

  const int ierr = nc_get_vara(ncid, varid, start_buffer, count_buffer, ip);

  if (start_buffer != NULL) {
    free(start_buffer);
  }
  if (count_buffer != NULL) {
    free(count_buffer);
  }
  return ierr;
}

int pfnc_put_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 const double_complex *data) {
  return pfnc_put_vara_double_complex(ncid, varid, indexp, coord_one, data);
}

int pfnc_get_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 double_complex *data) {
  return pfnc_get_vara_double_complex(ncid, varid, indexp, coord_one, data);
}

int pfnc_put_vara_float_complex(int ncid, int varid, const size_t *startp,
                                const size_t *countp, const float _Complex *op) {
  if (!pfnc_is_complex(ncid, varid)) {
    return NC_EBADTYPE;
  }

  // TODO: handle start/count/stride correctly for dimension convention
  // TODO: handle converting different float sizes

  return nc_put_vara(ncid, varid, startp, countp, op);
}

int pfnc_get_vara_float_complex(int ncid, int varid, const size_t *startp,
                                const size_t *countp, float_complex *ip) {
  if (!pfnc_is_complex(ncid, varid)) {
    return NC_EBADTYPE;
  }

  // TODO: handle converting different float sizes

  // Check if we can get away without fudging count/start sizes
  if (((startp == NULL) && (countp == NULL)) ||
      !pfnc_has_complex_dimension(ncid, varid)) {
    return nc_get_vara(ncid, varid, startp, countp, ip);
  }

  // The real variable has a complex dimension, but we're pretending
  // it doesn't, so now we need start/count arrays of the real size

  int numdims = 0;
  {
    const int ierr = nc_inq_varndims(ncid, varid, &numdims);
    if (ierr != NC_NOERR) {
      return ierr;
    }
  }

  // Copy start/count buffers, appending an extra element for the
  // complex dimension. This dimension starts at 0 and has 2 elements
  size_t *start_buffer = copy_complex_dim_size_t_array(startp, numdims, 0);
  size_t *count_buffer = copy_complex_dim_size_t_array(countp, numdims, 2);

  const int ierr = nc_get_vara(ncid, varid, start_buffer, count_buffer, ip);

  if (start_buffer != NULL) {
    free(start_buffer);
  }
  if (count_buffer != NULL) {
    free(count_buffer);
  }
  return ierr;
}

int pfnc_put_var1_float_complex(int ncid, int varid, const size_t *indexp,
                                const float_complex *data) {
  return pfnc_put_vara_float_complex(ncid, varid, indexp, coord_one, data);
}

int pfnc_get_var1_float_complex(int ncid, int varid, const size_t *indexp,
                                float_complex *data) {
  return pfnc_get_vara_float_complex(ncid, varid, indexp, coord_one, data);
}

int pfnc_inq_varndims(int ncid, int varid, int *ndimsp) {
  const int ierr = nc_inq_varndims(ncid, varid, ndimsp);
  if (ierr != NC_NOERR) {
    return ierr;
  }
  // Pretend that variable has one less dimension than it does
  if (pfnc_has_complex_dimension(ncid, varid)) {
    *ndimsp -= 1;
  }
  return NC_NOERR;
}

int pfnc_inq_vardimid(int ncid, int varid, int *dimidsp) {
  // Tricky bit: if variable has complex dimension, and user used
  // pfnc_inq_varndims, then dimidsp is one smaller than netCDF thinks
  // it should be. So we'll have to allocate our own array of the
  // correct size and copy out of that.

  // This buffer will point to either the user's array, or our own one
  int *buffer = dimidsp;
  bool allocated_internal_buffer = false;
  int numdims = 0;

  if (pfnc_has_complex_dimension(ncid, varid)) {
    const int ierr = nc_inq_varndims(ncid, varid, &numdims);
    if (ierr != NC_NOERR) {
      return ierr;
    }
    buffer = (int *)malloc(sizeof(int) * (size_t)numdims);
    allocated_internal_buffer = true;
  }

  int ierr = nc_inq_vardimid(ncid, varid, buffer);
  if (ierr != NC_NOERR) {
    goto cleanup;
  }

  if (allocated_internal_buffer) {
    if (numdims <= 0) {
      // This should never happen
      goto cleanup;
    }
    const size_t other_dims = (size_t)(numdims - 1);
    for (size_t i = 0; i < other_dims; i++) {
      dimidsp[i] = buffer[i];
    }
  }

cleanup:
  if (allocated_internal_buffer) {
    free(buffer);
  }
  return ierr;
}

int pfnc_get_vara(int ncid, int varid, const size_t *startp, const size_t *countp,
                  void *ip) {
  if (pfnc_is_complex(ncid, varid)) {
    nc_type base_type;
    const int ierr = pfnc_inq_var_complex_base_type(ncid, varid, &base_type);
    if (ierr != NC_NOERR) {
      return ierr;
    }
    switch (base_type) {
    case NC_DOUBLE:
      return pfnc_get_vara_double_complex(ncid, varid, startp, countp, ip);
    case NC_FLOAT:
      return pfnc_get_vara_float_complex(ncid, varid, startp, countp, ip);
    default:
      return NC_EBADTYPE;
    }
  }

  return nc_get_vara(ncid, varid, startp, countp, ip);
}
