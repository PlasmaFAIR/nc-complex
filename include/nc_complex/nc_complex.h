#ifndef PLASMA_FAIR_NC_COMPLEX
#define PLASMA_FAIR_NC_COMPLEX

#include <netcdf.h>

#ifdef _MSC_VER
#include <complex.h>
typedef _Dcomplex double_complex;
typedef _Fcomplex float_complex;
#else
#if defined(__cplusplus) && defined(__clang__)
#include <complex>
using double_complex = std::complex<double>;
using float_complex = std::complex<float>;
#else
typedef double _Complex double_complex;
typedef float _Complex float_complex;
#endif
#endif

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
#include <complex>

inline double_complex *cpp_to_c_complex(std::complex<double> *data) {
  return reinterpret_cast<double_complex *>(data);
}

inline std::complex<double> *c_to_cpp_complex(double_complex *data) {
  return reinterpret_cast<std::complex<double> *>(data);
}

inline float_complex *cpp_to_c_complex(std::complex<float> *data) {
  return reinterpret_cast<float_complex *>(data);
}

inline std::complex<float> *c_to_cpp_complex(float_complex *data) {
  return reinterpret_cast<std::complex<float> *>(data);
}

extern "C" {
#endif

/// Datatype for float complex, for use with `pfnc_def_var`
///
/// Uses complex compound datatype with netCDF4 format, and complex dimension otherwise
#define PFNC_FLOAT_COMPLEX (NC_FIRSTUSERTYPEID - 4)
/// Datatype for float complex, for use with `pfnc_def_var`
///
/// Always use a complex dimension, regardless of file format
#define PFNC_FLOAT_COMPLEX_DIM (NC_FIRSTUSERTYPEID - 3)
/// Datatype for double complex, for use with `pfnc_def_var`
///
/// Uses complex compound datatype with netCDF4 format, and complex dimension otherwise
#define PFNC_DOUBLE_COMPLEX (NC_FIRSTUSERTYPEID - 2)
/// Datatype for double complex, for use with `pfnc_def_var`
///
/// Always use a complex dimension, regardless of file format
#define PFNC_DOUBLE_COMPLEX_DIM (NC_FIRSTUSERTYPEID - 1)

/// Return true if variable is complex
bool pfnc_var_is_complex(int ncid, int varid);
/// Return true if variable is complex and uses a compound datatype
bool pfnc_var_is_complex_type(int ncid, int varid);
/// Return true if variable is complex and has a complex dimension
/// (assumed to be the last dimension)
bool pfnc_var_has_complex_dimension(int ncid, int varid);

/// Return true if dimension is complex
bool pfnc_is_complex_dim(int ncid, int dim_id);

/// Create complex datatype if it doesn't already exist
int pfnc_get_double_complex_typeid(int ncid, nc_type *nc_typeid);
int pfnc_get_float_complex_typeid(int ncid, nc_type *nc_typeid);

/// Get complex dimension, creating one if it doesn't already exist
int pfnc_get_complex_dim(int ncid, int *nc_dim);

/// Get the base numerical type of a complex type
///
/// Returns the type of the components for a compound type, or the
/// type of an element for a dimension type.
int pfnc_complex_base_type(int ncid, int nc_typeid, int *base_type_id);

/// Get the base numerical type of a complex variable
int pfnc_inq_var_complex_base_type(int ncid, int varid, int *nc_typeid);

int pfnc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, const double_complex *op);

int pfnc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, double_complex *ip);

int pfnc_put_vars_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, const ptrdiff_t *stridep,
                                 const double_complex *op);

int pfnc_get_vars_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, const ptrdiff_t *stridep,
                                 double_complex *ip);

int pfnc_put_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 const double_complex *data);
int pfnc_get_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 double_complex *data);

int pfnc_put_vara_float_complex(int ncid, int varid, const size_t *startp,
                                const size_t *countp, const float_complex *op);

int pfnc_get_vara_float_complex(int ncid, int varid, const size_t *startp,
                                const size_t *countp, float_complex *ip);

int pfnc_put_vars_float_complex(int ncid, int varid, const size_t *startp,
                                const size_t *countp, const ptrdiff_t *stridep,
                                const float_complex *op);

int pfnc_get_vars_float_complex(int ncid, int varid, const size_t *startp,
                                const size_t *countp, const ptrdiff_t *stridep,
                                float_complex *ip);

int pfnc_put_var1_float_complex(int ncid, int varid, const size_t *indexp,
                                const float_complex *data);
int pfnc_get_var1_float_complex(int ncid, int varid, const size_t *indexp,
                                float_complex *data);

// Custom shims for lying about dimensional variables

/// Extension to `nc_def_var` that also accepts `PFNC_FLOAT_COMPLEX`,
/// `PFNC_FLOAT_COMPLEX_DIM`, `PFNC_DOUBLE_COMPLEX`, and `PFNC_DOUBLE_COMPLEX_DIM`
int pfnc_def_var(int ncid, const char *name, nc_type xtype, int ndims,
                 const int *dimidsp, int *varidp);

int pfnc_inq_var(int ncid, int varid, char *name, nc_type *xtypep, int *ndimsp,
                 int *dimidsp, int *nattsp);
int pfnc_inq_varndims(int ncid, int varid, int *ndimsp) {
  return pfnc_inq_var(ncid, varid, NULL, NULL, ndimsp, NULL, NULL);
}
int pfnc_inq_vardimid(int ncid, int varid, int *dimidsp) {
  return pfnc_inq_var(ncid, varid, NULL, NULL, NULL, dimidsp, NULL);
}

int pfnc_def_var_chunking(int ncid, int varid, int storage, const size_t *chunksizesp);
int pfnc_inq_var_chunking(int ncid, int varid, int *storagep, size_t *chunksizesp);

int pfnc_get_vara(int ncid, int varid, const size_t *startp, const size_t *countp,
                  void *ip);
int pfnc_get_vars(int ncid, int varid, const size_t *startp, const size_t *countp,
                  const ptrdiff_t *stridep, void *ip);

// TODO: pfnc_libvers

#ifdef __cplusplus
}
#endif

#endif
