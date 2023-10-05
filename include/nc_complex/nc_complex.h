#ifndef PLASMA_FAIR_NC_COMPLEX
#define PLASMA_FAIR_NC_COMPLEX

#include <netcdf.h>

#ifdef _MSC_VER
#include <complex.h>
typedef _Dcomplex double_complex;
#else
#if defined(__cplusplus) && defined(__clang__)
#include <complex>
using double_complex = std::complex<double>;
#else
typedef double _Complex double_complex;
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

extern "C" {
#endif

/// Return true if variable is complex
bool pfnc_is_complex(int ncid, int varid);
/// Return true if variable is complex and uses a compound datatype
bool pfnc_is_complex_type(int ncid, int varid);
/// Return true if variable is complex and has a complex dimension
/// (assumed to be the last dimension)
bool pfnc_has_complex_dimension(int ncid, int varid);

/// Create complex datatype if it doesn't already exist
int pfnc_get_double_complex_typeid(int ncid, nc_type *nc_typeid);

int pfnc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, const double_complex *op);

int pfnc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, double_complex *ip);

int pfnc_put_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 const double_complex *data);
int pfnc_get_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 double_complex *data);


// Custom shims for lying about dimensional variables
int pfnc_inq_varndims(int ncid, int varid, int *ndimsp);
int pfnc_inq_vardimid(int ncid, int varid, int *dimidsp);

int pfnc_get_vara(int ncid, int varid, const size_t *startp, const size_t *countp,
                  void *ip);
// TODO: pfnc_libvers

#ifdef __cplusplus
}
#endif

#endif
