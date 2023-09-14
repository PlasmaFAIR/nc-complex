#ifndef PLASMA_FAIR_NC_COMPLEX
#define PLASMA_FAIR_NC_COMPLEX

#ifdef _MSC_VER
#include <complex.h>
typedef _Dcomplex double_complex;
#else
typedef double _Complex double_complex;
#endif

#include <complex.h>
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

int pfnc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, const double_complex *op);

int pfnc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                                 const size_t *countp, double_complex *ip);

int pfnc_put_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 const double_complex *data);
int pfnc_get_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                 double_complex *data);

#ifdef __cplusplus
}
#endif

#endif
