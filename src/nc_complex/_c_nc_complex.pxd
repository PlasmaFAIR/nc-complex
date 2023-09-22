# Some declarations from netCDF that we need directly and can't use
# via the netCDF4 Python API
cdef extern from "netcdf.h":
  cdef enum:
    NC_NOERR

  const char* nc_strerror(int ncerr)


# Our API declarations
cdef extern from "nc_complex/nc_complex.h":
  int pfnc_get_double_complex_typeid(int ncid, int *nc_typeid)

  int pfnc_put_vara_double_complex(int ncid, int varid, const size_t *startp,
                                   const size_t *countp, const double complex *op)
  int pfnc_get_vara_double_complex(int ncid, int varid, const size_t *startp,
                                   const size_t *countp, double complex *ip)

  int pfnc_put_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                   const double complex *data)
  int pfnc_get_var1_double_complex(int ncid, int varid, const size_t *indexp,
                                   double complex *data)
