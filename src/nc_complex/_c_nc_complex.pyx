"""
Cython helper functions for the C nc_complex library functions
"""



def _check_netcdf_error(extra_msg: str, ierr: int):
    # print netcdf error message, raise error.
    if ierr == NC_NOERR:
        return

    err_str = (<const char *>nc_strerror(ierr)).decode('ascii')
    raise RuntimeError(f"{extra_msg}: {err_str}")


cpdef int double_complex_typeid(int ncid):
    cdef int nc_typeid
    ierr = pfnc_get_double_complex_typeid(ncid, &nc_typeid)
    _check_netcdf_error("When attempting to get typeid for double complex", ierr)
    return nc_typeid


cpdef bint is_complex(int ncid, int varid):
    return <bint>pfnc_is_complex(ncid, varid)

cpdef bint is_complex_type(int ncid, int varid):
    return <bint>pfnc_is_complex_type(ncid, varid)

cpdef bint has_complex_dimension(int ncid, int varid):
    return <bint>pfnc_has_complex_dimension(ncid, varid)

