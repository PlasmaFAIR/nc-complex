#include <netcdf.h>
#include <stddef.h>
#include <stdio.h>

#include "nc_complex/nc_complex.h"

#define len_x 3
#define len_complex 2

void print_complex(const char* note, double_complex z) {
    printf("%s %f%+f*i\n", note, creal(z), cimag(z));
}

int main(void) {
    const char* filename = "nc_complex_simple_example_without_nc_complex.nc";

    int ncid = 0;
    nc_create(filename, NC_NETCDF4 | NC_CLOBBER, &ncid);

    int x_dim_id = 0;
    nc_def_dim(ncid, "x", len_x, &x_dim_id);

    int complex_dim_id = 0;
    nc_def_dim(ncid, "complex", len_complex, &complex_dim_id);

    int var_id = 0;
    int dim_ids[2] = {x_dim_id, len_complex};
    nc_def_var(ncid, "data_complex", NC_DOUBLE, 1, dim_ids, &var_id);

    double_complex data[len_x] = {CMPLX(0., 1.), CMPLX(2., 3.), CMPLX(4., 5)};
    nc_put_var(ncid, var_id, data);

    double_complex data_out[len_x];
    size_t starts[2] = {0, 0};
    nc_get_vara(ncid, var_id, starts, NULL, data_out);

    for (size_t i = 0; i < len_x; i++) {
        print_complex("", data_out[i]);
    }

    return 0;
}
