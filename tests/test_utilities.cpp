#include "test_utilities.h"

namespace fs = std::filesystem;
namespace plasmafair::nc_complex::testing {

/// Create a netCDF file with a variety of complex conventions

int create_file(const fs::path& filename) {
    fs::remove(filename);

    int ncid = 0;
    PFNC_CHECK(nc_create(filename.string().c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

    int x_dim_id = 0;
    PFNC_CHECK(nc_def_dim(ncid, "x", len_x, &x_dim_id));

    int ri_dim_id = 0;
    PFNC_CHECK(nc_def_dim(ncid, "ri", len_ri, &ri_dim_id));

    const std::array<int, 2> dimids{{x_dim_id, ri_dim_id}};
    {
        ////////////////////
        // Double variables

        int var_ri = 0;
        PFNC_CHECK(nc_def_var(ncid, "data_ri", NC_DOUBLE, 2, dimids.data(), &var_ri));

        int type_id = 0;
        PFNC_CHECK(nc_def_compound(
            ncid, sizeof(std::complex<double>), "my_double_complex", &type_id
        ));
        PFNC_CHECK(nc_insert_compound(ncid, type_id, "r", 0, NC_DOUBLE));
        PFNC_CHECK(nc_insert_compound(ncid, type_id, "i", sizeof(double), NC_DOUBLE));
        int var_struct_id = 0;
        const std::array<int, 1> dim_struct_ids{{x_dim_id}};
        PFNC_CHECK(nc_def_var(
            ncid, "data_struct", type_id, 1, dim_struct_ids.data(), &var_struct_id
        ));

        int long_names_type_id = 0;
        PFNC_CHECK(nc_def_compound(
            ncid,
            sizeof(std::complex<double>),
            "long_names_double_complex",
            &long_names_type_id
        ));
        PFNC_CHECK(nc_insert_compound(ncid, long_names_type_id, "Real", 0, NC_DOUBLE));
        PFNC_CHECK(nc_insert_compound(
            ncid, long_names_type_id, "Imag", sizeof(double), NC_DOUBLE
        ));
        int var_long_names_id = 0;
        const std::array<int, 1> dim_long_names_ids{{x_dim_id}};
        PFNC_CHECK(nc_def_var(
            ncid,
            "data_long_names",
            long_names_type_id,
            1,
            dim_long_names_ids.data(),
            &var_long_names_id
        ));

        PFNC_CHECK(nc_put_var(ncid, var_ri, double_data.data()));
        PFNC_CHECK(nc_put_var(ncid, var_struct_id, double_data.data()));
        PFNC_CHECK(nc_put_var(ncid, var_long_names_id, double_data.data()));
    }

    {
        ////////////////////
        // Float variables

        int var_ri = 0;
        PFNC_CHECK(
            nc_def_var(ncid, "data_ri_float", NC_FLOAT, 2, dimids.data(), &var_ri)
        );

        int type_id = 0;
        PFNC_CHECK(nc_def_compound(
            ncid, sizeof(std::complex<float>), "my_float_complex", &type_id
        ));
        PFNC_CHECK(nc_insert_compound(ncid, type_id, "r", 0, NC_FLOAT));
        PFNC_CHECK(nc_insert_compound(ncid, type_id, "i", sizeof(float), NC_FLOAT));
        int var_struct_id = 0;
        const std::array<int, 1> dim_struct_ids{{x_dim_id}};
        PFNC_CHECK(nc_def_var(
            ncid, "data_struct_float", type_id, 1, dim_struct_ids.data(), &var_struct_id
        ));

        int long_names_type_id = 0;
        PFNC_CHECK(nc_def_compound(
            ncid,
            sizeof(std::complex<float>),
            "long_names_float_complex",
            &long_names_type_id
        ));
        PFNC_CHECK(nc_insert_compound(ncid, long_names_type_id, "Real", 0, NC_FLOAT));
        PFNC_CHECK(nc_insert_compound(
            ncid, long_names_type_id, "Imag", sizeof(float), NC_FLOAT
        ));
        int var_long_names_id = 0;
        const std::array<int, 1> dim_long_names_ids{{x_dim_id}};
        PFNC_CHECK(nc_def_var(
            ncid,
            "data_long_names_float",
            long_names_type_id,
            1,
            dim_long_names_ids.data(),
            &var_long_names_id
        ));

        PFNC_CHECK(nc_put_var(ncid, var_ri, float_data.data()));
        PFNC_CHECK(nc_put_var(ncid, var_struct_id, float_data.data()));
        PFNC_CHECK(nc_put_var(ncid, var_long_names_id, float_data.data()));
    }

    PFNC_CHECK(nc_close(ncid));

    return 0;
}

}  // namespace plasmafair::nc_complex::testing
