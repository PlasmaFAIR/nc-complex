#include <catch2/catch_test_macros.hpp>
#include <netcdf.h>

#include <array>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>

#include "nc_complex/nc_complex.h"
#include "test_utilities.h"

namespace fs = std::filesystem;
using namespace plasmafair::nc_complex::testing;

TEST_CASE("Read test file") {
    using namespace std::string_literals;

    const auto test_file = test_directory() / "test_read.nc";

    // Only create the file for the first test case
    static bool first_run = true;
    if (first_run) {
        if (const auto res = create_file(test_file)) {
            const std::string error = nc_strerror(res);
            throw std::runtime_error("Couldn't create file: "s + error);
        }
        first_run = false;
    }

    int ncid = 0;
    REQUIRE_NETCDF(nc_open(test_file.string().c_str(), NC_NOWRITE, &ncid));

    SECTION("Using netCDF API") {
        SECTION("Reading (double) dimensional variable") {
            std::array<std::complex<double>, len_x> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri", &var_ri_id));
            REQUIRE_NETCDF(nc_get_var(ncid, var_ri_id, data_ri_out.data()));
            REQUIRE(data_ri_out == double_data);
        }

        SECTION("Reading (double) structure variable") {
            std::array<std::complex<double>, len_x> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct", &var_struct_id));
            REQUIRE_NETCDF(nc_get_var(ncid, var_struct_id, data_struct_out.data()));
            REQUIRE(data_struct_out == double_data);
        }

        SECTION("Reading (double) structure variable with long names") {
            std::array<std::complex<double>, len_x> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_long_names", &var_long_names_id));
            REQUIRE_NETCDF(
                nc_get_var(ncid, var_long_names_id, data_long_names_out.data())
            );
            REQUIRE(data_long_names_out == double_data);
        }

        SECTION("Reading (float) dimensional variable") {
            std::array<std::complex<float>, len_x> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri_float", &var_ri_id));
            REQUIRE_NETCDF(nc_get_var(ncid, var_ri_id, data_ri_out.data()));
            REQUIRE(data_ri_out == float_data);
        }

        SECTION("Reading (float) structure variable") {
            std::array<std::complex<float>, len_x> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct_float", &var_struct_id));
            REQUIRE_NETCDF(nc_get_var(ncid, var_struct_id, data_struct_out.data()));
            REQUIRE(data_struct_out == float_data);
        }

        SECTION("Reading (float) structure variable with long names") {
            std::array<std::complex<float>, len_x> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(
                nc_inq_varid(ncid, "data_long_names_float", &var_long_names_id)
            );
            REQUIRE_NETCDF(
                nc_get_var(ncid, var_long_names_id, data_long_names_out.data())
            );
            REQUIRE(data_long_names_out == float_data);
        }
    }

    SECTION("Using nc_complex untyped API") {
        constexpr std::array<size_t, 1> starts = {0};
        constexpr std::array<size_t, 1> counts = {len_x};

        SECTION("Reading (double) dimensional variable") {
            std::array<std::complex<double>, len_x> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri", &var_ri_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_ri_id));
            REQUIRE_NETCDF(pfnc_get_vara(
                ncid, var_ri_id, starts.data(), counts.data(), to_c_complex(data_ri_out)
            ));

            int var_ri_ndims = 0;
            REQUIRE_NETCDF(pfnc_inq_varndims(ncid, var_ri_id, &var_ri_ndims));
            REQUIRE(var_ri_ndims == 1);
            REQUIRE(data_ri_out == double_data);
        }

        SECTION("Reading (double) structure variable") {
            std::array<std::complex<double>, len_x> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct", &var_struct_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_struct_id));
            REQUIRE_NETCDF(pfnc_get_vara(
                ncid,
                var_struct_id,
                starts.data(),
                counts.data(),
                to_c_complex(data_struct_out)
            ));
            REQUIRE(data_struct_out == double_data);
        }

        SECTION("Reading (double) structure variable with long names") {
            std::array<std::complex<double>, len_x> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_long_names", &var_long_names_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_long_names_id));
            REQUIRE_NETCDF(pfnc_get_vara(
                ncid,
                var_long_names_id,
                starts.data(),
                counts.data(),
                to_c_complex(data_long_names_out)
            ));
            REQUIRE_NETCDF(
                nc_get_var(ncid, var_long_names_id, data_long_names_out.data())
            );
            REQUIRE(data_long_names_out == double_data);
        }

        SECTION("Reading (float) dimensional variable") {
            std::array<std::complex<float>, len_x> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri_float", &var_ri_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_ri_id));
            REQUIRE_NETCDF(pfnc_get_vara(
                ncid, var_ri_id, starts.data(), counts.data(), to_c_complex(data_ri_out)
            ));

            int var_ri_ndims = 0;
            REQUIRE_NETCDF(pfnc_inq_varndims(ncid, var_ri_id, &var_ri_ndims));
            REQUIRE(var_ri_ndims == 1);
            REQUIRE(data_ri_out == float_data);
        }

        SECTION("Reading (float) structure variable") {
            std::array<std::complex<float>, len_x> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct_float", &var_struct_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_struct_id));
            REQUIRE_NETCDF(pfnc_get_vara(
                ncid,
                var_struct_id,
                starts.data(),
                counts.data(),
                to_c_complex(data_struct_out)
            ));
            REQUIRE(data_struct_out == float_data);
        }

        SECTION("Reading (float) structure variable with long names") {
            std::array<std::complex<float>, len_x> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(
                nc_inq_varid(ncid, "data_long_names_float", &var_long_names_id)
            );
            REQUIRE(pfnc_var_is_complex(ncid, var_long_names_id));
            REQUIRE_NETCDF(pfnc_get_vara(
                ncid,
                var_long_names_id,
                starts.data(),
                counts.data(),
                to_c_complex(data_long_names_out)
            ));
            REQUIRE_NETCDF(
                nc_get_var(ncid, var_long_names_id, data_long_names_out.data())
            );
            REQUIRE(data_long_names_out == float_data);
        }
    }

    SECTION("Using nc_complex typed API") {
        SECTION("Reading (double) dimensional variable with simple get_var") {
            std::array<std::complex<double>, len_x> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri", &var_ri_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_ri_id));
            REQUIRE_NETCDF(
                pfnc_get_var_double_complex(ncid, var_ri_id, to_c_complex(data_ri_out))
            );

            int var_ri_ndims = 0;
            REQUIRE_NETCDF(pfnc_inq_varndims(ncid, var_ri_id, &var_ri_ndims));
            REQUIRE(var_ri_ndims == 1);
            REQUIRE(data_ri_out == double_data);
        }

        SECTION("Reading (double) dimensional variable") {
            std::array<std::complex<double>, len_x> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri", &var_ri_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_ri_id));
            REQUIRE_NETCDF(pfnc_get_vara_double_complex(
                ncid, var_ri_id, zeros.data(), nullptr, to_c_complex(data_ri_out)
            ));

            int var_ri_ndims = 0;
            REQUIRE_NETCDF(pfnc_inq_varndims(ncid, var_ri_id, &var_ri_ndims));
            REQUIRE(var_ri_ndims == 1);
            REQUIRE(data_ri_out == double_data);
        }

        SECTION("Reading (double) structure variable") {
            std::array<std::complex<double>, len_x> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct", &var_struct_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_struct_id));
            REQUIRE_NETCDF(pfnc_get_vara_double_complex(
                ncid,
                var_struct_id,
                zeros.data(),
                nullptr,
                to_c_complex(data_struct_out)
            ));
            REQUIRE(data_struct_out == double_data);
        }

        SECTION("Reading (double) structure variable with long names") {
            std::array<std::complex<double>, len_x> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_long_names", &var_long_names_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_long_names_id));
            REQUIRE_NETCDF(pfnc_get_vara_double_complex(
                ncid,
                var_long_names_id,
                zeros.data(),
                nullptr,
                to_c_complex(data_long_names_out)
            ));
            REQUIRE(data_long_names_out == double_data);
        }

        SECTION("Reading (float) dimensional variable") {
            std::array<std::complex<float>, len_x> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri_float", &var_ri_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_ri_id));
            REQUIRE_NETCDF(pfnc_get_vara_float_complex(
                ncid, var_ri_id, zeros.data(), nullptr, to_c_complex(data_ri_out)
            ));

            int var_ri_ndims = 0;
            REQUIRE_NETCDF(pfnc_inq_varndims(ncid, var_ri_id, &var_ri_ndims));
            REQUIRE(var_ri_ndims == 1);
            REQUIRE(data_ri_out == float_data);
        }

        SECTION("Reading (float) structure variable") {
            std::array<std::complex<float>, len_x> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct_float", &var_struct_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_struct_id));
            REQUIRE_NETCDF(pfnc_get_vara_float_complex(
                ncid,
                var_struct_id,
                zeros.data(),
                nullptr,
                to_c_complex(data_struct_out)
            ));
            REQUIRE(data_struct_out == float_data);
        }

        SECTION("Reading (float) structure variable with long names") {
            std::array<std::complex<float>, len_x> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(
                nc_inq_varid(ncid, "data_long_names_float", &var_long_names_id)
            );
            REQUIRE(pfnc_var_is_complex(ncid, var_long_names_id));
            REQUIRE_NETCDF(pfnc_get_vara_float_complex(
                ncid,
                var_long_names_id,
                zeros.data(),
                nullptr,
                to_c_complex(data_long_names_out)
            ));
            REQUIRE(data_long_names_out == float_data);
        }
    }

    SECTION("Using nc_complex typed API with counts and strides") {
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        SECTION("Reading (double) dimensional variable") {
            std::array<std::complex<double>, len_x / 2> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri", &var_ri_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_ri_id));
            REQUIRE_NETCDF(pfnc_get_vars_double_complex(
                ncid,
                var_ri_id,
                starts.data(),
                counts.data(),
                strides.data(),
                to_c_complex(data_ri_out)
            ));

            int var_ri_ndims = 0;
            REQUIRE_NETCDF(pfnc_inq_varndims(ncid, var_ri_id, &var_ri_ndims));
            REQUIRE(var_ri_ndims == 1);
            REQUIRE(data_ri_out == double_strided_data);
        }

        SECTION("Reading (double) structure variable") {
            std::array<std::complex<double>, len_x / 2> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct", &var_struct_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_struct_id));
            REQUIRE_NETCDF(pfnc_get_vars_double_complex(
                ncid,
                var_struct_id,
                starts.data(),
                counts.data(),
                strides.data(),
                to_c_complex(data_struct_out)
            ));
            REQUIRE(data_struct_out == double_strided_data);
        }

        SECTION("Reading (double) structure variable with long names") {
            std::array<std::complex<double>, len_x / 2> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_long_names", &var_long_names_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_long_names_id));
            REQUIRE_NETCDF(pfnc_get_vars_double_complex(
                ncid,
                var_long_names_id,
                starts.data(),
                counts.data(),
                strides.data(),
                to_c_complex(data_long_names_out)
            ));
            REQUIRE(data_long_names_out == double_strided_data);
        }

        SECTION("Reading (float) dimensional variable") {
            std::array<std::complex<float>, len_x / 2> data_ri_out;
            int var_ri_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri_float", &var_ri_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_ri_id));
            REQUIRE_NETCDF(pfnc_get_vars_float_complex(
                ncid,
                var_ri_id,
                starts.data(),
                counts.data(),
                strides.data(),
                to_c_complex(data_ri_out)
            ));

            int var_ri_ndims = 0;
            REQUIRE_NETCDF(pfnc_inq_varndims(ncid, var_ri_id, &var_ri_ndims));
            REQUIRE(var_ri_ndims == 1);
            REQUIRE(data_ri_out == float_strided_data);
        }

        SECTION("Reading (float) structure variable") {
            std::array<std::complex<float>, len_x / 2> data_struct_out;
            int var_struct_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct_float", &var_struct_id));
            REQUIRE(pfnc_var_is_complex(ncid, var_struct_id));
            REQUIRE_NETCDF(pfnc_get_vars_float_complex(
                ncid,
                var_struct_id,
                starts.data(),
                counts.data(),
                strides.data(),
                to_c_complex(data_struct_out)
            ));
            REQUIRE(data_struct_out == float_strided_data);
        }

        SECTION("Reading (float) structure variable with long names") {
            std::array<std::complex<float>, len_x / 2> data_long_names_out;
            int var_long_names_id = 0;
            REQUIRE_NETCDF(
                nc_inq_varid(ncid, "data_long_names_float", &var_long_names_id)
            );
            REQUIRE(pfnc_var_is_complex(ncid, var_long_names_id));
            REQUIRE_NETCDF(pfnc_get_vars_float_complex(
                ncid,
                var_long_names_id,
                starts.data(),
                counts.data(),
                strides.data(),
                to_c_complex(data_long_names_out)
            ));
            REQUIRE(data_long_names_out == float_strided_data);
        }
    }

    REQUIRE_NETCDF(nc_close(ncid));
}

TEST_CASE("Write complex-structure variable") {
    const auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    const auto full_filename = test_dir / "test_write_struct.nc";
    fs::remove(full_filename);

    int ncid = 0;
    REQUIRE_NETCDF(nc_create(full_filename.string().c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

    int x_dim_id = 0;
    REQUIRE_NETCDF(nc_def_dim(ncid, "x", len_x, &x_dim_id));

    int type_id{};
    REQUIRE_NETCDF(pfnc_get_double_complex_typeid(ncid, &type_id));

    SECTION("Check getting type is idempotent") {
        int type_id2{};
        REQUIRE_NETCDF(pfnc_get_double_complex_typeid(ncid, &type_id2));
        REQUIRE(type_id == type_id2);
    }

    SECTION("Check base type of compound type") {
        int base_type_id{};
        REQUIRE_NETCDF(pfnc_complex_base_type(ncid, type_id, &base_type_id));
        REQUIRE(base_type_id == NC_DOUBLE);
    }

    int var_struct_id = 0;
    const std::array<int, 1> dim_struct_ids{{x_dim_id}};
    REQUIRE_NETCDF(nc_def_var(
        ncid, "data_struct", type_id, 1, dim_struct_ids.data(), &var_struct_id
    ));

    SECTION("Check base type of complex variable") {
        int base_type_id{};
        REQUIRE_NETCDF(
            pfnc_inq_var_complex_base_type(ncid, var_struct_id, &base_type_id)
        );
        REQUIRE(base_type_id == NC_DOUBLE);
    }

    REQUIRE_NETCDF(
        pfnc_put_var_double_complex(ncid, var_struct_id, to_c_complex(double_data))
    );

    SECTION("Reading structure variable") {
        std::array<std::complex<double>, len_x> data_struct_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_struct_id, zeros.data(), nullptr, to_c_complex(data_struct_out)
        ));
        REQUIRE(data_struct_out == double_data);
    }

    SECTION("Writing slice") {
        constexpr std::array<std::complex<double>, len_x / 2> slice_data_in = {
            {{-2., -3.}, {-6., -7.}, {-10., -11.}}
        };
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_put_vars_double_complex(
            ncid,
            var_struct_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_in)
        ));
        REQUIRE_NETCDF(nc_sync(ncid));

        std::array<std::complex<double>, len_x> data_struct_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_struct_id, zeros.data(), nullptr, to_c_complex(data_struct_out)
        ));

        constexpr std::array<std::complex<double>, len_x> expected_sliced_data = {
            {{0., 1.}, {-2., -3.}, {4., 5.}, {-6., -7.}, {8., 9.}, {-10., -11.}}
        };

        REQUIRE(data_struct_out == expected_sliced_data);
    }

    SECTION("Reading slice") {
        constexpr std::array<std::complex<double>, len_x / 2> expected_sliced_data = {
            {{2., 3.}, {6., 7.}, {10., 11.}}
        };

        std::array<std::complex<double>, len_x / 2> slice_data_out{};
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_get_vars_double_complex(
            ncid,
            var_struct_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_out)
        ));
        REQUIRE(slice_data_out == expected_sliced_data);
    }

    REQUIRE_NETCDF(nc_close(ncid));
}

TEST_CASE("Write complex-dimensioned variable") {
    const auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    const auto full_filename = test_dir / "test_write_dim.nc";
    fs::remove(full_filename);

    int ncid = 0;
    REQUIRE_NETCDF(nc_create(full_filename.string().c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

    int x_dim_id = 0;
    REQUIRE_NETCDF(nc_def_dim(ncid, "x", len_x, &x_dim_id));

    int complex_dim_id = 0;
    REQUIRE_NETCDF(pfnc_get_complex_dim(ncid, &complex_dim_id));

    SECTION("Check getting dimension is idempotent") {
        int ri_dim_id2 = 0;
        REQUIRE_NETCDF(pfnc_get_complex_dim(ncid, &ri_dim_id2));
        REQUIRE(complex_dim_id == ri_dim_id2);
    }

    SECTION("Check other complex dimension names are recognised") {
        int ri_dim_id = 0;
        REQUIRE_NETCDF(nc_def_dim(ncid, "ri", 2, &ri_dim_id));
        REQUIRE(pfnc_is_complex_dim(ncid, ri_dim_id));

        int complex_dim_id2 = 0;
        REQUIRE_NETCDF(nc_def_dim(ncid, "complex", 2, &complex_dim_id2));
        REQUIRE(pfnc_is_complex_dim(ncid, complex_dim_id2));
    }

    int var_id = 0;
    const std::array<int, 2> dim_ids{{x_dim_id, complex_dim_id}};
    REQUIRE_NETCDF(nc_def_var(ncid, "data_dim", NC_DOUBLE, 2, dim_ids.data(), &var_id));

    SECTION("Check base type of complex variable") {
        int base_type_id{};
        REQUIRE_NETCDF(pfnc_inq_var_complex_base_type(ncid, var_id, &base_type_id));
        REQUIRE(base_type_id == NC_DOUBLE);
    }

    SECTION("Chunking") {
        constexpr std::array<std::size_t, 1> chunk_sizes{{len_x / 2}};
        REQUIRE_NETCDF(
            pfnc_def_var_chunking(ncid, var_id, NC_CHUNKED, chunk_sizes.data())
        );

        std::array<std::size_t, 1> chunk_sizes_out{};
        REQUIRE_NETCDF(
            pfnc_inq_var_chunking(ncid, var_id, nullptr, chunk_sizes_out.data())
        );
        REQUIRE(chunk_sizes_out == chunk_sizes);
    }

    REQUIRE_NETCDF(nc_put_var(ncid, var_id, double_data.data()));

    SECTION("Reading variable") {
        std::array<std::complex<double>, len_x> data_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_id, zeros.data(), nullptr, to_c_complex(data_out)
        ));
        REQUIRE(data_out == double_data);
    }

    SECTION("Writing slice") {
        constexpr std::array<std::complex<double>, len_x / 2> slice_data_in = {
            {{-2., -3.}, {-6., -7.}, {-10., -11.}}
        };
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_put_vars_double_complex(
            ncid,
            var_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_in)
        ));
        REQUIRE_NETCDF(nc_sync(ncid));

        std::array<std::complex<double>, len_x> data_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_id, zeros.data(), nullptr, to_c_complex(data_out)
        ));

        constexpr std::array<std::complex<double>, len_x> expected_sliced_data = {
            {{0., 1.}, {-2., -3.}, {4., 5.}, {-6., -7.}, {8., 9.}, {-10., -11.}}
        };

        REQUIRE(data_out == expected_sliced_data);
    }

    SECTION("Reading slice") {
        constexpr std::array<std::complex<double>, len_x / 2> expected_sliced_data = {
            {{2., 3.}, {6., 7.}, {10., 11.}}
        };

        std::array<std::complex<double>, len_x / 2> slice_data_out{};
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_get_vars_double_complex(
            ncid,
            var_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_out)
        ));
        REQUIRE(slice_data_out == expected_sliced_data);
    }

    REQUIRE_NETCDF(nc_close(ncid));
}

TEST_CASE("Write custom-type variable") {
    const auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    const auto full_filename = test_dir / "test_write_custom_type.nc";
    fs::remove(full_filename);

    int ncid = 0;
    REQUIRE_NETCDF(nc_create(full_filename.string().c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

    int x_dim_id = 0;
    REQUIRE_NETCDF(nc_def_dim(ncid, "x", len_x, &x_dim_id));

    int var_id = 0;
    const std::array<int, 1> dim_ids{{x_dim_id}};
    REQUIRE_NETCDF(
        pfnc_def_var(ncid, "data", PFNC_DOUBLE_COMPLEX, 1, dim_ids.data(), &var_id)
    );

    SECTION("Check base type of complex variable") {
        int base_type_id{};
        REQUIRE_NETCDF(pfnc_inq_var_complex_base_type(ncid, var_id, &base_type_id));
        REQUIRE(base_type_id == NC_DOUBLE);
    }

    REQUIRE_NETCDF(nc_put_var(ncid, var_id, double_data.data()));

    SECTION("Reading variable") {
        std::array<std::complex<double>, len_x> data_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_id, zeros.data(), nullptr, to_c_complex(data_out)
        ));
        REQUIRE(data_out == double_data);
    }

    SECTION("Writing slice") {
        constexpr std::array<std::complex<double>, len_x / 2> slice_data_in = {
            {{-2., -3.}, {-6., -7.}, {-10., -11.}}
        };
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_put_vars_double_complex(
            ncid,
            var_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_in)
        ));
        REQUIRE_NETCDF(nc_sync(ncid));

        std::array<std::complex<double>, len_x> data_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_id, zeros.data(), nullptr, to_c_complex(data_out)
        ));

        constexpr std::array<std::complex<double>, len_x> expected_sliced_data = {
            {{0., 1.}, {-2., -3.}, {4., 5.}, {-6., -7.}, {8., 9.}, {-10., -11.}}
        };

        REQUIRE(data_out == expected_sliced_data);
    }

    SECTION("Reading slice") {
        constexpr std::array<std::complex<double>, len_x / 2> expected_sliced_data = {
            {{2., 3.}, {6., 7.}, {10., 11.}}
        };

        std::array<std::complex<double>, len_x / 2> slice_data_out{};
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_get_vars_double_complex(
            ncid,
            var_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_out)
        ));
        REQUIRE(slice_data_out == expected_sliced_data);
    }

    REQUIRE_NETCDF(nc_close(ncid));
}

TEST_CASE("Write custom-type variable (netCDF3)") {
    const auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    const auto full_filename = test_dir / "test_write_custom_type_netcdf3.nc";
    fs::remove(full_filename);

    int ncid = 0;
    REQUIRE_NETCDF(nc_create(full_filename.string().c_str(), NC_CLOBBER, &ncid));

    int x_dim_id = 0;
    REQUIRE_NETCDF(nc_def_dim(ncid, "x", len_x, &x_dim_id));

    int var_id = 0;
    const std::array<int, 1> dim_ids{{x_dim_id}};
    REQUIRE_NETCDF(
        pfnc_def_var(ncid, "data", PFNC_DOUBLE_COMPLEX, 1, dim_ids.data(), &var_id)
    );

    SECTION("Check base type of complex variable") {
        int base_type_id{};
        REQUIRE_NETCDF(pfnc_inq_var_complex_base_type(ncid, var_id, &base_type_id));
        REQUIRE(base_type_id == NC_DOUBLE);
    }

    REQUIRE_NETCDF(nc_enddef(ncid));

    REQUIRE_NETCDF(nc_put_var(ncid, var_id, double_data.data()));

    SECTION("Reading variable") {
        std::array<std::complex<double>, len_x> data_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_id, zeros.data(), nullptr, to_c_complex(data_out)
        ));
        REQUIRE(data_out == double_data);
    }

    SECTION("Writing slice") {
        constexpr std::array<std::complex<double>, len_x / 2> slice_data_in = {
            {{-2., -3.}, {-6., -7.}, {-10., -11.}}
        };
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_put_vars_double_complex(
            ncid,
            var_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_in)
        ));
        REQUIRE_NETCDF(nc_sync(ncid));

        std::array<std::complex<double>, len_x> data_out;
        REQUIRE_NETCDF(pfnc_get_vara_double_complex(
            ncid, var_id, zeros.data(), nullptr, to_c_complex(data_out)
        ));

        constexpr std::array<std::complex<double>, len_x> expected_sliced_data = {
            {{0., 1.}, {-2., -3.}, {4., 5.}, {-6., -7.}, {8., 9.}, {-10., -11.}}
        };

        REQUIRE(data_out == expected_sliced_data);
    }

    SECTION("Reading slice") {
        constexpr std::array<std::complex<double>, len_x / 2> expected_sliced_data = {
            {{2., 3.}, {6., 7.}, {10., 11.}}
        };

        std::array<std::complex<double>, len_x / 2> slice_data_out{};
        constexpr std::array<std::size_t, 1> starts = {{1}};
        constexpr std::array<std::size_t, 1> counts = {{len_x / 2}};
        constexpr std::array<std::ptrdiff_t, 1> strides = {{2}};

        REQUIRE_NETCDF(pfnc_get_vars_double_complex(
            ncid,
            var_id,
            starts.data(),
            counts.data(),
            strides.data(),
            to_c_complex(slice_data_out)
        ));
        REQUIRE(slice_data_out == expected_sliced_data);
    }

    REQUIRE_NETCDF(nc_close(ncid));
}

TEST_CASE("Write custom-type variable with pre-existing type") {
    const auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    const auto full_filename = test_dir / "test_write_preexisting_type.nc";
    fs::remove(full_filename);

    int ncid = 0;
    REQUIRE_NETCDF(nc_create(full_filename.string().c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

    int x_dim_id = 0;
    REQUIRE_NETCDF(nc_def_dim(ncid, "x", len_x, &x_dim_id));

    nc_type complex_type_id;
    REQUIRE_NETCDF(nc_def_compound(
        ncid, sizeof(std::complex<double>), "custom_complex", &complex_type_id
    ));
    REQUIRE_NETCDF(nc_insert_compound(ncid, complex_type_id, "r", 0, NC_DOUBLE));
    REQUIRE_NETCDF(
        nc_insert_compound(ncid, complex_type_id, "i", sizeof(double), NC_DOUBLE)
    );

    int var_id = 0;
    const std::array<int, 1> dim_ids{{x_dim_id}};
    REQUIRE_NETCDF(
        pfnc_def_var(ncid, "data", PFNC_DOUBLE_COMPLEX, 1, dim_ids.data(), &var_id)
    );

    int base_type_id{};
    REQUIRE_NETCDF(pfnc_inq_var_complex_base_type(ncid, var_id, &base_type_id));
    REQUIRE(base_type_id == NC_DOUBLE);

    nc_type var_type_id{};
    REQUIRE_NETCDF(nc_inq_vartype(ncid, var_id, &var_type_id));
    REQUIRE(var_type_id == complex_type_id);

    REQUIRE_NETCDF(nc_put_var(ncid, var_id, double_data.data()));

    std::array<std::complex<double>, len_x> data_out;
    REQUIRE_NETCDF(pfnc_get_vara_double_complex(
        ncid, var_id, zeros.data(), nullptr, to_c_complex(data_out)
    ));
    REQUIRE(data_out == double_data);

    REQUIRE_NETCDF(nc_close(ncid));
}
