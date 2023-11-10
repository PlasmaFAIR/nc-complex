#include <catch2/catch_test_macros.hpp>
#include <netcdf.h>

#include <array>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <ostream>
#include <stdexcept>
#include <string>

#include "nc_complex/nc_complex.h"

namespace fs = std::filesystem;

constexpr int len_x = 6;
constexpr int len_ri = 2;

constexpr std::array<std::complex<double>, len_x> double_data = {
    {{0., 1.}, {2., 3.}, {4., 5.}, {6., 7.}, {8., 9.}, {10., 11.}}
};
constexpr std::array<std::complex<float>, len_x> float_data = {
    {{0.f, 1.f}, {2.f, 3.f}, {4.f, 5.f}, {6.f, 7.f}, {8.f, 9.f}, {10.f, 11.f}}
};

constexpr std::array<std::complex<double>, len_x / 2> double_strided_data = {
    {{2., 3.}, {6., 7.}, {10., 11.}}
};
constexpr std::array<std::complex<float>, len_x / 2> float_strided_data = {
    {{2.f, 3.f}, {6.f, 7.f}, {10.f, 11.f}}
};

constexpr std::array<size_t, NC_MAX_VAR_DIMS> zeros = {0};

#define PFNC_CHECK(func)                                        \
    do {                                                        \
        if (const auto res = (func)) {                          \
            printf(                                             \
                "Bailing out in file %s, line %d, error:%s.\n", \
                __FILE__,                                       \
                __LINE__,                                       \
                nc_strerror(res)                                \
            );                                                  \
            return res;                                         \
        }                                                       \
    } while (0)

constexpr auto pfnc_complex_dir = "nc_complex_tests";

auto test_directory() {
    auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    return test_dir;
}

struct NetCDFResult {
    explicit NetCDFResult(int result) : ierr(result) {}
    operator bool() const { return ierr == NC_NOERR; }
    operator int() const { return ierr; }

private:
    int ierr{-1};
};

std::ostream& operator<<(std::ostream& os, NetCDFResult const& value) {
    os << nc_strerror(value);
    return os;
}

/// Helper macro for `REQUIRE` that takes a call to a netCDF function and checks it
/// returns ok
#define REQUIRE_NETCDF(funccall) REQUIRE(NetCDFResult{funccall})

template <std::size_t N>
inline double_complex* to_c_complex(std::array<std::complex<double>, N>& data) {
    return reinterpret_cast<double_complex*>(data.data());
}

template <std::size_t N>
inline const double_complex* to_c_complex(
    const std::array<std::complex<double>, N>& data
) {
    return reinterpret_cast<const double_complex*>(data.data());
}

template <std::size_t N>
inline float_complex* to_c_complex(std::array<std::complex<float>, N>& data) {
    return reinterpret_cast<float_complex*>(data.data());
}

/// Create a netCDF file with a variety of complex conventions
int create_file(const fs::path& filename) {
    fs::remove(filename);

    int ncid = 0;
    PFNC_CHECK(nc_create(filename.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

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

TEST_CASE("Read test file") {
    using namespace std::string_literals;

    const auto test_file = test_directory() / "test_read.nc";

    if (const auto res = create_file(test_file)) {
        const std::string error = nc_strerror(res);
        throw std::runtime_error("Couldn't create file: "s + error);
    }

    int ncid = 0;
    REQUIRE_NETCDF(nc_open(test_file.c_str(), NC_NOWRITE, &ncid));

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
}

TEST_CASE("Write complex-structure variable") {
    const auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    const auto full_filename = test_dir / "test_write_struct.nc";
    fs::remove(full_filename);

    int ncid = 0;
    REQUIRE_NETCDF(nc_create(full_filename.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

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

    REQUIRE_NETCDF(nc_put_var(ncid, var_struct_id, double_data.data()));

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
    REQUIRE_NETCDF(nc_create(full_filename.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

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
    REQUIRE_NETCDF(nc_create(full_filename.c_str(), NC_NETCDF4 | NC_CLOBBER, &ncid));

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
    REQUIRE_NETCDF(nc_create(full_filename.c_str(), NC_CLOBBER, &ncid));

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
