// -*- mode: c++ -*-

#include <netcdf.h>

#include <array>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <ostream>
#include <string_view>

#include "nc_complex/nc_complex.h"

namespace plasmafair::nc_complex::testing {

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

inline auto test_directory() {
    namespace fs = std::filesystem;
    auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    return test_dir;
}

/// Ensure the test directory exists and `filename` does *not* exist in it
inline auto test_file(std::string_view filename) {
    namespace fs = std::filesystem;
    const auto test_dir = fs::temp_directory_path() / pfnc_complex_dir;
    fs::create_directory(test_dir);
    const auto full_filename = test_dir / filename;
    fs::remove(full_filename);
    return full_filename.string();
}

struct NetCDFResult {
    explicit NetCDFResult(int result) : ierr(result) {}
    operator bool() const { return ierr == NC_NOERR; }
    operator int() const { return ierr; }

private:
    int ierr{-1};
};

inline std::ostream& operator<<(std::ostream& os, NetCDFResult const& value) {
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

NC_COMPLEX_EXPORT int create_file(const std::filesystem::path& filename);

}  // namespace plasmafair::nc_complex::testing
