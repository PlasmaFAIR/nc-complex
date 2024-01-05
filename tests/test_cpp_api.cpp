#include <catch2/catch_test_macros.hpp>

#include <complex>
#include <cstddef>
#include <filesystem>
#include <netcdf>
#include <stdexcept>
#include <string>
#include <vector>

#include "nc_complex/nc_complex_cpp.h"
#include "test_utilities.h"

using namespace nc_complex::testing;
using namespace std::string_literals;

TEST_CASE("Read test file") {
    const auto test_file = test_directory() / "test_read_cpp.nc";

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

    SECTION("Use NcVar directly") {
        REQUIRE_NETCDF(nc_open(test_file.c_str(), NC_NOWRITE, &ncid));

        SECTION("Reading (double) dimensional variable") {
            int var_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_ri", &var_id));

            const nc_complex::NcVar var{ncid, var_id};

            REQUIRE(var.isComplex());
            REQUIRE(!var.isComplexType());
            REQUIRE(var.hasComplexDim());

            std::array<std::complex<double>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == double_data);
            REQUIRE(var.getDimCount() == 1);
            auto dims = var.getDims();
            REQUIRE(dims[0].getName() == "x");
        }

        SECTION("Reading (double) structure variable") {
            int var_id = 0;
            REQUIRE_NETCDF(nc_inq_varid(ncid, "data_struct", &var_id));

            const nc_complex::NcVar var{ncid, var_id};

            REQUIRE(var.isComplex());
            REQUIRE(var.isComplexType());
            REQUIRE(!var.hasComplexDim());

            std::array<std::complex<double>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == double_data);
            REQUIRE(var.getDimCount() == 1);
            auto dims = var.getDims();
            REQUIRE(dims[0].getName() == "x");
        }

        REQUIRE_NETCDF(nc_close(ncid));
    }

    SECTION("Use NcGroup") {
        const nc_complex::NcFile nc_file{test_file, nc_complex::NcFile::FileMode::read};
        REQUIRE_FALSE(nc_file.isNull());

        const auto vars = nc_file.getVars();
        REQUIRE(vars.size() == 6);

        SECTION("Reading (double) structure variable") {
            const auto var = nc_file.getVar("data_ri");

            REQUIRE(var.isComplex());
            REQUIRE(!var.isComplexType());
            REQUIRE(var.hasComplexDim());

            std::array<std::complex<double>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == double_data);
            REQUIRE(var.getDimCount() == 1);
            auto dims = var.getDims();
            REQUIRE(dims[0].getName() == "x");
        }

        SECTION("Reading (double) structure variable") {
            const auto var = nc_file.getVar("data_struct");

            REQUIRE(var.isComplex());
            REQUIRE(var.isComplexType());
            REQUIRE(!var.hasComplexDim());

            std::array<std::complex<double>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == double_data);
        }

        SECTION("Reading (double) structure variable") {
            const auto var = nc_file.getVar("data_long_names");

            REQUIRE(var.isComplex());
            REQUIRE(var.isComplexType());
            REQUIRE(!var.hasComplexDim());

            std::array<std::complex<double>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == double_data);
        }

        SECTION("Reading (float) structure variable") {
            const auto var = nc_file.getVar("data_ri_float");

            REQUIRE(var.isComplex());
            REQUIRE(!var.isComplexType());
            REQUIRE(var.hasComplexDim());

            std::array<std::complex<float>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == float_data);
            REQUIRE(var.getDimCount() == 1);
            auto dims = var.getDims();
            REQUIRE(dims[0].getName() == "x");
        }

        SECTION("Reading (float) structure variable") {
            const auto var = nc_file.getVar("data_struct_float");

            REQUIRE(var.isComplex());
            REQUIRE(var.isComplexType());
            REQUIRE(!var.hasComplexDim());

            std::array<std::complex<float>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == float_data);
        }

        SECTION("Reading (float) structure variable") {
            const auto var = nc_file.getVar("data_long_names_float");

            REQUIRE(var.isComplex());
            REQUIRE(var.isComplexType());
            REQUIRE(!var.hasComplexDim());

            std::array<std::complex<float>, len_x> data_out;
            var.getVar(data_out.data());

            REQUIRE(data_out == float_data);
        }
    }
}

TEST_CASE("Write complex-structure variable") {
    const auto full_filename = test_file("test_write_struct_cpp.nc");

    nc_complex::NcFile nc_file{full_filename, nc_complex::NcFile::FileMode::newFile};

    const nc_complex::NcDim x_dim = nc_file.addDim("x", len_x);

    auto var = nc_file.addVar("data_struct", nc_complex::NcDoubleComplex{}, x_dim);

    REQUIRE(var.getDimCount() == 1);
    REQUIRE(var.isComplex());

    var.putVar(double_data.data());
    nc_file.sync();

    std::array<std::complex<double>, len_x> data_out;
    var.getVar(data_out.data());
    REQUIRE(data_out == double_data);

    SECTION("Writing slice") {
        constexpr std::array<std::complex<double>, len_x / 2> slice_data_in = {
            {{-2., -3.}, {-6., -7.}, {-10., -11.}}
        };
        var.putVar({1}, {len_x / 2}, {2}, slice_data_in.data());
        nc_file.sync();

        std::array<std::complex<double>, len_x> data_struct_out;
        var.getVar(data_struct_out.data());

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

        var.getVar({1}, {len_x / 2}, {2}, slice_data_out.data());

        REQUIRE(slice_data_out == expected_sliced_data);
    }
}

TEST_CASE("Write complex-dimensioned variable") {
    const auto full_filename = test_file("test_write_dim_cpp.nc");

    nc_complex::NcFile nc_file{full_filename, nc_complex::NcFile::FileMode::newFile};

    const nc_complex::NcDim x_dim = nc_file.addDim("x", len_x);

    auto var = nc_file.addVar("data_ri", nc_complex::NcDoubleComplexDim{}, x_dim);

    REQUIRE(var.getDimCount() == 1);
    REQUIRE(var.isComplex());

    var.putVar(double_data.data());

    std::array<std::complex<double>, len_x> data_out;
    var.getVar(data_out.data());
    REQUIRE(data_out == double_data);

    SECTION("Writing slice") {
        constexpr std::array<std::complex<double>, len_x / 2> slice_data_in = {
            {{-2., -3.}, {-6., -7.}, {-10., -11.}}
        };
        var.putVar({1}, {len_x / 2}, {2}, slice_data_in.data());
        nc_file.sync();

        std::array<std::complex<double>, len_x> data_struct_out;
        var.getVar(data_struct_out.data());

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

        var.getVar({1}, {len_x / 2}, {2}, slice_data_out.data());

        REQUIRE(slice_data_out == expected_sliced_data);
    }
}

TEST_CASE("Base library functionality") {
    // These are a bunch of tests copied from the base netcdf-cxx4 test suite

    using namespace nc_complex;

    SECTION("Classic format") {
        const auto full_filename = test_file("test_base_classic.nc");

        {
            const NcFile ncFile{full_filename, NcFile::replace, NcFile::classic};

            [[maybe_unused]] const NcDim dim1 = ncFile.addDim("dim1", 11);
            [[maybe_unused]] const NcDim dim2 = ncFile.addDim("dim2");
            [[maybe_unused]] const NcDim dim3 = ncFile.addDim("dim3", 13);

            const NcVar var_gw = ncFile.addVar("George_Washington", ncInt, dim1);
            constexpr std::array arr = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
            var_gw.putVar(arr.data());

            const std::vector<NcDim> dimArray{dim2, dim1};
            [[maybe_unused]] const NcVar varA1_3 =
                ncFile.addVar("varA1_3", ncInt, dimArray);
        }

        {
            const NcFile ncFile{full_filename, NcFile::read};
            REQUIRE(ncFile.getVarCount() == 2);
        }

        {
            const NcFile ncFile{full_filename, NcFile::write};
            ncFile.putAtt("name", "value");
        }
    }

    SECTION("NetCDF4") {
        // Just tests a subset of the interface

        const auto filename = test_file("test_base_var_group_structure.nc");
        const NcFile ncFile{filename, NcFile::replace};

        const NcGroup groupA = ncFile.addGroup("groupA");
        const NcGroup groupA0 = ncFile.addGroup("groupA0");
        const NcGroup groupB = groupA.addGroup("groupB");
        const NcGroup groupC = groupA.addGroup("groupC");

        const NcDim dim1 = ncFile.addDim("dim1", 2);
        const NcDim dim2 = ncFile.addDim("dim2");
        const NcDim dim3 = ncFile.addDim("dim3", 3);
        const NcDim dim4 = groupA0.addDim("dim4", 5);

        ncFile.addVar("var_1", ncByte);
        groupA.addVar("varA_1", ncChar, dim1);
        groupA.addVar("varA_2", ncShort, {dim1, dim2});
        groupA0.addVar("varA0_1", ncInt, {dim1, dim2});
        groupA0.addVar("varA0_2", ncFloat, {dim1, dim3, dim4});
        groupA0.addVar("varA0_3", ncDouble, {dim1, dim3, dim4, dim2});
        groupB.addVar("varB_1", ncUbyte, dim1);
        groupB.addVar("varB_2", ncUshort, dim1);
        groupB.addVar("varB_3", ncUint, dim1);
        groupB.addVar("varB_4", ncInt64, dim1);
        groupC.addVar("varC_1", ncUint64, dim1);
        groupC.addVar("varC_2", ncString, dim1);
        groupC.addVar("varC_3", "byte", "dim1");
        groupC.addVar("varC_4", "double", {"dim2"s, "dim3"s});
        groupC.addVar("varC_5", ncByte, dim1);

        REQUIRE(ncFile.getVarCount() == 1);
        REQUIRE(ncFile.getVarCount(NcGroup::Current) == 1);
        REQUIRE(ncFile.getVarCount(NcGroup::Parents) == 0);
        REQUIRE(ncFile.getVarCount(NcGroup::Children) == 14);
        REQUIRE(ncFile.getVarCount(NcGroup::ParentsAndCurrent) == 1);
        REQUIRE(ncFile.getVarCount(NcGroup::ChildrenAndCurrent) == 15);
        REQUIRE(ncFile.getVarCount(NcGroup::All) == 15);

        REQUIRE(groupA.getVarCount() == 2);
        REQUIRE(groupA.getVarCount(NcGroup::Current) == 2);
        REQUIRE(groupA.getVarCount(NcGroup::Parents) == 1);
        REQUIRE(groupA.getVarCount(NcGroup::Children) == 9);
        REQUIRE(groupA.getVarCount(NcGroup::ParentsAndCurrent) == 3);
        REQUIRE(groupA.getVarCount(NcGroup::ChildrenAndCurrent) == 11);
        REQUIRE(groupA.getVarCount(NcGroup::All) == 12);

        const auto file_vars = ncFile.getVars();
        REQUIRE(file_vars.size() == 1);
        REQUIRE(file_vars.find("var_1") != file_vars.end());
        REQUIRE(file_vars.find("varA_1") == file_vars.end());
        REQUIRE(file_vars.find("varA0_1") == file_vars.end());
        REQUIRE(file_vars.find("varB_3") == file_vars.end());
        REQUIRE(file_vars.find("varC_5") == file_vars.end());

        REQUIRE(ncFile.getVars(NcGroup::Parents).empty());

        constexpr std::size_t N = 10;
        const auto int_attribute = fill_sequence<int, N>();
        const NcGroupAtt att1_1 =
            ncFile.putAtt("att1_1", ncByte, N, int_attribute.data());
        const NcGroupAtt att1_2 =
            ncFile.putAtt("att1_2", ncInt, N, int_attribute.data());
        const NcGroupAtt att1_3 =
            ncFile.putAtt("att1_3", ncDouble, N, int_attribute.data());
        const NcGroupAtt att2_1 =
            groupA.putAtt("att2_1", ncUshort, N, int_attribute.data());
        const NcGroupAtt att2_2 =
            groupA.putAtt("att2_2", ncFloat, N, int_attribute.data());

        REQUIRE(ncFile.getAttCount() == 3);
        REQUIRE(groupA.getAttCount(NcGroup::ParentsAndCurrent) == 5);

        REQUIRE(att1_1.getAttLength() == N);
        REQUIRE(att1_2.getName() == "att1_2");
        REQUIRE(att1_3.getType() == ncDouble);
        REQUIRE(att2_1.getParentGroup() == groupA);

        const auto expected_float_att = fill_sequence<float, N>();
        std::vector<float> float_att_out(N);
        att2_2.getValues(float_att_out.data());
        REQUIRE(float_att_out == expected_float_att);
    }
}
