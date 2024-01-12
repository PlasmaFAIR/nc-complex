#include "nc_complex/nc_complex_cpp.h"

#include <netcdf.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "nc_complex/nc_complex.h"

using namespace nc_complex::exceptions;

#define PFNC_CPP_CHECK(retcode) ncCheck(retcode, __FILE__, __LINE__)
#define PFNC_CPP_NullType(message) NcNullType(message, __FILE__, __LINE__)
#define PFNC_CPP_NullDim(message) NcNullDim(message, __FILE__, __LINE__)

namespace nc_complex {

std::multimap<std::string, NcGroup> NcGroup::getGroups(NcGroup::GroupLocation location
) const {
    const auto base_groups = netCDF::NcGroup::getGroups(location);
    std::multimap<std::string, NcGroup> complex_groups;
    for (const auto& group : base_groups) {
        complex_groups.emplace(group);
    }
    return complex_groups;
}
std::set<NcGroup> NcGroup::getGroups(
    const std::string& name, NcGroup::GroupLocation location
) const {
    const auto base_groups = netCDF::NcGroup::getGroups(name, location);
    std::set<NcGroup> complex_groups;
    for (const auto& group : base_groups) {
        complex_groups.emplace(group);
    }
    return complex_groups;
}
NcGroup NcGroup::getGroup(const std::string& name, NcGroup::GroupLocation location)
    const {
    return NcGroup{netCDF::NcGroup::getGroup(name, location)};
}
NcGroup NcGroup::addGroup(const std::string& name) const {
    return NcGroup{netCDF::NcGroup::addGroup(name)};
}

std::multimap<std::string, NcVar> NcGroup::getVars(NcGroup::Location location) const {
    const auto base_vars = netCDF::NcGroup::getVars(location);
    std::multimap<std::string, NcVar> complex_vars;
    for (const auto& var : base_vars) {
        complex_vars.emplace(var.first, var.second);
    }
    return complex_vars;
}

std::set<NcVar> NcGroup::getVars(const std::string& name, NcGroup::Location location)
    const {
    const auto base_vars = netCDF::NcGroup::getVars(name, location);
    std::set<NcVar> complex_vars;
    for (const auto& var : base_vars) {
        complex_vars.emplace(var);
    }
    return complex_vars;
}

NcVar NcGroup::getVar(const std::string& name, NcGroup::Location location) const {
    return NcVar{netCDF::NcGroup::getVar(name, location)};
}

NcVar NcGroup::addVar(
    const std::string& name, const std::string& typeName, const std::string& dimName
) const {
    return addVar(
        name,
        getType(typeName, NcGroup::ParentsAndCurrent),
        getDim(dimName, NcGroup::ParentsAndCurrent)
    );
}
NcVar NcGroup::addVar(const std::string& name, const NcType& ncType, const NcDim& ncDim)
    const {
    return addVar(name, ncType, std::vector{ncDim});
}
NcVar NcGroup::addVar(
    const std::string& name, const NcAtomicType& ncType, const NcDim& ncDim
) const {
    return addVar(name, ncType, std::vector{ncDim});
}

NcVar NcGroup::addVar(const std::string& name, const NcType& ncType) const {
    return NcGroup::addVar(name, ncType, std::vector<NcDim>{});
}

NcVar NcGroup::addVar(const std::string& name, const NcAtomicType& ncType) const {
    return NcGroup::addVar(name, ncType, std::vector<NcDim>{});
}

/// Check the given dimension is defined either in group or one of its
/// direct ancestors.
/// Returns the dimension ID
auto checkDimensionDefinedInGroupAncestors(const NcGroup& group, const NcDim& dim) {
    const auto result = group.getDim(dim.getName(), NcGroup::ParentsAndCurrent);
    if (result.isNull()) {
        throw PFNC_CPP_NullDim(
            "Attempt to invoke NcGroup::addVar failed: NcDim must be defined in either the current group or a parent group"
        );
    }
}

NcVar NcGroup::addVar(
    const std::string& name,
    const std::string& typeName,
    const std::vector<std::string>& dimNames
) const {
    std::vector<NcDim> dims;
    dims.reserve(dimNames.size());
    std::transform(
        dimNames.begin(),
        dimNames.end(),
        std::back_inserter(dims),
        [this](const auto& dim_name) {
            return getDim(dim_name, NcGroup::ParentsAndCurrent);
        }
    );
    return addVar(name, getType(typeName), dims);
}

/// Checks all dimensions in input vector are non-null and are defined
/// in the current group, or one of its direct ancestors.
/// Returns a vector of dimension IDs
auto checkDimensionsAreValid(const NcGroup& group, const std::vector<NcDim>& dims) {
    std::vector<int> dim_ids;
    dim_ids.reserve(dims.size());
    for (const auto& dim : dims) {
        if (dim.isNull()) {
            throw PFNC_CPP_NullDim(
                "Attempt to invoke NcGroup::addVar with a Null NcDim object"
            );
        }
        checkDimensionDefinedInGroupAncestors(group, dim);
        dim_ids.push_back(dim.getId());
    }
    return dim_ids;
}

NcVar NcGroup::addVar(
    const std::string& name, const NcType& ncType, const std::vector<NcDim>& ncDimVector
) const {
    ncCheckDefineMode(myId);

    // check NcType object is valid
    if (ncType.isNull()) {
        throw PFNC_CPP_NullType(
            "Attempt to invoke NcGroup::addVar with a Null NcType object"
        );
    }
    const NcType tmpType(getType(ncType.getName(), NcGroup::ParentsAndCurrent));
    if (tmpType.isNull()) {
        throw PFNC_CPP_NullType(
            "Attempt to invoke NcGroup::addVar failed: NcType must be defined in either the current group or a parent group"
        );
    }

    const auto dimIds = checkDimensionsAreValid(*this, ncDimVector);

    int varId;
    const int* dimIdsPtr = dimIds.empty() ? nullptr : dimIds.data();
    PFNC_CPP_CHECK(pfnc_def_var(
        myId,
        name.c_str(),
        tmpType.getId(),
        static_cast<int>(dimIds.size()),
        dimIdsPtr,
        &varId
    ));

    return NcVar(*this, varId);
}

NcVar NcGroup::addVar(
    const std::string& name,
    const NcAtomicType& ncType,
    const std::vector<NcDim>& ncDimVector
) const {
    ncCheckDefineMode(myId);
    const auto dimIds = checkDimensionsAreValid(*this, ncDimVector);

    int varId;
    const int* dimIdsPtr = dimIds.empty() ? nullptr : dimIds.data();
    PFNC_CPP_CHECK(pfnc_def_var(
        myId,
        name.c_str(),
        ncType.getId(),
        static_cast<int>(dimIds.size()),
        dimIdsPtr,
        &varId
    ));
    // return an NcVar object for this new variable
    return NcVar(*this, varId);
}

void NcVar::getVar(std::complex<double>* dataValues) const {
    PFNC_CPP_CHECK(
        pfnc_get_var_double_complex(groupId(), getId(), cpp_to_c_complex(dataValues))
    );
}
void NcVar::getVar(const std::vector<size_t>& index, std::complex<double>* dataValues)
    const {
    PFNC_CPP_CHECK(pfnc_get_var1_double_complex(
        groupId(), getId(), index.data(), cpp_to_c_complex(dataValues)
    ));
}
void NcVar::getVar(
    const std::vector<size_t>& startp,
    const std::vector<size_t>& countp,
    std::complex<double>* dataValues
) const {
    getVar(startp, countp, {}, dataValues);
}
void NcVar::getVar(
    const std::vector<size_t>& startp,
    const std::vector<size_t>& countp,
    const std::vector<ptrdiff_t>& stridep,
    std::complex<double>* dataValues
) const {
    PFNC_CPP_CHECK(pfnc_get_vars_double_complex(
        groupId(),
        getId(),
        startp.data(),
        countp.data(),
        stridep.empty() ? nullptr : stridep.data(),
        cpp_to_c_complex(dataValues)
    ));
}

void NcVar::getVar(std::complex<float>* dataValues) const {
    PFNC_CPP_CHECK(
        pfnc_get_var_float_complex(groupId(), getId(), cpp_to_c_complex(dataValues))
    );
}
void NcVar::getVar(const std::vector<size_t>& index, std::complex<float>* dataValues)
    const {
    PFNC_CPP_CHECK(pfnc_get_var1_float_complex(
        groupId(), getId(), index.data(), cpp_to_c_complex(dataValues)
    ));
}
void NcVar::getVar(
    const std::vector<size_t>& startp,
    const std::vector<size_t>& countp,
    std::complex<float>* dataValues
) const {
    getVar(startp, countp, {}, dataValues);
}
void NcVar::getVar(
    const std::vector<size_t>& startp,
    const std::vector<size_t>& countp,
    const std::vector<ptrdiff_t>& stridep,
    std::complex<float>* dataValues
) const {
    PFNC_CPP_CHECK(pfnc_get_vars_float_complex(
        groupId(),
        getId(),
        startp.data(),
        countp.data(),
        stridep.empty() ? nullptr : stridep.data(),
        cpp_to_c_complex(dataValues)
    ));
}

void NcVar::putVar(const std::complex<double>* dataValues) const {
    ncCheckDataMode(groupId());

    PFNC_CPP_CHECK(
        pfnc_put_var_double_complex(groupId(), getId(), cpp_to_c_complex(dataValues))
    );
}

void NcVar::putVar(
    const std::vector<size_t>& index, const std::complex<double>* dataValues
) const {
    ncCheckDataMode(groupId());
    PFNC_CPP_CHECK(pfnc_put_var1_double_complex(
        groupId(), getId(), index.data(), cpp_to_c_complex(dataValues)
    ));
}

void NcVar::putVar(
    const std::vector<size_t>& startp,
    const std::vector<size_t>& countp,
    const std::complex<double>* dataValues
) const {
    putVar(startp, countp, {}, dataValues);
}

void NcVar::putVar(
    const std::vector<size_t>& startp,
    const std::vector<size_t>& countp,
    const std::vector<ptrdiff_t>& stridep,
    const std::complex<double>* dataValues
) const {
    ncCheckDataMode(groupId());

    PFNC_CPP_CHECK(pfnc_put_vars_double_complex(
        groupId(),
        getId(),
        startp.data(),
        countp.data(),
        stridep.empty() ? nullptr : stridep.data(),
        cpp_to_c_complex(dataValues)
    ));
}

void NcVar::putVar(const std::complex<float>* dataValues) const {
    ncCheckDataMode(groupId());

    PFNC_CPP_CHECK(
        pfnc_put_var_float_complex(groupId(), getId(), cpp_to_c_complex(dataValues))
    );
}

bool NcVar::isComplex() const {
    return pfnc_var_is_complex(groupId(), getId());
}
bool NcVar::isComplexType() const {
    return pfnc_var_is_complex_type(groupId(), getId());
}
bool NcVar::hasComplexDim() const {
    return pfnc_var_has_complex_dimension(groupId(), getId());
}

int NcVar::getDimCount() const {
    int dim_count;
    PFNC_CPP_CHECK(pfnc_inq_varndims(groupId(), getId(), &dim_count));
    return dim_count;
}

std::vector<NcDim> NcVar::getDims() const {
    const auto dim_count = static_cast<std::size_t>(getDimCount());
    if (dim_count == 0) {
        return {};
    }

    std::vector<int> dim_ids(dim_count);
    PFNC_CPP_CHECK(pfnc_inq_vardimid(groupId(), getId(), dim_ids.data()));

    std::vector<NcDim> dims(dim_count);
    const auto parent = getParentGroup();
    std::transform(
        dim_ids.begin(),
        dim_ids.end(),
        dims.begin(),
        [&parent](auto dim_id) {
            return NcDim{parent, dim_id};
        }
    );
    return dims;
}

NcDim NcVar::getDim(int i) const {
    const auto dims = getDims();
    if (i < 0 || static_cast<std::size_t>(i) >= dims.size()) {
        throw exceptions::NcException("Index out of range", __FILE__, __LINE__);
    }
    return dims[static_cast<std::size_t>(i)];
}

NcFile::NcFile(const std::string& filePath, FileMode fMode, FileFormat fFormat) {
    open(filePath, fMode, fFormat);
}

NcFile::~NcFile() {
    // destructor may be called due to an exception being thrown
    // hence throwing an exception from within a destructor
    // causes undefined behaviour! so just printing a warning message
    try {
        close();
    } catch (const NcException& e) {
        std::cerr << e.what() << "\n";
    }
}

void NcFile::open(const std::string& filePath, FileMode fMode, FileFormat fFormat) {
    const int format = [fFormat]() {
        switch (fFormat) {
        case nc_complex::NcFile::FileFormat::classic:
            return 0;
        case nc_complex::NcFile::FileFormat::classic64:
            return NC_64BIT_OFFSET;
        case nc_complex::NcFile::FileFormat::nc4:
            return NC_NETCDF4;
        case nc_complex::NcFile::FileFormat::nc4classic:
            return NC_NETCDF4 | NC_CLASSIC_MODEL;
        default:
            throw exceptions::NcException(
                "Bad format in NcFile constructor", __FILE__, __LINE__
            );
        }
    }();

    switch (fMode) {
    case nc_complex::NcFile::FileMode::write:
        PFNC_CPP_CHECK(nc_open(filePath.c_str(), format | NC_WRITE, &myId));
        break;
    case nc_complex::NcFile::FileMode::read:
        PFNC_CPP_CHECK(nc_open(filePath.c_str(), format | NC_NOWRITE, &myId));
        break;
    case nc_complex::NcFile::FileMode::newFile:
        PFNC_CPP_CHECK(nc_create(filePath.c_str(), format | NC_NOCLOBBER, &myId));
        break;
    case nc_complex::NcFile::FileMode::replace:
        PFNC_CPP_CHECK(nc_create(filePath.c_str(), format | NC_CLOBBER, &myId));
        break;
    }

    nullObject = false;
}

void NcFile::close() {
    if (!nullObject) {
        PFNC_CPP_CHECK(nc_close(myId));
    }

    nullObject = true;
}

void NcFile::sync() {
    PFNC_CPP_CHECK(nc_sync(myId));
}
// Set fill mode for netCDF dataset open for writing and return current fill mode
void NcFile::set_Fill(int fillmode, int* old_modep) {
    PFNC_CPP_CHECK(nc_set_fill(myId, fillmode, old_modep));
}

// Put open netCDF dataset into define mode
void NcFile::redef() {
    PFNC_CPP_CHECK(nc_redef(myId));
}

// Leave define mode, used for classic model
void NcFile::enddef() {
    PFNC_CPP_CHECK(nc_enddef(myId));
}

}  // namespace nc_complex
