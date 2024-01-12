#include <complex>
#include <cstddef>
#include <netcdf>
#include <vector>

#include "nc_complex.h"
#include "nc_complex/nc_complex_export.h"

namespace nc_complex {

// Bring in everything we aren't wrapping from the netCDF namespace into our
// one. Then all names can be used from our namespace
using netCDF::NcAtomicType;
using netCDF::NcCompoundType;
using netCDF::NcOpaqueType;
using netCDF::NcType;
using netCDF::NcVlenType;

using netCDF::ncByte;
using netCDF::ncChar;
using netCDF::ncDouble;
using netCDF::ncFloat;
using netCDF::ncInt;
using netCDF::ncInt64;
using netCDF::ncShort;
using netCDF::ncString;
using netCDF::ncUbyte;
using netCDF::ncUint;
using netCDF::ncUint64;
using netCDF::ncUshort;

using netCDF::NcDim;

using netCDF::NcAtt;
using netCDF::NcGroupAtt;
using netCDF::NcVarAtt;

using netCDF::ncCheck;
using netCDF::ncCheckDataMode;
using netCDF::ncCheckDefineMode;

namespace exceptions {
using namespace netCDF::exceptions;
}

class NcVar;

class NC_COMPLEX_EXPORT NcGroup : public netCDF::NcGroup {
public:
    using netCDF::NcGroup::NcGroup;

    // Convert base class to specialised class. We need this so we can
    // wrap methods on `netCDF::NcGroup` and return this class instead
    // of `netCDF::NcVar`.
    explicit NcGroup(const netCDF::NcGroup& other) : netCDF::NcGroup(other) {}

    NcGroup getParentGroup() const {
        return NcGroup{netCDF::NcGroup::getParentGroup()};
    }
    std::multimap<std::string, NcGroup> getGroups(
        NcGroup::GroupLocation location = ChildrenGrps
    ) const;
    std::set<NcGroup> getGroups(
        const std::string& name, NcGroup::GroupLocation location = ChildrenGrps
    ) const;
    NcGroup getGroup(
        const std::string& name, NcGroup::GroupLocation location = ChildrenGrps
    ) const;
    NcGroup addGroup(const std::string& name) const;

    std::multimap<std::string, NcVar> getVars(NcGroup::Location location = Current)
        const;
    std::set<NcVar> getVars(
        const std::string& name, NcGroup::Location location = Current
    ) const;
    NcVar getVar(const std::string& name, NcGroup::Location location = Current) const;

    NcVar addVar(const std::string& name, const NcType& ncType) const;
    NcVar addVar(const std::string& name, const NcAtomicType& ncType) const;

    NcVar addVar(
        const std::string& name, const std::string& typeName, const std::string& dimName
    ) const;
    NcVar addVar(const std::string& name, const NcType& ncType, const NcDim& ncDim)
        const;
    NcVar addVar(
        const std::string& name, const NcAtomicType& ncType, const NcDim& ncDim
    ) const;

    NcVar addVar(
        const std::string& name,
        const std::string& typeName,
        const std::vector<std::string>& dimNames
    ) const;
    NcVar addVar(
        const std::string& name,
        const NcType& ncType,
        const std::vector<NcDim>& ncDimVector
    ) const;
    NcVar addVar(
        const std::string& name,
        const NcAtomicType& ncType,
        const std::vector<NcDim>& ncDimVector
    ) const;
};

// We have to completely reimplement `NcFile` so that it inherits from *our*
// `NcGroup`

class NC_COMPLEX_EXPORT NcFile : public NcGroup {
public:
    enum FileMode {
        read,  //!< File exists, open read-only.
        write,  //!< File exists, open for writing.
        replace,  //!< Create new file, even if already exists.
        newFile  //!< Create new file, fail if already exists.
    };

    enum FileFormat {
        classic,  //!< Classic format, classic data model
        classic64,  //!< 64-bit offset format, classic data model
        nc4,  //!< (default) netCDF-4/HDF5 format, enhanced data model
        nc4classic  //!< netCDF-4/HDF5 format, classic data model
    };

    NcFile() = default;
    NcFile(
        const std::string& filePath,
        FileMode fMode,
        FileFormat fFormat = FileFormat::nc4
    );

    ~NcFile() override;

    /// Do not allow definition of NcFile involving copying any NcFile or NcGroup.
    /// Because the destructor closes the file and releases al resources such
    /// an action could leave NcFile objects in an invalid state
    NcFile& operator=(const NcGroup& rhs) = delete;
    NcFile& operator=(const NcFile& rhs) = delete;
    NcFile(const NcGroup& rhs) = delete;
    NcFile(const NcFile& rhs) = delete;

    NcFile& operator=(NcFile&& rhs) = delete;
    NcFile(NcFile&& rhs) = delete;

    void open(
        const std::string& filePath,
        FileMode fMode,
        FileFormat fFormat = FileFormat::nc4
    );
    void close();

    void sync();
    void set_Fill(int fillmode, int* old_modep);
    void redef();
    void enddef();
};

class NC_COMPLEX_EXPORT NcFloatComplex : public NcAtomicType {
public:
    NcFloatComplex() : NcAtomicType(PFNC_FLOAT_COMPLEX) {}
};

class NC_COMPLEX_EXPORT NcDoubleComplex : public NcAtomicType {
public:
    NcDoubleComplex() : NcAtomicType(PFNC_DOUBLE_COMPLEX) {}
};

class NC_COMPLEX_EXPORT NcFloatComplexDim : public NcAtomicType {
public:
    NcFloatComplexDim() : NcAtomicType(PFNC_FLOAT_COMPLEX_DIM) {}
};

class NC_COMPLEX_EXPORT NcDoubleComplexDim : public NcAtomicType {
public:
    NcDoubleComplexDim() : NcAtomicType(PFNC_DOUBLE_COMPLEX_DIM) {}
};

class NC_COMPLEX_EXPORT NcVar : public netCDF::NcVar {
public:
    using netCDF::NcVar::NcVar;

    // Convert base class to specialised class. We need this so we can
    // wrap methods on `netCDF::NcGroup` and return this class instead
    // of `netCDF::NcVar`.
    explicit NcVar(const netCDF::NcVar& other) : netCDF::NcVar(other) {}

    using netCDF::NcVar::getVar;
    using netCDF::NcVar::putVar;

    void getVar(std::complex<double>* dataValues) const;
    void getVar(const std::vector<size_t>& index, std::complex<double>* dataValues)
        const;
    void getVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        std::complex<double>* dataValues
    ) const;
    void getVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        const std::vector<ptrdiff_t>& stridep,
        std::complex<double>* dataValues
    ) const;

    void getVar(std::complex<float>* dataValues) const;
    void getVar(const std::vector<size_t>& index, std::complex<float>* dataValues)
        const;
    void getVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        std::complex<float>* dataValues
    ) const;
    void getVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        const std::vector<ptrdiff_t>& stridep,
        std::complex<float>* dataValues
    ) const;

    void putVar(const std::complex<double>* dataValues) const;
    void putVar(
        const std::vector<size_t>& index, const std::complex<double>* dataValues
    ) const;
    void putVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        const std::complex<double>* dataValues
    ) const;
    void putVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        const std::vector<ptrdiff_t>& stridep,
        const std::complex<double>* dataValues
    ) const;

    void putVar(const std::complex<float>* dataValues) const;
    void putVar(const std::vector<size_t>& index, const std::complex<float>* dataValues)
        const;
    void putVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        const std::complex<float>* dataValues
    ) const;
    void putVar(
        const std::vector<size_t>& startp,
        const std::vector<size_t>& countp,
        const std::vector<ptrdiff_t>& stridep,
        const std::complex<float>* dataValues
    ) const;

    /// Returns true if the variable is complex
    bool isComplex() const;
    /// Returns true if the variable is complex and uses a compound datatype
    bool isComplexType() const;
    /// Returns true if the variable is complex and uses a complex dimension
    bool hasComplexDim() const;

    /*! The the number of dimensions. */
    int getDimCount() const;

    /*! Gets the i'th NcDim object. */
    netCDF::NcDim getDim(int i) const;

    /*! Gets the set of NcDim objects. */
    std::vector<netCDF::NcDim> getDims() const;

private:
    int groupId() const { return this->getParentGroup().getId(); }
};

}  // namespace nc_complex
