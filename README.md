nc-complex
==========

`nc-complex` is a library for reading and writing complex numbers
using netCDF.

The problem
-----------

Complex numbers are widely used in the physical sciences (`[citation
needed]`). As of 2023, netCDF, a very popular file format for research
software, is lacking a native complex number datatype. This means that
software authors wishing to write complex data to netCDF have to roll
their own, with the result that there are several competing methods of
doing so.

One reason why netCDF is lacking a complex datatype is because its
main, modern backing file format, HDF5, is also lacking one (although
note that another modern backing format, Zarr, does). This has been
requested since at least 2010, and doesn't appear to be high on the
maintainers' priority list. Oh well.

The aim of `nc-complex` is to smooth over these differences and
present a single interface for reading and modifying complex data
stored in all common conventions, and writing to a single blessed
representation.

Conventions for Complex Numbers
-------------------------------

There are several existing conventions for storing complex numbers,
which boil down to three main flavours:

1. A compound datatype with fields for the real and imaginary
   components;
2. An extra dimension with length 2;
3. Separate variables for the real and imaginary components.

Each flavour has several variations, mainly around the names used for
the real/imaginary components, and whether the numbers are stored in
Cartesian/rectangular or polar convention (using magnitude and phase
angle).

You can find examples of all three flavours (and their variations) in
the wild, although given the nature of research software, it's
difficult to accurately survey the field and determine which flavour
is most used.

However, we can look at how complex numbers are natively represented
in programming languages with official netCDF bindings, as well as
HDF5 libraries, and other file formats commonly used in research
software.

Native Complex Numbers in Programming Languages
-----------------------------------------------

Of the various official language bindings for netCDF, four have native
complex number types: C, C++, Fortran, Python. There are also
unofficial bindings for other languages, which also have native
complex numbers, such as Julia, Matlab, and R.

In all of these languages, complex numbers have the same byte
representation which is equivalent to an array or struct of two
floating point numbers, the real part followed by the imaginary
part. Note that this doesn't seem to always be defined as such at the
language standard level, but does appear to always be true in
practice. All languages do explicitly use the Cartesian/rectangular
form, rather than polar, although they may have conversion functions.

[Cppreference][cpp_memcpy_example] demonstrates this layout:

```C
float a[4] = {1, 2, 3, 4};
float complex z1, z2;
memcpy(&z1, a, sizeof z1); // z1 becomes 1.0 + 2.0i
memcpy(&z2, a+2, sizeof z2); // z2 becomes 3.0 + 4.0i
```

The actual implementation in terms of an array, struct, or native type
varies between languages. For example, C++ defines a `std::complex`
type in the `<complex>` header, while Fortran has a native type `complex`.

Regardless of the implementation, most languages call the two
components `real` and `imag`. These may be struct members (or member
functions), as in C++ and Fortran, or free functions, as in Julia or C
(although in C they are called `creal` and `cimag`).

There is a lot more variation in the name of the imaginary unit (`i`
literal suffix in C and C++, `im` in Julia, `j` in Python), but this
is not important for us here.

Complex Numbers in other IO Libraries
-------------------------------------

We've looked at how complex numbers are represented in memory, now
let's consider how other IO libraries store them on disk. We're only
interested in libraries that aim for some degree of portability, and
particularly those used in the scientific community.

Let's first look at the Python and Julia bindings for HDF5:
[h5py][h5py] and [HDF5.jl][hdf5jl]. These both store complex numbers
using a compound datatype (flavour 1 above), with configurable field
names that default to `r`, `i`. HDF5.jl is explicitly modelled on
h5py. There's also [hdf5-rust][hdf5-rust] which also follows this
convention.

Matlab's `.mat` files are just plain HDF5 files, and since
version 2006b `[citation needed]` they also use compound types,
although with the field names `real` and `imag`.

[H5netcdf][h5netcdf] is "a Python implementation of the netCDF4 file
format base on h5py". Writing complex numbers uses h5py's datatype,
although an explicit `invalid_netcdf=True` flag must be passed when
opening a file in order to do so. The reason for this is that h5py
doesn't commit the datatype to the HDF5 file, instead it's stored in
the variable directly resulting in an unnamed (or "transient")
datatype. Previous versions of netCDF were unable to read variables
using such datatypes, although this was recently fixed due to work by
PlasmaFAIR, and the next version of netCDF will be able to read files
with complex numbers written by h5netcdf/h5py.

[Zarr][zarr] is "a format for the storage of chunked, compressed,
N-dimensional arrays", and can store files on disk, in a Zip file, or
in cloud storage. It's also a possible backing store for netCDF. Zarr
uses the [numpy dtypes][dtypes] system and can natively handle complex
numbers.

[ADIOS2][adios2] is a "unified high-performance IO framework" for the
HPC community, and can handle things like coupling code or streaming
data across a network, as well as writing to file. ADIOS2 supports
several backend engines, including HDF5, as well as its own
binary-packed `.bp` format. For most of its backends, complex numbers
seem to be streamed as-is, but for HDF5 they use a compound type with
fields `r` and `i`.

[TREXIO][trexio] is a file format for quantum chemistry calculations,
and is backed either by plain text or HDF5. This uses the split
variable scheme, with the real part having an undecorated name, and
the imaginary part an `"_im"` suffix.

Conventions used in applications
--------------------------------

It's difficult to accurately assess exactly what approaches are used
by the community. The following survey represents a couple of hours of
effort searching on GitHub for variations on "complex number" and
"netCDF" or "HDF5". Searching for the specific netCDF functions is
difficult as it mostly returns either forks or netCDF vendored/bundled
into another package.

Of course, this only covers open source software, as it's even more
difficult to find proprietary software that uses netCDF or HDF5 and
also supports complex numbers, let alone work out what convention it
might use.

### 1. Compound types

- [GDAL](https://gdal.org) (`"r"`, `"i"`)
    - Note: big ecosystem built on this!
    - Interface over many, many file types, but this applies to netCDF
      and HDF5
- [QCoDeS](https://github.com/QCoDeS/Qcodes) (`"r"`, `"i"`)
    - Uses h5netcdf
- [deal.ii](https://github.com/dealii/dealii) (`"r"`, `"i"`)
    - Actually writes to HDF5, netCDF output deprecated
- [DCA++](https://github.com/CompFUSE/DCA) (`"r"`, `"i"`)
    - HDF5
- [Armadillo](https://gitlab.com/conradsnicta/armadillo-code)
  (`"real"`, `"imag"`)
    - HDF5
- [FlexiBLAS](https://github.com/mpimd-csc/flexiblas) (`"real"`,
  `"imag"`)
    - HDF5
- [Eigen-HDF5](https://github.com/garrison/eigen3-hdf5) (`"r"`, `"i"`)
    - HDF5
- [Yardl](https://github.com/microsoft/yardl/) (`"real"`,
  `"imaginary"`)
    - HDF5 and JSON
- [pyuvdata](https://github.com/RadioAstronomySoftwareGroup/pyuvdata/)
  (`"r"`, `"i"`)
    - HDF5

### 2. Dimensions

- [GS2](https://bitbucket.org/gyrokinetics/gs2/) (`"ri"`)
- [Genray](https://github.com/compxco/genray) (`"two"`)
- [capytaine](https://github.com/capytaine/capytaine) (`"complex"`, with
  values `"re"`, `"im"`)
- [swot-hydrology-toolbox](https://github.com/CNES/swot-hydrology-toolbox)
  (`"depth"`)

### 3. Split variables

- [Nansat](https://github.com/nansencenter/nansat), (`"*_real"`, `"*_imag"`)
- [mafredo](https://github.com/RubendeBruin/mafredo) (`"real"`, `"imag"`
  variables in separate group per variable)
- [SONAR-netcdf4](https://github.com/ices-publications/SONAR-netCDF4) (`"*_r"`,
  `"*_i"`)
    - This is actually a schema itself, rather than a single application
- [cgenie.muffin](https://github.com/derpycode/cgenie.muffin) (`"*_re"`,
  `"*_im"`)
- [Abinit](https://github.com/abinit/abinit) (`"*_real"`, `"*_imag"`)
- [kineticj](https://github.com/ORNL-Fusion/kineticj) (`"*_re"`, `"*_im"`)
- [EDGI_PCA](https://github.com/cmdupuis3/EDGI_PCA) (`"*_re"`, `"*_im"`)
- [ITensor](https://github.com/ITensor/ITensor) (`"r"`, `"i"`
  variables in separate group per variable)

### 4. Other representations

- [STELLOPT](https://github.com/PrincetonUniversity/STELLOPT) using
  [EZcdf](https://w3.pppl.gov/ntcc/rib_repositories_NTCC_catalog_Asset/EZcdf.html),
  which doubles the length of the leading dimension

### Observations

Overall, I found 22 pieces of software (not including the explicit IO
libraries and wrappers) using 13 different representations. Packages
using HDF5 directly seem to exclusively use compound datatypes, and
almost all name the fields `r` and `i`. I only found four packages
using dimensions, and each one used a different name for it (including
the surprising `"depth"`!). Eight packages used two separate
variables, with five different naming conventions, including two which
use a group containing two variables for the separate components to
represent the whole variable. A single package used an apparently
unique representation, doubling the length of the fastest dimension.

Perhaps an important point to note, I struggled to find _any_ uses of
netCDF compound types in production use. Searches for
`nf90_def_compound` (and ruling out any forks or bundled copies of
netCDF-Fortran) turn up no results at all, while I found it almost
impossible to find any real results for `nc_def_compound`, wading
through the many bundled copies of netCDF-C. It's not immediately
obvious why this is case, especially compared to uses of compound
types in HDF5 code.

[cpp_memcpy_example]: https://en.cppreference.com/w/c/language/arithmetic_types#Complex_floating_types
[h5py]: https://docs.h5py.org/en/stable/index.html
[hdf5jl]: https://juliaio.github.io/HDF5.jl/stable/
[h5netcdf]: https://h5netcdf.org
[zarr]: https://zarr.readthedocs.io/en/stable/index.html
[dtypes]: https://numpy.org/doc/stable/reference/arrays.dtypes.html
[hdf5-rust]: https://github.com/aldanor/hdf5-rust
[adios2]: https://github.com/ornladios/ADIOS2/tree/c503940b06020fcbc1731424cca42f7f76dd4512
