import numpy as np
from numpy import ma
from numpy.lib.stride_tricks import as_strided


def _safecast(a, b):
    """Check to see if array a can be safely cast to array b.  A
    little less picky than numpy.can_cast.
    """
    try:
        return ((a == b) | (np.isnan(a) & np.isnan(b))).all()
    except:
        try:
            return (a == b).all()  # string arrays.
        except:
            return False


def _sortbylist(A, B):
    """sort one list (A) using the values from another list (B)"""
    return [A[i] for i in sorted(range(len(A)), key=B.__getitem__)]


def _find_dim(grp, dimname):
    """Find Dimension instance given group and name.  look in current
    group, and parents.

    """
    group = grp
    dim = None
    while 1:
        try:
            dim = group.dimensions[dimname]
            break
        except KeyError:
            try:
                group = group.parent
            except AttributeError:
                raise ValueError(
                    f"cannot find dimension {dimname} in this group or parent groups"
                )
    if dim is None:
        raise KeyError(
            f"dimension {dimname} not defined in group {grp.path} or any group in its family tree"
        )

    return dim


def _walk_grps(topgrp):
    """Iterate through all (sub-) groups of topgrp, similar to os.walktree."""
    yield topgrp.groups.values()
    for grp in topgrp.groups.values():
        yield from _walk_grps(grp)


def _quantize(data, least_significant_digit):
    """
    quantize data to improve compression. data is quantized using
    around(scale*data)/scale, where scale is 2**bits, and bits is determined
    from the least_significant_digit. For example, if
    least_significant_digit=1, bits will be 4.
    """
    precision = pow(10.0, -least_significant_digit)
    exp = np.log10(precision)
    if exp < 0:
        exp = int(np.floor(exp))
    else:
        exp = int(np.ceil(exp))
    bits = np.ceil(np.log2(pow(10.0, -exp)))
    scale = pow(2.0, bits)
    datout = np.around(scale * data) / scale
    if ma.isMA(datout):
        datout.set_fill_value(data.fill_value)
        return datout
    else:
        return datout


def _StartCountStride(
    elem, shape, dimensions=None, grp=None, datashape=None, put=False, use_get_vars=True
):
    """Return start, count, stride and indices needed to store/extract data
    into/from a netCDF variable.

    This function is used to convert a slicing expression into a form that is
    compatible with the nc_get_vars function. Specifically, it needs
    to interpret integers, slices, Ellipses, and 1-d sequences of integers
    and booleans.

    Numpy uses "broadcasting indexing" to handle array-valued indices.
    "Broadcasting indexing" (a.k.a "fancy indexing") treats all multi-valued
    indices together to allow arbitrary points to be extracted. The index
    arrays can be multidimensional, and more than one can be specified in a
    slice, as long as they can be "broadcast" against each other.
    This style of indexing can be very powerful, but it is very hard
    to understand, explain, and implement (and can lead to hard to find bugs).
    Most other python packages and array processing
    languages (such as netcdf4-python, xray, biggus, matlab and fortran)
    use "orthogonal indexing" which only allows for 1-d index arrays and
    treats these arrays of indices independently along each dimension.

    The implementation of "orthogonal indexing" used here requires that
    index arrays be 1-d boolean or integer. If integer arrays are used,
    the index values must be sorted and contain no duplicates.

    In summary, slicing netcdf4-python variable objects with 1-d integer or
    boolean arrays is allowed, but may give a different result than slicing a
    numpy array.

    Numpy also supports slicing an array with a boolean array of the same
    shape. For example x[x>0] returns a 1-d array with all the positive values of x.
    This is also not supported in netcdf4-python, if x.ndim > 1.

    Orthogonal indexing can be used in to select netcdf variable slices
    using the dimension variables. For example, you can use v[lat>60,lon<180]
    to fetch the elements of v obeying conditions on latitude and longitude.
    Allow for this sort of simple variable subsetting is the reason we decided to
    deviate from numpy's slicing rules.

    This function is used both by the __setitem__ and __getitem__ method of
    the Variable class.

    Parameters
    ----------
    elem : tuple of integer, slice, ellipsis or 1-d boolean or integer
    sequences used to slice the netCDF Variable (Variable[elem]).
    shape : tuple containing the current shape of the netCDF variable.
    dimensions : sequence
      The name of the dimensions.
      __setitem__.
    grp  : netCDF Group
      The netCDF group to which the variable being set belongs to.
    datashape : sequence
      The shape of the data that is being stored. Only needed by __setitem__
    put : True|False (default False).  If called from __setitem__, put is True.

    Returns
    -------
    start : ndarray (..., n)
      A starting indices array of dimension n+1. The first n
      dimensions identify different independent data chunks. The last dimension
      can be read as the starting indices.
    count : ndarray (..., n)
      An array of dimension (n+1) storing the number of elements to get.
    stride : ndarray (..., n)
      An array of dimension (n+1) storing the steps between each datum.
    indices : ndarray (..., n)
      An array storing the indices describing the location of the
      data chunk in the target/source array (__getitem__/__setitem__).

    Notes:

    netCDF data is accessed via the function:
       nc_get_vars(grpid, varid, start, count, stride, data)

    Assume that the variable has dimension n, then

    start is a n-tuple that contains the indices at the beginning of data chunk.
    count is a n-tuple that contains the number of elements to be accessed.
    stride is a n-tuple that contains the step length between each element.

    """
    # Adapted from pycdf (http://pysclint.sourceforge.net/pycdf)
    # by Andre Gosselin..
    # Modified by David Huard to handle efficiently fancy indexing with
    # sequences of integers or booleans.

    nDims = len(shape)
    if nDims == 0:
        nDims = 1
        shape = (1,)

    # is there an unlimited dimension? (only defined for __setitem__)
    if put:
        hasunlim = False
        unlimd = {}
        if dimensions:
            for i in range(nDims):
                dimname = dimensions[i]
                # is this dimension unlimited?
                # look in current group, and parents for dim.
                dim = _find_dim(grp, dimname)
                unlimd[dimname] = dim.isunlimited()
                if unlimd[dimname]:
                    hasunlim = True
    else:
        hasunlim = False

    # When a single array or (non-tuple) sequence of integers is given
    # as a slice, assume it applies to the first dimension,
    # and use ellipsis for remaining dimensions.
    if np.iterable(elem):
        if type(elem) == np.ndarray or (
            type(elem) != tuple and np.array([_is_int(e) for e in elem]).all()
        ):
            elem = [elem]
            for n in range(len(elem) + 1, nDims + 1):
                elem.append(slice(None, None, None))
    else:  # Convert single index to sequence
        elem = [elem]

    # ensure there is at most 1 ellipse
    #  we cannot use elem.count(Ellipsis), as with fancy indexing would occur
    #  np.array() == Ellipsis which gives ValueError: The truth value of an
    #  array with more than one element is ambiguous. Use a.any() or a.all()
    if sum(1 for e in elem if e is Ellipsis) > 1:
        raise IndexError("At most one ellipsis allowed in a slicing expression")

    # replace boolean arrays with sequences of integers.
    newElem = []
    IndexErrorMsg = "only integers, slices (`:`), ellipsis (`...`), and 1-d integer or boolean arrays are valid indices"
    i = 0
    for e in elem:
        # string-like object try to cast to int
        # needs to be done first, since strings are iterable and
        # hard to distinguish from something castable to an iterable numpy array.
        if isinstance(e, (str, bytes)):
            try:
                e = int(e)
            except ValueError:
                raise IndexError(IndexErrorMsg)
        ea = np.asarray(e)
        # Raise error if multidimensional indexing is used.
        if ea.ndim > 1:
            raise IndexError("Index cannot be multidimensional")
        # set unlim to True if dimension is unlimited and put==True
        # (called from __setitem__)
        if hasunlim and put and dimensions:
            try:
                dimname = dimensions[i]
                unlim = unlimd[dimname]
            except IndexError:  # more slices than dimensions (issue 371)
                unlim = False
        else:
            unlim = False
        # convert boolean index to integer array.
        if np.iterable(ea) and ea.dtype.kind == "b":
            # check that boolean array not too long
            if not unlim and shape[i] != len(ea):
                raise IndexError(
                    "Boolean array must have the same shape as the data along this dimension"
                )
            ea = np.flatnonzero(ea)
        # an iterable (non-scalar) integer array.
        if np.iterable(ea) and ea.dtype.kind == "i":
            # convert negative indices in 1d array to positive ones.
            ea = np.where(ea < 0, ea + shape[i], ea)
            if np.any(ea < 0):
                raise IndexError("integer index out of range")
            # if unlim, let integer index be longer than current dimension
            # length.
            if ea.shape != (0,):
                elen = shape[i]
                if unlim:
                    elen = max(ea.max() + 1, elen)
                if ea.max() + 1 > elen:
                    raise IndexError("integer index exceeds dimension size")
            newElem.append(ea)
        # integer scalar
        elif ea.dtype.kind == "i":
            newElem.append(e)
        # slice or ellipsis object
        elif type(e) == slice or type(e) == type(Ellipsis):
            if (
                not use_get_vars
                and type(e) == slice
                and e.step not in [None, -1, 1]
                and dimensions is not None
                and grp is not None
            ):
                # convert strided slice to integer sequence if possible
                # (this will avoid nc_get_vars, which is slow - issue #680).
                start = e.start if e.start is not None else 0
                step = e.step
                if e.stop is None and dimensions is not None and grp is not None:
                    stop = len(_find_dim(grp, dimensions[i]))
                else:
                    stop = e.stop
                    if stop < 0:
                        stop = len(_find_dim(grp, dimensions[i])) + stop
                try:
                    ee = np.arange(start, stop, e.step)
                    if len(ee) > 0:
                        e = ee
                except:
                    pass
            newElem.append(e)
        else:  # castable to a scalar int, otherwise invalid
            try:
                e = int(e)
                newElem.append(e)
            except:
                raise IndexError(IndexErrorMsg)
        if type(e) == type(Ellipsis):
            i += 1 + nDims - len(elem)
        else:
            i += 1
    elem = newElem

    # replace Ellipsis and integer arrays with slice objects, if possible.
    newElem = []
    for e in elem:
        ea = np.asarray(e)
        # Replace ellipsis with slices.
        if type(e) == type(Ellipsis):
            # The ellipsis stands for the missing dimensions.
            newElem.extend((slice(None, None, None),) * (nDims - len(elem) + 1))
        # Replace sequence of indices with slice object if possible.
        elif np.iterable(e) and len(e) > 1:
            start = e[0]
            stop = e[-1] + 1
            step = e[1] - e[0]
            try:
                ee = range(start, stop, step)
            except ValueError:  # start, stop or step is not valid for a range
                ee = False
            if ee and len(e) == len(ee) and (e == np.arange(start, stop, step)).all():
                # don't convert to slice unless abs(stride) == 1
                # (nc_get_vars is very slow, issue #680)
                if not use_get_vars and step not in [1, -1]:
                    newElem.append(e)
                else:
                    newElem.append(slice(start, stop, step))
            else:
                newElem.append(e)
        elif np.iterable(e) and len(e) == 1:
            newElem.append(slice(e[0], e[0] + 1, 1))
        else:
            newElem.append(e)
    elem = newElem

    # If slice doesn't cover all dims, assume ellipsis for rest of dims.
    if len(elem) < nDims:
        for n in range(len(elem) + 1, nDims + 1):
            elem.append(slice(None, None, None))

    # make sure there are not too many dimensions in slice.
    if len(elem) > nDims:
        raise ValueError(
            "slicing expression exceeds the number of dimensions of the variable"
        )

    # Compute the dimensions of the start, count, stride and indices arrays.
    # The number of elements in the first n dimensions corresponds to the
    # number of times the _get method will be called.
    sdim = []
    for i, e in enumerate(elem):
        # at this stage e is a slice, a scalar integer, or a 1d integer array.
        # integer array:  _get call for each True value
        if np.iterable(e):
            sdim.append(len(e))
        # Scalar int or slice, just a single _get call
        else:
            sdim.append(1)

    # broadcast data shape when assigned to full variable (issue #919)
    try:
        fullslice = elem.count(slice(None, None, None)) == len(elem)
    except:  # fails if elem contains a numpy array.
        fullslice = False
    if fullslice and datashape and put and not hasunlim:
        datashape = broadcasted_shape(shape, datashape)

    # pad datashape with zeros for dimensions not being sliced (issue #906)
    # only used when data covers slice over subset of dimensions
    if (
        datashape
        and len(datashape) != len(elem)
        and len(datashape) == sum(1 for e in elem if type(e) == slice)
    ):
        datashapenew = ()
        i = 0
        for e in elem:
            if type(e) != slice and not np.iterable(e):  # scalar integer slice
                datashapenew = datashapenew + (0,)
            else:  # slice object
                datashapenew = datashapenew + (datashape[i],)
                i += 1
        datashape = datashapenew

    # Create the start, count, stride and indices arrays.

    sdim.append(max(nDims, 1))
    start = np.empty(sdim, dtype=np.intp)
    count = np.empty(sdim, dtype=np.intp)
    stride = np.empty(sdim, dtype=np.intp)
    indices = np.empty(sdim, dtype=object)

    for i, e in enumerate(elem):
        ea = np.asarray(e)

        # set unlim to True if dimension is unlimited and put==True
        # (called from __setitem__). Note: grp and dimensions must be set.
        if hasunlim and put and dimensions:
            dimname = dimensions[i]
            unlim = unlimd[dimname]
        else:
            unlim = False

        #    SLICE    #
        if type(e) == slice:
            # determine length parameter for slice.indices.

            # shape[i] can be zero for unlim dim that hasn't been written to
            # yet.
            # length of slice may be longer than current shape
            # if dimension is unlimited (and we are writing, not reading).
            if unlim and e.stop is not None and e.stop > shape[i]:
                length = e.stop
            elif unlim and e.stop is None and datashape != ():
                try:
                    if e.start is None:
                        length = datashape[i]
                    else:
                        length = e.start + datashape[i]
                except IndexError:
                    raise IndexError("shape of data does not conform to slice")
            else:
                if unlim and datashape == () and len(dim) == 0:
                    # writing scalar along unlimited dimension using slicing
                    # syntax (var[:] = 1, when var.shape = ())
                    length = 1
                else:
                    length = shape[i]

            beg, end, inc = e.indices(length)
            n = len(range(beg, end, inc))

            start[..., i] = beg
            count[..., i] = n
            stride[..., i] = inc
            indices[..., i] = slice(None)

        #    ITERABLE    #
        elif np.iterable(e) and np.array(e).dtype.kind in "i":  # Sequence of integers
            start[..., i] = np.apply_along_axis(lambda x: e * x, i, np.ones(sdim[:-1]))
            indices[..., i] = np.apply_along_axis(
                lambda x: np.arange(sdim[i]) * x, i, np.ones(sdim[:-1], int)
            )

            count[..., i] = 1
            stride[..., i] = 1

        #   all that's left is SCALAR INTEGER    #
        else:
            if e >= 0:
                start[..., i] = e
            elif e < 0 and (-e <= shape[i]):
                start[..., i] = e + shape[i]
            else:
                raise IndexError("Index out of range")

            count[..., i] = 1
            stride[..., i] = 1
            indices[..., i] = -1  # Use -1 instead of 0 to indicate that
            # this dimension shall be squeezed.

    return start, count, stride, indices  # , out_shape


def _out_array_shape(count):
    """Return the output array shape given the count array created by getStartCountStride"""

    s = list(count.shape[:-1])
    out = []

    for i, n in enumerate(s):
        if n == 1 and count.size > 0:
            c = count[..., i].ravel()[0]  # All elements should be identical.
            out.append(c)
        else:
            out.append(n)
    return out


def _is_int(a):
    try:
        return int(a) == a
    except ValueError:
        return False


def _tostr(s):
    try:
        return str(s)
    except ValueError:
        return s


def broadcasted_shape(shp1, shp2):
    # determine shape of array of shp1 and shp2 broadcast against one another.
    x = np.array([1])
    # trick to define array with certain shape that doesn't allocate all the
    # memory.
    a = as_strided(x, shape=shp1, strides=[0] * len(shp1))
    b = as_strided(x, shape=shp2, strides=[0] * len(shp2))
    return np.broadcast(a, b).shape
