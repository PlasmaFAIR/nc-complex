from __future__ import annotations

import numpy as np

import posixpath
import netCDF4
from netCDF4.utils import _find_dim

from ._c_nc_complex import double_complex_typeid

from typing import Union


class Variable(netCDF4.Variable):
    def __init__(
        self,
        group: Dataset,
        name: str,
        datatype: Union[str, netCDF4.CompoundType],
        *args,
        **kwargs,
    ):
        # If datatype is complex dtype, then create a new Python object representing it
        is_complex = datatype == "c16"
        if is_complex:
            # Either get existing typeid or create new one
            datatype_id = double_complex_typeid(group._grpid)
            # This is just the Python wrapper around the datatype
            datatype = netCDF4.CompoundType(
                group, "f8, f8", "complex", typeid=datatype_id
            )
            # Now something ugly: we have to lie about the datatype here in
            # order to convince Variable it has a numpy complex dtype so it
            # casts things correctly -- we can't change it directly, because
            # Variable.__setattr__ is locked down. We also have to change it
            # _back_ later because some string-conversion function falls over
            # when it tries to iterate over the non-existent fields of
            # np.dtype("c16")
            old_datatype_dtype = datatype.dtype
            datatype.dtype = np.dtype("c16")

        super().__init__(group, name, datatype, *args, **kwargs)

        # Don't forget to restore the old dtype
        if is_complex:
            datatype.dtype = old_datatype_dtype

    def __getitem__(self, key):
        data = super().__getitem__(key)
        # handle converting to complex

        return data


class Dataset(netCDF4.Dataset):
    def createVariable(self, varname, datatype, dimensions=(), *args, **kwargs):
        dirname, varname = posixpath.split(posixpath.normpath(varname))

        if not dirname:
            group = self
        else:
            group = self.createGroup(dirname)

        if isinstance(dimensions, (str, bytes, netCDF4.Dimension)):
            dimensions = (dimensions,)

        dimensions = tuple(
            _find_dim(group, d) if isinstance(d, (str, bytes)) else d
            for d in dimensions
        )
        # create variable.
        group.variables[varname] = Variable(
            group, varname, datatype, dimensions, *args, **kwargs
        )
        return group.variables[varname]
