from __future__ import annotations

import numpy as np

import posixpath
from .utils import _find_dim
from ._pf_netCDF4 import (
    Variable as BaseVariable,
    Dataset as BaseDataset,
    CompoundType,
    Dimension,
)
from ._c_nc_complex import double_complex_typeid, is_complex, has_complex_dimension

from typing import Union, cast


def dtype_is_complex(dtype: Union[str, CompoundType, np.dtype]) -> bool:
    if dtype == "c16":
        return True

    return False


class Variable(BaseVariable):
    def __init__(
        self,
        group: Dataset,
        name: str,
        datatype: Union[str, CompoundType],
        *args,
        **kwargs,
    ):
        # If datatype is complex dtype, then create a new Python object representing it
        datatype_is_complex = dtype_is_complex(datatype)
        if datatype_is_complex:
            # Either get existing typeid or create new one
            datatype_id = double_complex_typeid(group._grpid)
            # This is just the Python wrapper around the datatype
            datatype = CompoundType(group, "f8, f8", "complex", typeid=datatype_id)
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
        if datatype_is_complex:
            datatype = cast(CompoundType, datatype)
            datatype.dtype = old_datatype_dtype

    def __getitem__(self, key):
        data = super().__getitem__(key)
        # handle converting to complex

        return data

    @property
    def is_complex(self):
        return is_complex(self.group()._grpid, self._varid)


class Dataset(BaseDataset):
    def __init__(self, filename, *args, **kwargs):
        super().__init__(filename, *args, **kwargs)

        # Recreate variables using our class
        for name, var in self.variables.items():
            var_is_complex = is_complex(self._grpid, var._varid)

            datatype = "c16" if var_is_complex else var.dtype
            dims = (
                var.dimensions
                if has_complex_dimension(self._grpid, var._varid)
                else var.dimensions[:-1]
            )

            dimensions = (_find_dim(self, dim) for dim in dims)

            self.variables[name] = Variable(
                self,
                var.name,
                datatype,
                dimensions,
                id=var._varid,
                endian=var.endian(),
            )

    def createVariable(self, varname, datatype, dimensions=(), *args, **kwargs):
        dirname, varname = posixpath.split(posixpath.normpath(varname))

        if not dirname:
            group = self
        else:
            group = self.createGroup(dirname)

        if isinstance(dimensions, (str, bytes, Dimension)):
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
