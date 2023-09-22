import posixpath
import netCDF4
from netCDF4.utils import _find_dim

from ._c_nc_complex import double_complex_typeid


class Variable(netCDF4.Variable):
    def __init__(self, group, name, datatype, *args, **kwargs):
        # handle datatype is complex
        if datatype == "c16":
            datatype = double_complex_typeid(group._grpid)

        super().__init__(group, name, datatype, *args, **kwargs)

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
