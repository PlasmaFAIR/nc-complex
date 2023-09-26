import nc_complex

import numpy as np

complex_array = np.array([0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j, 0.25 + 0.75j], dtype="c16")


def test_read(tmp_path):
    filename = str(tmp_path / "test.nc")

    with nc_complex.Dataset(filename, "w") as f:
        f.createDimension("x", size=len(complex_array))
        f.createDimension("ri", size=2)
        c_ri = f.createVariable("data_dim", np.float64, ("x", "ri"))
        as_dim_array = np.vstack((complex_array.real, complex_array.imag)).T
        c_ri[:] = as_dim_array

        np_dt = np.dtype("f8, f8")
        nc_dt = f.createCompoundType(np_dt, "nc_complex")
        c_struct = f.createVariable("data_struct", nc_dt, ("x",))
        as_struct_array = np.array(
            [(r, i) for r, i in zip(complex_array.real, complex_array.imag)],
            dtype="f8, f8",
        )
        c_struct[:] = as_struct_array


def test_write(tmp_path):
    filename = str(tmp_path / "test.nc")
    with nc_complex.Dataset(filename, "w") as f:
        f.createDimension("x", size=len(complex_array))
        var = f.createVariable("data", "c16", ("x",))
        var[:] = complex_array

    with nc_complex.Dataset(filename, "r") as f:
        assert "data" in f.variables
        assert np.array_equal(var, complex_array)


if __name__ == "__main__":
    import tempfile
    import pathlib

    with tempfile.TemporaryDirectory() as tmp_path:
        test_write(pathlib.Path(tmp_path))
