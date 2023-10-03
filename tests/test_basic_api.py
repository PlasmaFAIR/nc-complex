import nc_complex
import netCDF4

import numpy as np

complex_array = np.array([0 + 0j, 1 + 0j, 0 + 1j, 1 + 1j, 0.25 + 0.75j], dtype="c16")


def test_read_dim(tmp_path):
    filename = str(tmp_path / "test.nc")

    with netCDF4.Dataset(filename, "w") as f:
        f.createDimension("x", size=len(complex_array))
        f.createDimension("ri", size=2)
        c_ri = f.createVariable("data_dim", np.float64, ("x", "ri"))
        as_dim_array = np.vstack((complex_array.real, complex_array.imag)).T
        c_ri[:] = as_dim_array

    with nc_complex.Dataset(filename, "r") as f:
        assert "data_dim" in f.variables
        assert f["data_dim"].is_complex
        breakpoint()
        data = f["data_dim"][:]

    assert np.array_equal(data, complex_array)


def test_read_struct(tmp_path):
    filename = str(tmp_path / "test.nc")

    with netCDF4.Dataset(filename, "w") as f:
        np_dt = np.dtype([("r", np.float64), ("i", np.float64)])
        f.createDimension("x", size=len(complex_array))
        nc_dt = f.createCompoundType(np_dt, "nc_complex")
        c_struct = f.createVariable("data_struct", nc_dt, ("x",))
        as_struct_array = np.array(
            [(r, i) for r, i in zip(complex_array.real, complex_array.imag)],
            dtype=np_dt,
        )
        c_struct[:] = as_struct_array

    with nc_complex.Dataset(filename, "r") as f:
        assert "data_struct" in f.variables
        assert f["data_struct"].is_complex
        data = f["data_struct"][:]

    assert np.array_equal(data, complex_array)


def test_write(tmp_path):
    filename = str(tmp_path / "test.nc")
    with nc_complex.Dataset(filename, "w") as f:
        f.createDimension("x", size=len(complex_array))
        complex_var = f.createVariable("complex_data", "c16", ("x",))
        complex_var[:] = complex_array
        assert complex_var.is_complex

        real_var = f.createVariable("real_data", "f8", ("x",))
        real_var[:] = np.arange(len(complex_array))

    with nc_complex.Dataset(filename, "r") as f:
        assert "complex_data" in f.variables
        assert f["complex_data"].is_complex
        assert not f["real_data"].is_complex
        assert np.array_equal(f["complex_data"], complex_array)


if __name__ == "__main__":
    import tempfile
    import pathlib

    with tempfile.TemporaryDirectory() as tmp_path:
        test_read_dim(pathlib.Path(tmp_path))
