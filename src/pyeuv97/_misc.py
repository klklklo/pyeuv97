import functools
import xarray as xr
from importlib_resources import files


@functools.cache
def _read_coeffs(file):
    return xr.open_dataset(files('pyeuv97._coeffs').joinpath(file))


def get_euv97_coeffs():
    return _read_coeffs('prepared_intervals.nc').copy(), _read_coeffs('prepared_lines.nc')