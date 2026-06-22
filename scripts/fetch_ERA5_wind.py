#!/usr/bin/env python3
"""
Fetch long-term mean wind profile at a given location from ERA5 (CDS API).

Usage
-----
    python fetch_ERA5_wind.py <lat> <lon> <cache_dir>

Output
------
    JSON on stdout:  {"heights": [10, 100], "v_mean": [<ws10>, <ws100>]}
    The downloaded NetCDF is cached in <cache_dir> so subsequent calls
    with the same lat/lon are instantaneous.

Requirements
------------
    pip install cdsapi netCDF4 numpy

CDS API key
-----------
    Register at https://cds.climate.copernicus.eu, then create ~/.cdsapirc:
        url: https://cds.climate.copernicus.eu/api/v2
        key: <UID>:<API-KEY>
"""

import sys
import os
import json
import numpy as np

try:
    import cdsapi
except ImportError:
    sys.exit('ERROR: cdsapi not installed. Run: pip install cdsapi')

try:
    import netCDF4 as nc
except ImportError:
    sys.exit('ERROR: netCDF4 not installed. Run: pip install netCDF4')


def fetch(lat, lon, cache_dir):
    os.makedirs(cache_dir, exist_ok=True)
    nc_file = os.path.join(cache_dir, f'era5_{lat:.3f}_{lon:.3f}.nc')

    if not os.path.exists(nc_file):
        # ERA5 single-levels monthly means, 2000-2023
        # Variables: 10m and 100m wind components (u, v)
        c = cdsapi.Client()  # reads ~/.cdsapirc automatically
        c.retrieve(
            'reanalysis-era5-single-levels-monthly-means',
            {
                'product_type': ['monthly_averaged_reanalysis'],
                'variable': [
                    '10m_u_component_of_wind',
                    '10m_v_component_of_wind',
                    '100m_u_component_of_wind',
                    '100m_v_component_of_wind',
                ],
                'year':  [str(y) for y in range(2000, 2024)],
                'month': [f'{m:02d}' for m in range(1, 13)],
                'time':  '00:00',
                # Bounding box: ±0.25° around the point (one ERA5 grid cell)
                'area':  [lat + 0.25, lon - 0.25, lat - 0.25, lon + 0.25],
                'data_format': 'netcdf',
                'download_format': 'unarchived',
            },
        ).download(nc_file)

    # Compute long-term mean wind speed (scalar average of monthly means)
    ds = nc.Dataset(nc_file)
    ws10  = float(np.mean(np.sqrt(ds['u10'][:] ** 2 + ds['v10'][:] ** 2)))
    ws100 = float(np.mean(np.sqrt(ds['u100'][:] ** 2 + ds['v100'][:] ** 2)))
    ds.close()

    return {'heights': [10, 100], 'v_mean': [ws10, ws100]}


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit('Usage: fetch_ERA5_wind.py <lat> <lon> [<cache_dir>]')

    lat_in   = float(sys.argv[1])
    lon_in   = float(sys.argv[2])
    cache_in = sys.argv[3] if len(sys.argv) > 3 else '.'

    result = fetch(lat_in, lon_in, cache_in)
    print(json.dumps(result))
