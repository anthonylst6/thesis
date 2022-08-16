# importing the module
import tracemalloc
 
# code or function for which memory
# has to be monitored
# Import libraries

from datetime import datetime
import glob
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import xskillscore as xs
import hvplot.xarray
import metpy.calc as mpcalc
import metpy.plots as mpplots
from matplotlib.patheffects import withStroke
from metpy.io import parse_metar_file
from metpy.units import pandas_dataframe_to_unit_arrays
from siphon.catalog import TDSCatalog
import re
import dask
import gc
from dask.distributed import Client

# starting the monitoring
tracemalloc.start()

ds_era5 = xr.open_dataset("../data_archived/sa_era5-slv_month-hour_1980-2021_LARGE.nc",
                          chunks = {'time': '500MB'},
                          engine = "netcdf4")

ds_era5["msl"].isel(longitude = slice(0, 5), latitude = slice(0, 5)).groupby("time").mean().compute()

ds_era5["msl"].isel(longitude = slice(0, 20), latitude = slice(0, 20)).groupby("time").mean().compute()

ds_era5["msl"].groupby("time").mean().compute()

# displaying the memory
print(tracemalloc.get_traced_memory())
 
# stopping the library
tracemalloc.stop()
