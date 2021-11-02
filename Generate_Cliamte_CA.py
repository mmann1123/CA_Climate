#%%
import xclim
from xclim.testing import open_dataset
import xarray as xr

#%%

ds = open_dataset("ERA5/daily_surface_cancities_1990-1993.nc")
ds.tas
# %%
gdd = xclim.atmos.growing_degree_days(tas=ds.tas, thresh="10.0 degC", freq="YS")
gdd
# %%

gdd.sel(time="1990-01-01")
#%%


import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# %%
ds6h = xr.tutorial.open_dataset("air_temperature")
daily_ds = ds6h.resample(time="D").mean(keep_attrs=True)
with xclim.set_options(cf_compliance="log"):
    gdd = xclim.atmos.growing_degree_days(
        tas="air", thresh="10.0 degC", freq="MS", ds=daily_ds
    )

    # variable names default to xclim names, so we can even do this:
    renamed_daily_ds = daily_ds.rename(air="tas")
    gdd = xclim.atmos.growing_degree_days(
        thresh="10.0 degC", freq="MS", ds=renamed_daily_ds
    )

# %%
gdd.sel(time="2013-07").plot()
# %%
plt.imshow(gdd.sel(time="2013-07")[0].data)
# %%


import xarray as xr
import geowombat as gw
import os, sys

sys.path.append("/home/mmann1123/Documents/github/xr_fresh/")
from xr_fresh.feature_calculators import *
from xr_fresh.backends import Cluster
from xr_fresh.extractors import extract_features
from glob import glob
from datetime import datetime
import matplotlib.pyplot as plt
from xr_fresh.utils import *
import logging
import warnings
import xarray as xr
from numpy import where
from xr_fresh import feature_calculators
from itertools import chain
from geowombat.backends import concat as gw_concat

_logger = logging.getLogger(__name__)
from numpy import where
from xr_fresh.utils import xarray_to_rasterio
import pandas as pd
from pathlib import Path

#%% ppt time series meher growing season May to Feb

files = "/mnt/space/Dropbox/Ethiopia_data/Precip/original/"
band_name = "ppt"
file_glob = f"{files}/*.tif"
strp_glob = f"{files}precipitation_%Y%m%d_5000m.tif"


#%%  Precip


def list_2_years(file_path, year: str):
    fs = sorted(glob(file_path + "*" + str(year) + "*.tif")) + sorted(
        glob(file_path + "*" + str(year + 1) + "*.tif")
    )
    ds = sorted(datetime.strptime(string, strp_glob) for string in fs)
    return fs, ds


f_list = sorted(glob(file_glob))

dates = sorted(datetime.strptime(string, strp_glob) for string in f_list)


#%%
import glob
import pandas as pd
import xarray as xr


def time_index_from_filenames(filenames):
    """helper function to create a pandas DatetimeIndex
       Filename example: 20150520_0164.tif"""
    return pd.DatetimeIndex([pd.Timestamp(f[:8]) for f in filenames])


filenames = glob.glob("*.tif")
time = xr.Variable("time", time_index_from_filenames(filenames))
chunks = {"x": 5490, "y": 5490, "band": 1}
da = xr.concat([xr.open_rasterio(f, chunks=chunks) for f in filenames], dim=time)


#%%
from datetime import datetime

print(datetime.now())

# open xarray
with gw.config.update(ref_res=f_list[0], ref_bounds=f_list[0]):

    for year in sorted(list(set([x.year for x in dates])))[:]:
        file_year, date_year = list_2_years(files, year)
        with gw.open(
            file_year, band_names=[band_name], time_names=date_year, nodata=-9999,
        ) as ds:
            # print(ds)

            ds = ds.chunk(
                {"time": -1, "band": 1, "y": "auto", "x": "auto"}
            )  # rechunk to time

            # move dates back 2 months so year ends feb 29,
            # so month range now May = month 3,
            # feb of following year = month 12
            ds = ds.assign_coords(
                time=(pd.Series(ds.time.values) - pd.DateOffset(months=2)).values
            )

            # start cluster
            cluster = Cluster()
            cluster.start_large_object()

            # generate features
            ds_year = ds.sel(
                time=slice(str(year) + "-03-01", str(year) + "-07-29")
            )  # full '-03-01' to -12-29' - may to sept year+'-03-01', year+'-07-29'
            ds_year = ds_year.interpolate_na(dim="time")
            ds_year = ds_year.chunk(
                {"time": -1, "band": 1, "y": "auto", "x": "auto"}
            )  # rechunk to time

            print(ds_year)

# %%
