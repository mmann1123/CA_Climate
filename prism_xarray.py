#%%
# https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Opening_GeoTIFFs_NetCDFs.html

import os
import xarray as xr
from glob import glob
from datetime import datetime
import xclim

# from xclim.testing import open_dataset
import geopandas as gpd
import pandas as pd

#%%
# example = open_dataset("ERA5/daily_surface_cancities_1990-1993.nc")

#%% Get all tmax data

mxtmp_path = (
    r"/mnt/space/Dropbox/USA_Data/prism/Daily_rasters_temp_maxmin/Maximum_temperature"
)
band_name = "mxtmp"
file_glob = f"{mxtmp_path}/**/*_bil.bil"
strp_glob = f"PRISM_tmax_stable_4kmD2_%Y%m%d_bil.bil"

mxtmp_f_list = sorted(glob(file_glob, recursive=True))
mxtmp_dates = sorted(
    datetime.strptime(os.path.basename(string), strp_glob) for string in mxtmp_f_list
)

gdd_out_path = r"/mnt/space/Dropbox/USA_Data/prism/xclim_indicators/growing_degree_days"

# f_list = f_list[: 365 * 10]
# dates = dates[: 365 * 10]

# states = gpd.read_file("/mnt/space/Dropbox/USA_Data/states/cb_2018_us_state_20m.shp")
# CA = states[states.STUSPS == "CA"]
# bounds = BoundingBox(
#     left=-124.409591, bottom=32.534156, right=-114.139055, top=42.009247
# )

#%%
# Load in and concatenate all individual GeoTIFFs
# Subset by coorindates
for year in range(1981, 2020):
    print(year)
    year_bool = (
        pd.Series(mxtmp_f_list)
        .str.startswith(
            os.path.join(mxtmp_path, "PRISM_tmax_stable_4kmD2_" + str(year))
        )
        .tolist()
    )

    mxtmp_year_files = [i for (i, v) in zip(mxtmp_f_list, year_bool) if v]
    mxtmp_year_dates = [i for (i, v) in zip(mxtmp_dates, year_bool) if v]

    ds = xr.concat(
        [xr.open_rasterio(i) for i in mxtmp_year_files], dim=mxtmp_year_dates
    )
    ds = ds.where((-125 < ds.x) & (ds.x < -115) & (31 < ds.y) & (ds.y < 43), drop=True)

    # Rename the variable to a more useful name
    ds = ds.rename({"concat_dim": "time"})

    # Covert our xarray.DataArray into a xarray.Dataset
    ds = ds.to_dataset("band")
    ds = ds.rename({1: band_name})
    ds.tmp.attrs = {"units": "degC"}

    gdd = xclim.atmos.growing_degree_days(tas=ds.tmp, thresh="10.0 degC", freq="YS")
    print(gdd)
    gdd.to_netcdf(os.path.join(gdd_out_path, "GDD_" + str(year) + ".nc"))


# %%
gdd.sel(time="1981-01-01").plot()

#%%
################################################

#%%
# Load in and concatenate all individual GeoTIFFs
# Subset by coorindates
for year in range(1981, 2020):
    print(year)
    year_bool = (
        pd.Series(f_list)
        .str.startswith(os.path.join(files, "PRISM_tmax_stable_4kmD2_" + str(year)))
        .tolist()
    )

    year_files = [i for (i, v) in zip(f_list, year_bool) if v]
    year_dates = [i for (i, v) in zip(dates, year_bool) if v]

    ds = xr.concat([xr.open_rasterio(i) for i in year_files], dim=year_dates)
    ds = ds.where((-125 < ds.x) & (ds.x < -115) & (31 < ds.y) & (ds.y < 43), drop=True)

    # Rename the variable to a more useful name
    ds = ds.rename({"concat_dim": "time"})

    # Covert our xarray.DataArray into a xarray.Dataset
    ds = ds.to_dataset("band")
    ds = ds.rename({1: "tmp"})
    ds.tmp.attrs = {"units": "degC"}

    gdd = xclim.atmos.heat_wave_frequency(tas=ds.tmp, thresh="10.0 degC", freq="YS")
    print(gdd)
    gdd.to_netcdf(os.path.join(out_path, "HWF_" + str(year) + ".nc"))


############################################################
# %% env pygisbookgw

import geowombat as gw
import os
import xarray as xr
from glob import glob
from datetime import datetime
import xclim
import geopandas as gpd
from rasterio.coords import BoundingBox

#%% Get all tmax data

states = gpd.read_file("/mnt/space/Dropbox/USA_Data/states/cb_2018_us_state_20m.shp")
CA = states[states.STUSPS == "CA"]

files = (
    r"/mnt/space/Dropbox/USA_Data/prism/Daily_rasters_temp_maxmin/Maximum_temperature"
)
band_name = "ppt"
file_glob = f"{files}/**/*_bil.bil"
strp_glob = f"PRISM_tmax_stable_4kmD2_%Y%m%d_bil.bil"

f_list = sorted(glob(file_glob, recursive=True))
dates = sorted(
    datetime.strptime(os.path.basename(string), strp_glob) for string in f_list
)

f_list = f_list[: 365 * 1]
dates = dates[: 365 * 1]

#%%

print(CA.bounds.values)
bounds = BoundingBox(
    left=-124.409591, bottom=32.534156, right=-114.139055, top=42.009247
)

with gw.config.update(ref_bounds=bounds):
    with gw.open(
        f_list,
        stack_dim="time",
        time_names=dates,
        band_names=["tmx"],
        nodata=-9999,
        num_workers=3,
    ) as src:
        print(src)
        ds = src.to_dataset("band")
        ds.tmx.attrs = {"units": "degC"}
        ds.to_netcdf("~/Desktop/test.nc")


#%%
ds = xr.open_dataset("~/Desktop/test.nc")
gdd = xclim.atmos.growing_degree_days(tas=ds.tmx, thresh="10.0 degC", freq="YS")
gdd


# %%
gdd.sel(time="1981-01-01").plot()
# %%
