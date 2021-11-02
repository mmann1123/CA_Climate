#%%
# https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Opening_GeoTIFFs_NetCDFs.html

import os
import xarray as xr
from glob import glob
from datetime import datetime
import xclim
from xclim.testing import open_dataset

#%%
# example = open_dataset("ERA5/daily_surface_cancities_1990-1993.nc")

#%% Get all tmax data

files = r"/mnt/space/Dropbox/USA_Data/Daily_rasters_temp_maxmin/Maximum_temperature"
band_name = "ppt"
file_glob = f"{files}/**/*_bil.bil"
strp_glob = f"PRISM_tmax_stable_4kmD2_%Y%m%d_bil.bil"

f_list = sorted(glob(file_glob, recursive=True))
dates = sorted(
    datetime.strptime(os.path.basename(string), strp_glob) for string in f_list
)

f_list = f_list[: 365 * 5]
dates = dates[: 365 * 5]

#%%
# Load in and concatenate all individual GeoTIFFs
ds = xr.concat([xr.open_rasterio(i) for i in f_list], dim=dates)

# Rename the variable to a more useful name
ds = ds.rename({"concat_dim": "time"})

# Covert our xarray.DataArray into a xarray.Dataset
ds = ds.to_dataset("band")
ds = ds.rename({1: "tmp"})
ds.tmp.attrs = {"units": "degC"}

print(ds)
# %%
gdd = xclim.atmos.growing_degree_days(tas=ds.tmp, thresh="10.0 degC", freq="YS")

# %%
import matplotlib.pyplot as plt

plt.plot(gdd.sel(time="1981-01-01").values)
plt.show()


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

f_list = f_list[: 365 * 5]
dates = dates[: 365 * 5]

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
        num_workers=2,
    ) as src:
        print(src)
        ds = src.to_dataset("band")
        ds.tmx.attrs = {"units": "degC"}
        gdd = xclim.atmos.growing_degree_days(tas=ds.tmx, thresh="10.0 degC", freq="YS")
        gdd.compute()


# %%
import matplotlib.pyplot as plt

fig, ax = plt.subplots(dpi=200)
plt.plot(gdd.sel(time="1981-01-01").values)
plt.show()

# %%
