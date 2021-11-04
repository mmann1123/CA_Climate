#%%  pygisbookgw
# https://docs.dea.ga.gov.au/notebooks/Frequently_used_code/Opening_GeoTIFFs_NetCDFs.html

import os
import xarray as xr
from glob import glob
from datetime import datetime
import xclim

# from xclim.testing import open_dataset
import geopandas as gpd
import pandas as pd
import rioxarray  # for writing out xarray to tif

#%%
# example = open_dataset("ERA5/daily_surface_cancities_1990-1993.nc")

################################################
# Load in and concatenate all individual GeoTIFFs
# Subset by coorindates
for year in range(1981, 2020):
    print(year)

    # get tmax
    year_bool = (
        pd.Series(mxtmp_f_list)
        .str.startswith(
            os.path.join(mxtmp_path, "PRISM_tmax_stable_4kmD2_" + str(year))
        )
        .tolist()
    )

    mxtmp_year_files = [i for (i, v) in zip(mxtmp_f_list, year_bool) if v]
    mxtmp_year_dates = [i for (i, v) in zip(mxtmp_dates, year_bool) if v]

    ds_mxtmp = xr.concat(
        [xr.open_rasterio(i) for i in mxtmp_year_files], dim=mxtmp_year_dates
    )
    ds_mxtmp = ds_mxtmp.where(
        (-125 < ds_mxtmp.x)
        & (ds_mxtmp.x < -113)
        & (31 < ds_mxtmp.y)
        & (ds_mxtmp.y < 43),
        drop=True,
    )

    # Rename the variable to a more useful name
    ds_mxtmp = ds_mxtmp.rename({"concat_dim": "time"})

    # Covert our xarray.DataArray into a xarray.Dataset
    ds = ds_mxtmp.to_dataset(dim="band")
    ds = ds.rename({1: "tmax"})
    ds.tmax.attrs = {"units": "degC"}

    # replace -9999 with nan
    ds = ds.where(ds["tmax"] != -9999.0)

    #############################################
    # get tmin

    year_bool = (
        pd.Series(mntmp_f_list)
        .str.startswith(
            os.path.join(mntmp_path, "PRISM_tmin_stable_4kmD2_" + str(year))
        )
        .tolist()
    )

    mntmp_year_files = [i for (i, v) in zip(mntmp_f_list, year_bool) if v]
    mntmp_year_dates = [i for (i, v) in zip(mntmp_dates, year_bool) if v]

    ds_mntmp = xr.concat(
        [xr.open_rasterio(i) for i in mntmp_year_files], dim=mntmp_year_dates
    )
    ds_mntmp = ds_mntmp.where(
        (-125 < ds_mntmp.x)
        & (ds_mntmp.x < -113)
        & (31 < ds_mntmp.y)
        & (ds_mntmp.y < 43),
        drop=True,
    )

    # Rename the variable to a more useful name
    ds_mntmp = ds_mntmp.rename({"concat_dim": "time"})

    ds["tmin"] = ds_mntmp
    ds.tmin.attrs = {"units": "degC"}

    # replace -9999 with nan
    ds = ds.where(ds["tmin"] != -9999.0)

    #############################################
    # add average temp
    ds["tmean"] = ds.tmin / ds.tmax
    ds.tmean.attrs = {"units": "degC"}

    #############################################
    # calculate metrics
    hwf = xclim.atmos.heat_wave_frequency(
        tasmin=ds.tmin, tasmax=ds.tmax, window=3, freq="YS"
    )
    print(hwf)
    # hwf.to_netcdf(os.path.join(hwf_out_path, "HWF_" + str(year) + ".nc"))
    hwf.rio.write_crs("epsg:4326", inplace=True)
    hwf.sel(band=1).rio.to_raster(
        os.path.join(hwf_out_path, "HWF_" + str(year) + ".tif")
    )

    gdd = xclim.atmos.growing_degree_days(tas=ds.tmean, thresh="4.0 degC", freq="YS")
    print(gdd)
    gdd.rio.write_crs("epsg:4326", inplace=True)
    gdd.sel(band=1).rio.to_raster(
        os.path.join(gdd_out_path, "GDD_" + str(year) + ".tif")
    )


#%%
############################################################
# # %% env pygisbookgw

# import geowombat as gw
# import os
# import xarray as xr
# from glob import glob
# from datetime import datetime
# import xclim
# import geopandas as gpd
# from rasterio.coords import BoundingBox

# #%% Get all tmax data

# states = gpd.read_file("/mnt/space/Dropbox/USA_Data/states/cb_2018_us_state_20m.shp")
# CA = states[states.STUSPS == "CA"]

# files = (
#     r"/mnt/space/Dropbox/USA_Data/prism/Daily_rasters_temp_maxmin/Maximum_temperature"
# )
# band_name = "ppt"
# file_glob = f"{files}/**/*_bil.bil"
# strp_glob = f"PRISM_tmax_stable_4kmD2_%Y%m%d_bil.bil"

# f_list = sorted(glob(file_glob, recursive=True))
# dates = sorted(
#     datetime.strptime(os.path.basename(string), strp_glob) for string in f_list
# )

# f_list = f_list[: 365 * 1]
# dates = dates[: 365 * 1]

# #%%

# print(CA.bounds.values)
# bounds = BoundingBox(
#     left=-124.409591, bottom=32.534156, right=-114.139055, top=42.009247
# )

# with gw.config.update(ref_bounds=bounds):
#     with gw.open(
#         f_list,
#         stack_dim="time",
#         time_names=dates,
#         band_names=["tmx"],
#         nodata=-9999,
#         num_workers=3,
#     ) as src:
#         print(src)
#         ds = src.to_dataset("band")
#         ds.tmx.attrs = {"units": "degC"}
#         ds.to_netcdf("~/Desktop/test.nc")


# #%%
# ds = xr.open_dataset("~/Desktop/test.nc")
# gdd = xclim.atmos.growing_degree_days(tas=ds.tmx, thresh="10.0 degC", freq="YS")
# gdd


# # %%
# gdd.sel(time="1981-01-01").plot()
# # %%

# %%
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

#%% Min temp data

mntmp_path = (
    r"/mnt/space/Dropbox/USA_Data/prism/Daily_rasters_temp_maxmin/Minimum_temperature"
)

band_name = "mntmp"
file_glob = f"{mntmp_path}/**/*_bil.bil"
strp_glob = f"PRISM_tmin_stable_4kmD2_%Y%m%d_bil.bil"

mntmp_f_list = sorted(glob(file_glob, recursive=True))
mntmp_dates = sorted(
    datetime.strptime(os.path.basename(string), strp_glob) for string in mntmp_f_list
)

hwf_out_path = r"/mnt/space/Dropbox/USA_Data/prism/xclim_indicators/heat_wave_frequency"


# f_list = f_list[: 365 * 10]
# dates = dates[: 365 * 10]

# states = gpd.read_file("/mnt/space/Dropbox/USA_Data/states/cb_2018_us_state_20m.shp")
# CA = states[states.STUSPS == "CA"]
# bounds = BoundingBox(
#     left=-124.409591, bottom=32.534156, right=-114.139055, top=42.009247
# # )

# #%% GROWING DEGREE DAYS
# # Load in and concatenate all individual GeoTIFFs
# # Subset by coorindates
# for year in range(1981, 2020):
#     print(year)
#     year_bool = (
#         pd.Series(mxtmp_f_list)
#         .str.startswith(
#             os.path.join(mxtmp_path, "PRISM_tmax_stable_4kmD2_" + str(year))
#         )
#         .tolist()
#     )

#     mxtmp_year_files = [i for (i, v) in zip(mxtmp_f_list, year_bool) if v]
#     mxtmp_year_dates = [i for (i, v) in zip(mxtmp_dates, year_bool) if v]

#     ds_mxtmp = xr.concat(
#         [xr.open_rasterio(i) for i in mxtmp_year_files], dim=mxtmp_year_dates
#     )
#     ds_mxtmp = ds_mxtmp.where(
#         (-125 < ds_mxtmp.x)
#         & (ds_mxtmp.x < -113)
#         & (31 < ds_mxtmp.y)
#         & (ds_mxtmp.y < 43),
#         drop=True,
#     )

#     # Rename the variable to a more useful name
#     ds_mxtmp = ds_mxtmp.rename({"concat_dim": "time"})

#     # Covert our xarray.DataArray into a xarray.Dataset
#     ds_mxtmp = ds_mxtmp.to_dataset("band")
#     ds_mxtmp = ds_mxtmp.rename({1: "tmp"})
#     ds_mxtmp.tmp.attrs = {"units": "degC"}

#     gdd = xclim.atmos.growing_degree_days(
#         tas=ds_mxtmp.tmp, thresh="10.0 degC", freq="YS"
#     )
#     print(gdd)
#     # gdd.to_netcdf()
#     gdd.rio.write_crs("epsg:4326", inplace=True)
#     gdd.rio.to_raster(os.path.join(gdd_out_path, "GDD_" + str(year) + ".tif"))
# # gdd.sel(time="1981-01-01").plot()
