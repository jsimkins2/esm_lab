#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 13:31:16 2021

@author: james
"""

import xesmf
import xarray as xr
import numpy as np


ds = xr.tutorial.open_dataset("ROMS_example.nc", chunks={"ocean_time": 1})
lon_centers = ds["lon_rho"].values
lat_centers = ds["lat_rho"].values

# To use conservative regidding, we need the cells corners. 
# Since they are not provided, we are creating some using a crude approximation. 
lon_corners = 0.25 * (
    lon_centers[:-1, :-1]
    + lon_centers[1:, :-1]
    + lon_centers[:-1, 1:]
    + lon_centers[1:, 1:]
)

lat_corners = 0.25 * (
    lat_centers[:-1, :-1]
    + lat_centers[1:, :-1]
    + lat_centers[:-1, 1:]
    + lat_centers[1:, 1:]
)

ds["lon_psi"] = xr.DataArray(data=lon_corners, dims=("eta_psi", "xi_psi"))
ds["lat_psi"] = xr.DataArray(data=lat_corners, dims=("eta_psi", "xi_psi"))

ds = ds.assign_coords({"lon_psi": ds["lon_psi"], "lat_psi": ds["lat_psi"]})

ds = ds.isel(
    eta_rho=slice(1, -10),
    xi_rho=slice(1, -10),
    eta_psi=slice(0, -9),
    xi_psi=slice(0, -9),
)

# We also need a coarse resolution grid. We’re going to build one by coarsening the ROMS dataset. 
# coarsen.mean() typically works as a nan-mean on the 10x10 blocks of the grid so the resulting land mask 
# looks like a flooded version of the original.
ds_coarse = xr.Dataset()

ds_coarse["zeta"] = xr.DataArray(
    ds["zeta"].coarsen(xi_rho=10, eta_rho=10).mean().values,
    dims=("ocean_time", "eta_rho", "xi_rho"),
)
# we want to subsample coordinates instead of coarsening them
ds_coarse["lon_rho"] = xr.DataArray(
    ds["lon_rho"].values[::10, ::10], dims=("eta_rho", "xi_rho")
)
ds_coarse["lon_psi"] = xr.DataArray(
    ds["lon_psi"].values[::10, ::10], dims=("eta_psi", "xi_psi")
)
ds_coarse["lat_rho"] = xr.DataArray(
    ds["lat_rho"].values[::10, ::10], dims=("eta_rho", "xi_rho")
)
ds_coarse["lat_psi"] = xr.DataArray(
    ds["lat_psi"].values[::10, ::10], dims=("eta_psi", "xi_psi")
)



# As usual, xESMF expects fixed variable names for longitude/latitude in cell centers and corners
    
ds["lon"] = ds["lon_rho"]
ds["lat"] = ds["lat_rho"]
ds["lon_b"] = ds["lon_psi"]
ds["lat_b"] = ds["lat_psi"]

ds_coarse["lon"] = ds_coarse["lon_rho"]
ds_coarse["lat"] = ds_coarse["lat_rho"]
ds_coarse["lon_b"] = ds_coarse["lon_psi"]
ds_coarse["lat_b"] = ds_coarse["lat_psi"]
    
    
    
regrid_nomask = xesmf.Regridder(ds, ds_coarse, method="conservative")

zeta_remapped = regrid_nomask(ds["zeta"])

zeta_remapped.isel(ocean_time=0).plot()

import matplotlib.pyplot as plt
plt.spy(regrid_nomask.weights)
plt.xlabel("input grid indices")
plt.ylabel("output grid indices")


