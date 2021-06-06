# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import xesmf as xe
import os


gridFile="/Users/james/Documents/Github/esm_lab/gridTools/nep7_grid/ocean_hgrid.nc"
#gridFile = "/Users/james/Downloads/gridFile.nc"
bathFile="/Users/james/Downloads/gebco_2020_netcdf/GEBCO_2020.nc"
gridGeoLoc = "corner"
bathGeoLoc = "center"
bathVarName = 'elevation'
gridLatName = None
gridLonName = None
bathLatName = None
bathLonName = None
coarsenInt = 8
grid = xr.open_dataset(gridFile)
gridGeoLoc = 'corner'
if gridGeoLoc == "center":
    if 'lat_centers' not in grid.variables:
        if gridLatName != None:
            grid = grid.rename({gridLatName: 'lat_centers'})
        elif 'y' in grid.variables:
            grid = grid.rename({'y': 'lat_centers'})
        elif 'lat' in grid.variables:
            grid = grid.rename({'lat': 'lat_centers'})
        elif 'latitude' in grid.variables:
            grid = grid.rename({'latitude': 'lat_centers'})
        else:
            print('Error: please define gridlatname')
    
            
    if 'lon_centers' not in grid.variables:
        if gridLonName != None:
            grid = grid.rename({gridLonName : 'lon_centers'})
        elif 'x' in grid.variables:
            grid = grid.rename({'x': 'lon_centers'})
        elif 'lon' in grid.variables:
            grid = grid.rename({'lon': 'lon_centers'})
        elif 'longitude' in grid.variables:
            grid= grid.rename({'longitude': 'lon_centers'})
        else:
            print('Error: Please define gridLonName')
    
    # fix longitudes from 0 to 360 for computation
    grid['lon_centers'].values =  np.where(grid['lon_centers'].values < 0., grid['lon_centers'].values + 360, grid['lon_centers'].values)
    lon_centers = grid['lon_centers'].values
    lat_centers = grid['lat_centers'].values
    
    # To use conservative regidding, we need the cells corners. 
    # Since they are not provided, we are creating some using a crude approximation. 
    lon_corners = 0.25 * (
        lon_centers[:-1, :-1]
        + lon_centers[1:, :-1]
        + lon_centers[:-1, 1:]
        + lon_centers[1:, 1:]
    )
    
    # trying to expand out the cornesr here by taking the difference of the last row/column, 
    # adding the difference to the last row/column, then plugging that new row/column into an expanded matrix called filler_array
    filler_array = np.zeros((lon_corners.shape[0] + 1, lon_corners.shape[1] + 1))
    filler_array[0:lon_corners.shape[0], 0:lon_corners.shape[1]] = lon_corners[:,:]
    extCol = np.append(np.diff(lon_corners[:,-1]), np.diff(lon_corners[:,-1])[-1])  + lon_corners[:,-1]
    extRow = np.append(np.diff(lon_corners[-1,:], axis=0), np.diff(lon_corners[-1,:], axis=0)[-1])  + lon_corners[-1,:]
    filler_array[0:extCol.shape[0], -1] = extCol
    filler_array[-1, 0:extRow.shape[0]] = extRow
    # there is one final corner that's not been interpolated, fill that cell here - have to double the difference to avoid repeating 
    # the same value for that final corner
    final_corner = np.diff(lon_corners[-1,:])[-1]*2  + lon_corners[-1,-1]
    filler_array[-1,-1] = final_corner
    # testing
    # plt.pcolormesh(filler_array)
    grid['lon_corners'] = xr.DataArray(data=filler_array, dims=("nyp", "nxp"))
    
    lat_corners = 0.25 * (
        lat_centers[:-1, :-1]
        + lat_centers[1:, :-1]
        + lat_centers[:-1, 1:]
        + lat_centers[1:, 1:]
    )
    
    
    # trying to expand out the cornesr here by taking the difference of the last row/column, 
    # adding the difference to the last row/column, then plugging that new row/column into an expanded matrix called filler_array
    filler_array = np.zeros((lat_corners.shape[0] + 1, lat_corners.shape[1] + 1))
    filler_array[0:lat_corners.shape[0], 0:lat_corners.shape[1]] = lat_corners[:,:]
    extCol = np.append(np.diff(lat_corners[:,-1]), np.diff(lat_corners[:,-1])[-1])  + lat_corners[:,-1]
    extRow = np.append(np.diff(lat_corners[-1,:], axis=0), np.diff(lat_corners[-1,:], axis=0)[-1])  + lat_corners[-1,:]
    filler_array[0:extCol.shape[0], -1] = extCol
    filler_array[-1, 0:extRow.shape[0]] = extRow
    # there is one final corner that's not been interpolated, fill that cell here - have to double the difference to avoid repeating 
    # the same value for that final corner
    final_corner = np.diff(lat_corners[-1,:])[-1]*2  + lat_corners[-1,-1]
    filler_array[-1,-1] = final_corner
    # testing
    #plt.pcolormesh(filler_array)
    # add formal lat corner to grid object
    grid['lat_corners'] = xr.DataArray(data=filler_array, dims=("nyp", "nxp"))
    
if gridGeoLoc == "corner":
    if 'lat_corners' not in grid.variables:
        if gridLatName != None:
            grid = grid.rename({gridLatName: 'lat_corners'})
        elif 'y' in grid.variables:
            grid = grid.rename({'y': 'lat_corners'})
        elif 'lat' in grid.variables:
            grid = grid.rename({'lat': 'lat_corners'})
        elif 'latitude' in grid.variables:
            grid = grid.rename({'latitude': 'lat_corners'})
        else:
            print('Error: please define gridlatname')
    
            
    if 'lon_corners' not in grid.variables:
        if gridLonName != None:
            grid = grid.rename({gridLonName : 'lon_corners'})
        elif 'x' in grid.variables:
            grid = grid.rename({'x': 'lon_corners'})
        elif 'lon' in grid.variables:
            grid = grid.rename({'lon': 'lon_corners'})
        elif 'longitude' in grid.variables:
            grid= grid.rename({'longitude': 'lon_corners'})
        else:
            print('Error: Please define gridLonName')

    # fix longitudes from 0 to 360 for computation
    #if "lon_corners" in grid.coords:
     #   grid = grid.assign_coords(lon_corners=(np.where(grid['lon_corners'].values < 0., grid['lon_corners'].values + 360, grid['lon_corners'].values)))
     #   grid = grid.swap_dims({'lon_corners' : 'nxp'})    
    #if "lon_corners" in grid.data_vars:
      #  grid['lon_corners'].values =  np.where(grid['lon_corners'].values < 0., grid['lon_corners'].values + 360, grid['lon_corners'].values)

    lon_corners = grid['lon_corners'].values
    lat_corners = grid['lat_corners'].values
    
    # To use conservative regidding, we need the cells centers. 
    # Since they are not provided, we are creating some using a crude approximation. 
    lon_centers = 0.25 * (
        lon_corners[:-1, :-1]
        + lon_corners[1:, :-1]
        + lon_corners[:-1, 1:]
        + lon_corners[1:, 1:]
    )
    
    lat_centers = 0.25 * (
        lat_corners[:-1, :-1]
        + lat_corners[1:, :-1]
        + lat_corners[:-1, 1:]
        + lat_corners[1:, 1:]
    )
    
    grid['lat_centers'] = xr.DataArray(data=lat_centers, dims=("ny", "nx"))
    grid['lon_centers'] = xr.DataArray(data=lon_centers, dims=("ny", "nx"))
    

# note here that we will automatically declare that the bath/topo grid is defined by the centers

bath = xr.open_dataset(bathFile)
bath = bath.rename_dims({"lon" : "nx"})
bath = bath.rename_dims({"lat" : "ny"})
# coarsen bath file down 
bath = bath.coarsen(nx=coarsenInt,ny=coarsenInt, boundary='pad').mean()

if 'lat_centers' not in bath.variables:
    if bathLatName != None:
        bath = bath.rename({bathLatName: 'lat_centers'})
    elif 'y' in bath.variables:
        bath = bath.rename({'y': 'lat_centers'})
    elif 'lat' in bath.variables:
        bath = bath.rename({'lat': 'lat_centers'})
    elif 'latitude' in bath.variables:
        bath = bath.rename({'latitude': 'lat_centers'})
    else:
        print('Error: please define gridlatname')

        
if 'lon_centers' not in bath.variables:
    if bathLonName != None:
        bath = bath.rename({bathLonName : 'lon_centers'})
    elif 'x' in bath.variables:
        bath = bath.rename({'x': 'lon_centers'})
    elif 'lon' in bath.variables:
        bath = bath.rename({'lon': 'lon_centers'})
    elif 'longitude' in grid.variables:
        bath = bath.rename({'longitude': 'lon_centers'})
    else:
        print('Error: Please define gridLonName')


# grab index location of grid extents
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# fix longitudes from 0 to 360 for computation
#if "lon_centers" in bath.coords:
#    bath = bath.assign_coords(lon_centers=(np.where(bath['lon_centers'].values < 0., bath['lon_centers'].values + 360, bath['lon_centers'].values)))
#    bath = bath.swap_dims({'lon_centers' : 'nx'})
#if "lon_centers" in bath.data_vars:
#    bath['lon_centers'].values =  np.where(bath['lon_centers'].values < 0., bath['lon_centers'].values + 360, bath['lon_centers'].values)


latMinInd = find_nearest(array = bath.lat_centers.values, value = np.min(grid.lat_centers.values))
latMaxInd = find_nearest(array = bath.lat_centers.values, value = np.max(grid.lat_centers.values))
lonMinInd = find_nearest(array = bath.lon_centers.values, value = np.min(grid.lon_centers.values))
lonMaxInd = find_nearest(array = bath.lon_centers.values, value = np.max(grid.lon_centers.values))

# after we go 0-360 for longitude, the indices need to be swapped so we have a slice
if lonMinInd > lonMaxInd:
    temp = lonMinInd
    lonMinInd = lonMaxInd
    lonMaxInd = temp
# slice the large bathymetry file down to extents PLUS COARSEN INT as corners array should be 1 GREATER than centers array
# calculate corners from this
# still not sure if we need to do this step, but just keeping jst in case

# slice the large bathymetry file down to the extents of the grid file + 1 on either side because we will slice down 2 points 
# after the corner points are calculated
bath = bath.isel(nx=slice(lonMinInd - 1, lonMaxInd + 1), ny=slice(latMinInd - 1, latMaxInd + 1))





lon_centers = bath['lon_centers'].values
lat_centers = bath['lat_centers'].values

# To use conservative regidding, we need the cells centers. 
# Since they are not provided, we are creating some using a crude approximation. 
# this works because our bathymetry file is a rectangular grid - 
# otherwise the indexing of this step would like it does in the grid center point creation
lon_corners = 0.25 * (
    lon_centers[:-1]
    + lon_centers[1:]
    + lon_centers[:-1]
    + lon_centers[1:]
)

lat_corners = 0.25 * (
    lat_centers[:-1]
    + lat_centers[1:]
    + lat_centers[:-1]
    + lat_centers[1:]
)

# trim down the centers so they are 1 less than the corner points we just calculated
bath = bath.isel(nx=slice(1,-1), ny=slice(1,-1))

# extract the topo values because we will need to add them back in later
# after we expand dimensions, the dimensions are added to teh elevation variable for soem reason
elev = bath[bathVarName].values

# add nxp and nyp dimensions for the lat/lon corners to latch onto
bath = bath.expand_dims({'nyp':(len(bath.ny) + 1)})
bath = bath.expand_dims({'nxp':(len(bath.nx) + 1)})

# add the lat/lon corners as data variables
bath['lat_corners'] = xr.DataArray(data=lat_corners, dims=("nyp"))
bath['lon_corners'] = xr.DataArray(data=lon_corners, dims=("nxp"))

# drop elevation and bring it back, this time constraining the dimensions to lat/lon centers
bath = bath.drop_vars(bathVarName)
bath[bathVarName] = (('ny', 'nx'), elev)


# I am not sure if we need this curvilinear portion or not - test this
# get 2D versions of the lat and lon variables
lon2d, lat2d = np.meshgrid(bath.lon_centers.values, bath.lat_centers.values)
lon2d_b, lat2d_b = np.meshgrid(bath.lon_corners.values, bath.lat_corners.values)
# assign 2d coordinates as lat/lon 
bath = bath.assign_coords({"lon" : (("ny", "nx"), lon2d)})
bath = bath.assign_coords({"lat" : (("ny", "nx"), lat2d)})
bath = bath.assign_coords({"lon_b" : (("nyp", "nxp"), lon2d_b)})
bath = bath.assign_coords({"lat_b" : (("nyp", "nxp"), lat2d_b)})

# rename for xesmf
grid["lon"] = grid["lon_centers"]
grid["lat"] = grid["lat_centers"]
grid["lon_b"] = grid["lon_corners"]
grid["lat_b"] = grid["lat_corners"]

#bath["lon"] = bath["lon_centers"]
#bath["lat"] = bath["lat_centers"]
#bath["lon_b"] = bath["lon_corners"]
#bath["lat_b"] = bath["lat_corners"]


# create xarray data array of the bathymetry variable name containing the topographic elevation data
dr = bath[bathVarName]

# create land mask xarray dataarray
lm_ds = bath[bathVarName].where(bath[bathVarName] < 0)
lm_ds = lm_ds.fillna(0)
lm_ds = lm_ds.where(lm_ds > -0.000001)
lm_ds = lm_ds.fillna(1)
lm_ds.name = 'ocean_fraction'

# regrid our bathymetry and landmask
regridder = xe.Regridder(bath, grid, method="conservative")
dr_out = regridder(dr)
lm_ds_out = regridder(lm_ds)


# coarsen the bathymetry and landmask fraction from supergrid to regular grid supergrid is True
dr_out = dr_out.coarsen(nx=2,ny=2, boundary='pad').mean()
lm_ds_out = lm_ds_out.coarsen(nx=2,ny=2, boundary='pad').mean()

# save our netCDF files
opath = os.path.dirname(gridFile)
dr_out.to_netcdf(opath + "/ocean_topog.nc")
lm_ds_out.to_netcdf(opath + "/ocean_mask.nc")
        
        
        

        
        
