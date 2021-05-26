# General imports and definitions
import os, sys, datetime, logging
import xesmf as xe
import xarray as xr
import numpy as np

class BathUtils:

    def __init__(self):
        # Constants
        self.PI_180 = np.pi/180.
        # Adopt GRS80 ellipse from proj
        self._default_Re = 6.378137e6
        self._default_ellps = 'GRS80'


    # Topography functions
    
    def regridTopo(gridFile, bathFile, geoLoc="corner", bathVarName="elevation", coarsenInt=1, superGrid=True, gridLatName=None, gridLonName=None, bathLatName=None, bathLonName=None):
        """Regrid topography file to the resolution of a given grid file.
        
        gridFile: Path to netCDF gridFile
        bathFile: Path to netCDF bathymetry file
        cellLoc: "center" or "corner" - cell location of geographic placement - for example, x/y are located at 'q' points, which are located at the corner of each cell
        bathVarName: varialbe name of the bathymetry/depth/topography variable within the bathymetry file
        coarsenInt: Integer value used to decrease resolution of a given bathymetry file - see `xarray.coarsen`
        superGrid: When true, this assumes the gridFile is a supergrid and the resulting bathymetry is coarsened to a regular grid.
        
        """
        
        grid = xr.open_dataset(gridFile)
        
        if 'lat' not in grid.variables:
            if 'y' in grid.variables:
                grid  = grid.assign_coords({"lat" : (("y", "x"), grid.y.values)})
            if 'latitude' in grid.variables:
                grid  = grid.assign_coords({"lat" : (("latitude", "longitude"), grid.latitude.values)})
            if gridLatName != None:
                grid  = grid.assign_coords({"lat" : ((gridLatName, gridLonName), grid[gridLatName].values)})

        if 'y' not in grid.variables:
            if 'lat' in grid.variables:
                grid  = grid.assign_coords({"y" : (("lat", "lon"), grid.lat.values)})
            if 'latitude' in grid.variables:
                grid  = grid.assign_coords({"y" : (("latitude", "longitude"), grid.latitude.values)})
            if gridLatName != None:
                grid  = grid.assign_coords({"y" : ((gridLatName, gridLonName), grid[gridLatName].values)})
                
        if 'lon' not in grid.variables:
            if 'x' in grid.variables:
                grid  = grid.assign_coords({"lon" : (("y", "x"), grid.x.values)})
            if 'longitude' in grid.variables:
                grid  = grid.assign_coords({"lon" : (("latitude", "longitude"), grid.longitude.values)})
            if gridLonName != None:
                grid  = grid.assign_coords({"lon" : ((gridLatName, gridLonName), grid[gridLonName].values)})

        if 'x' not in grid.variables:
            if 'lon' in grid.variables:
                grid  = grid.assign_coords({"x" : (("lat", "lon"), grid.lon.values)})
            if 'longitude' in grid.variables:
                grid  = grid.assign_coords({"x" : (("latitude", "longitude"), grid.longitude.values)})
            if gridLonName != None:
                grid  = grid.assign_coords({"x" : ((gridLatName, gridLonName), grid[gridLonName].values)})

        
        # open bathymetry file
        bathy = xr.open_dataset(bathFile)

        
        if 'lat' not in bathy.variables:
            if 'y' in bathy.variables:
                bathy  = bathy.assign_coords({"lat" : (("y"), bathy.y.values)})
                bathy.rename({'y': 'y_i'})
            if 'latitude' in bathy.variables:
                bathy  = bathy.assign_coords({"lat" : (("latitude", "longitude"), bathy.latitude.values)})
            if bathLatName != None:
                bathy  = bathy.assign_coords({"lat" : ((bathLatName, bathLonName), bathy[bathLatName].values)})
                
        if 'lon' not in bathy.variables:
            if 'x' in bathy.variables:
                bathy  = bathy.assign_coords({"lon" : (("x"), bathy.x.values)})
                bathy = bathy.rename({'x': 'x_i'})
            if 'longitude' in bathy.variables:
                bathy  = bathy.assign_coords({"lon" : (("latitude"), bathy.longitude.values)})
            if bathLonName != None:
                bathy  = bathy.assign_coords({"lon" : ((bathLatName, bathLonName), bathy[bathLonName].values)})
                

        bathy = bathy.rename({'lon': 'x', 'lat': 'y'})
        bathy = bathy.sel(y=slice(x=slice(np.min(grid.lon.values), np.max(grid.lon.values), np.min(grid.lat.values), np.max(grid.lat.values))))
        bathy = bathy.coarsen(x=coarsenInt,y=coarsenInt, boundary='pad').mean()
        
        
        # get 2D versions of the lat and lon variables
        lon2d, lat2d = np.meshgrid(bathy.x.values, bathy.y.values)
        
        # assign 2d coordinates as lat/lon 
        bathy = bathy.assign_coords({"lon" : (("x", "y"), lon2d)})
        bathy = bathy.assign_coords({"lat" : (("x", "y"), lat2d)})
        
        # create xarray data array of the bathymetry variable name containing the topographic elevation data
        dr = bathy[bathVarName]
        
        # create land mask file 
        lm_ds = bathy[bathVarName].where(bathy[bathVarName] < 0)
        lm_ds = lm_ds.fillna(1)
        lm_ds = bathy[bathVarName].where(bathy[bathVarName] > 0)
        lm_ds = lm_ds.fillna(0)
        
        # regrid our bathymetry and landmask
        regridder = xe.Regridder(bathy, grid, method="conservative")
        dr_out = regridder(dr)
        lm_ds_out = regridder(lm_ds)
        
        # coarsen the bathymetry and landmask fraction from supergrid to regular grid supergrid is True
        dr_out = dr_out.coarsen(lon=2,lat=2, boundary='pad').mean()
        lm_ds_out = lm_ds_out.coarsen(lon=2,lat=2, boundary='pad').mean()
        
        # save our netCDF files
        opath = os.path.dirname(bathFile)
        dr_out.to_netcdf(opath)
        lm_ds_out.to_netcdf(opath)
        
        return
    
    
    

