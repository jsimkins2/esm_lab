# General imports and definitions
import os, sys, datetime, logging
import xesmf as xe
import xarray as xr
import numpy as np

class TopoUtils:

    def __init__(self):
        # Constants
        self.PI_180 = np.pi/180.
        # Adopt GRS80 ellipse from proj
        self._default_Re = 6.378137e6
        self._default_ellps = 'GRS80'


    # Topography functions
    
    def regridTopo(gridFile, bathFile, gridGeoLoc = "corner", bathVarName="elevation", coarsenInt=10, superGrid=True, gridDimX=None, gridDimY=None, gridLatName=None, gridLonName=None, bathDimX=None, bathDimY=None, bathLatName=None, bathLonName=None):
        """Regrid topography file to the resolution of a given grid file.
        
        gridFile: Path to netCDF gridFile
        bathFile: Path to netCDF bathymetry file
        gridGeoLoc: "center" or "corner" - cell location of geographic placement of grid file - for example, x/y are located at 'q' points, which are located at the corner of each cell
        bathVarName: varialbe name of the bathymetry/depth/topography variable within the bathymetry file
        coarsenInt: Integer value used to decrease resolution of a given bathymetry file - see `xarray.coarsen`
        superGrid: When true, this assumes the gridFile is a supergrid and the resulting bathymetry is coarsened to a regular grid.
        gridDimX: The name of the dimension along the X axis of the grid file
        gridDimY: The name of the dimension along the Y axis of the grid file
        gridLatName: The name of the latitude variable of the grid file
        gridLonName: The name of the longitude variable of the grid file
        bathDimX: The name of the dimension along the X axis of the topography file
        bathDimY: The name of the dimension along the Y axis of the topography file
        bath
        
        """
        
        grid = xr.open_dataset(gridFile)

        if gridGeoLoc == "center":
            if 'nx' not in grid.dims:
                if gridDimX != None:
                    grid = grid.rename_dims({gridDimX : "nx"})
                elif 'lon' in grid.dims:
                    grid = grid.rename_dims({"lon" : "nx"})
                elif 'longitude' in grid.dims:
                    grid = grid.rename_dims({"longitude" : "nx"})
                elif 'x' in grid.dims:
                    grid = grid.rename_dims({"x" : "nx"})
                else:
                    print ('Error: plase define gridDimX')
            if 'ny' not in grid.dims:
                if gridDimX != None:
                    grid = grid.rename_dims({gridDimY : "ny"})
                elif 'lon' in grid.dims:
                    grid = grid.rename_dims({"lat" : "ny"})
                elif 'longitude' in grid.dims:
                    grid = grid.rename_dims({"latitude" : "ny"})
                elif 'y' in grid.dims:
                    grid = grid.rename_dims({"y" : "ny"})
                else:
                    print ('Error: plase define gridDimY')

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
            
            # create filler arrays to interpolate lon corner values onto.
            filler_array = np.zeros((lon_corners.shape[0] + 1, lon_corners.shape[1] + 1))
            filler_array[0:lon_corners.shape[0], 0:lon_corners.shape[1]] = lon_corners[:,:]
            extCol = np.append(np.diff(lon_corners[:,-1]), np.diff(lon_corners[:,-1])[-1])  + lon_corners[:,-1]
            extRow = np.append(np.diff(lon_corners[-1,:], axis=0), np.diff(lon_corners[-1,:], axis=0)[-1])  + lon_corners[-1,:]
            filler_array[0:extCol.shape[0], -1] = extCol
            filler_array[-1, 0:extRow.shape[0]] = extRow
            final_corner = np.diff(lon_corners[-1,:])[-1]*2  + lon_corners[-1,-1]
            filler_array[-1,-1] = final_corner

            grid['lon_corners'] = xr.DataArray(data=filler_array, dims=("nyp", "nxp"))
            
            lat_corners = 0.25 * (
                lat_centers[:-1, :-1]
                + lat_centers[1:, :-1]
                + lat_centers[:-1, 1:]
                + lat_centers[1:, 1:]
            )
            
            # create filler arrays to interpolate lat corner values onto.
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

            grid['lat_corners'] = xr.DataArray(data=filler_array, dims=("nyp", "nxp"))
            
        if gridGeoLoc == "corner":
            if 'nxp' not in grid.dims:
                if gridDimX != None:
                    grid = grid.rename_dims({gridDimX : "nxp"})
                elif 'lon' in grid.dims:
                    grid = grid.rename_dims({"lon" : "nxp"})
                elif 'longitude' in grid.dims:
                    grid = grid.rename_dims({"longitude" : "nxp"})
                elif 'x' in grid.dims:
                    grid = grid.rename_dims({"x" : "nxp"})
                else:
                    print ('Error: plase define gridDimX')
            if 'nyp' not in grid.dims:
                if gridDimX != None:
                    grid = grid.rename_dims({gridDimY : "nyp"})
                elif 'lon' in grid.dims:
                    grid = grid.rename_dims({"lat" : "nyp"})
                elif 'longitude' in grid.dims:
                    grid = grid.rename_dims({"latitude" : "nyp"})
                elif 'y' in grid.dims:
                    grid = grid.rename_dims({"y" : "nyp"})
                else:
                    print ('Error: plase define gridDimY')
                    
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


        # BATHYMETRY/TOPOGRAPHY XARRAY OBJECT ORGANIZATION
        bath = xr.open_dataset(bathFile)
        
        if 'nx' not in bath.dims:
            if bathDimX != None:
                bath = bath.rename_dims({bathDimX : "nx"})
            elif 'lon' in bath.dims:
                bath = bath.rename_dims({"lon" : "nx"})
            elif 'longitude' in bath.dims:
                bath = bath.rename_dims({"longitude" : "nx"})
            elif 'x' in bath.dims:
                bath = bath.rename_dims({"x" : "nx"})
            else:
                print ('Error: plase define bathDimX')
        if 'ny' not in bath.dims:
            if bathDimY != None:
                bath = bath.rename_dims({bathDimY : "ny"})
            elif 'lat' in bath.dims:
                bath = bath.rename_dims({"lat" : "ny"})
            elif 'latitude' in bath.dims:
                bath = bath.rename_dims({"latitude" : "ny"})
            elif 'y' in bath.dims:
                bath = bath.rename_dims({"y" : "ny"})
            else:
                print ('Error: plase define bathDimY')
                
        # coarsen bath file down based on coarsenInt
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

        latMinInd = find_nearest(array = bath.lat_centers.values, value = np.min(grid.lat_centers.values))
        latMaxInd = find_nearest(array = bath.lat_centers.values, value = np.max(grid.lat_centers.values))
        lonMinInd = find_nearest(array = bath.lon_centers.values, value = np.min(grid.lon_centers.values))
        lonMaxInd = find_nearest(array = bath.lon_centers.values, value = np.max(grid.lon_centers.values))
        
        if lonMinInd > lonMaxInd:
            temp = lonMinInd
            lonMinInd = lonMaxInd
            lonMaxInd = temp

        # slice the large bathymetry file down to the extents of the grid file + 1 on either side because we will slice down 2 points 
        # after the corner points are calculated
        bath = bath.isel(nx=slice(lonMinInd - 1, lonMaxInd + 1), ny=slice(latMinInd - 1, latMaxInd + 1))
        
        lon_centers = bath['lon_centers'].values
        lat_centers = bath['lat_centers'].values
        
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
        
        # extract the topo values and add them back later with proper dimensions
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

        
        # create xarray data array of the bathymetry variable name containing the topographic elevation data
        dr = bath[bathVarName]
        
        # create ocean fraction array where ocean cells are 1 and land cells are 0
        lm_ds = bath[bathVarName].where(bath[bathVarName] < 0)
        lm_ds = lm_ds.fillna(0)
        lm_ds = lm_ds.where(lm_ds > -0.000001)
        lm_ds = lm_ds.fillna(1)
        
        
        # regrid our bathymetry and land/ocean mask based on method
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
        
        return
    
    
    

