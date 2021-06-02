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
    
    def regridTopo(gridFile, bathFile, gridGeoLoc = "corner", bathGeoLoc = "center", bathVarName="elevation", coarsenInt=10, superGrid=True, gridLatName=None, gridLonName=None, bathLatName=None, bathLonName=None):
        """Regrid topography file to the resolution of a given grid file.
        
        gridFile: Path to netCDF gridFile
        bathFile: Path to netCDF bathymetry file
        gridGeoLoc: "center" or "corner" - cell location of geographic placement of grid file - for example, x/y are located at 'q' points, which are located at the corner of each cell
        bathGeoLoc: "center" or "corner" - cell location of geographic placement of bathymetry file - for example, x/y are located at 'q' points, which are located at the corner of each cell
        bathVarName: varialbe name of the bathymetry/depth/topography variable within the bathymetry file
        coarsenInt: Integer value used to decrease resolution of a given bathymetry file - see `xarray.coarsen`
        superGrid: When true, this assumes the gridFile is a supergrid and the resulting bathymetry is coarsened to a regular grid.
        
        """
        
        grid = xr.open_dataset(gridFile)

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

        
        # open bathymetry file
        bath = xr.open_dataset(bathFile)
        bathGeoLoc = "center"
        if bathGeoLoc == "center":
            # have to rename the dimensions to nx, ny
            # will need to revisit this in case dimensions are not named lat/lon
            bath = bath.rename_dims({"lon" : "nx"})
            bath = bath.rename_dims({"lat" : "ny"})
            
            # add lat/lon centers to bath variables
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
                    print('Error: please define bathlatname')
            
                    
            if 'lon_centers' not in bath.variables:
                if bathLonName != None:
                    bath = bath.rename({bathLonName : 'lon_centers'})
                elif 'x' in bath.variables:
                    bath = bath.rename({'x': 'lon_centers'})
                elif 'lon' in bath.variables:
                    bath = bath.rename({'lon': 'lon_centers'})
                elif 'longitude' in bath.variables:
                    bath= bath.rename({'longitude': 'lon_centers'})
                else:
                    print('Error: Please define bathLonName')
            
            # grab index location of grid extents
            def find_nearest(array, value):
                array = np.asarray(array)
                idx = (np.abs(array - value)).argmin()
                return idx
            
            latMinInd = find_nearest(array = bath.lat_centers.values, value = np.min(grid.lat_centers.values))
            latMaxInd = find_nearest(array = bath.lat_centers.values, value = np.max(grid.lat_centers.values))
            lonMinInd = find_nearest(array = bath.lon_centers.values, value = np.min(grid.lon_centers.values))
            lonMaxInd = find_nearest(array = bath.lon_centers.values, value = np.max(grid.lon_centers.values))
            
            # slice the large bathymetry file down to extents PLUS COARSEN INT as corners array should be 1 GREATER than centers array
            # calculate corners from this
            # still not sure if we need to do this step, but just keeping jst in case
            bathCo = bath.sel(nx=slice(lonMinInd, lonMaxInd + coarsenInt*2), ny=slice(latMinInd, latMaxInd + coarsenInt*2))
            bathCo = bathCo.coarsen(nx=coarsenInt,ny=coarsenInt, boundary='pad').mean()
            # slice the large bathymetry file down to the extents of the grid file
            bath = bath.isel(nx=slice(lonMinInd, lonMaxInd), ny=slice(latMinInd, latMaxInd))
            # if we want to coarsen the bath file to make it lighter on the machine, do so here
            # I'm not sure how coarsening will affect nxp/nyp arrays..hopefully adding coarsenInt solves that
            bath = bath.coarsen(nx=coarsenInt,ny=coarsenInt, boundary='pad').mean()
        
            # create dimensions of nxp nyp for the corners to latch onto
            # first extract elevation because expanding dimensions automatically adds them to elevation variable for some reason
            elev = bath[bathVarName].values
            bath = bath.expand_dims({'nyp':(len(bath.ny) + 1)})
            bath = bath.expand_dims({'nxp':(len(bath.nx) + 1)})
            ##### Calculate the Lat/Lon Corner Values
            lon_centers = bathCo['lon_centers'].values
            lat_centers = bathCo['lat_centers'].values
            
            # To use conservative regidding, we need the cells corners. 
            # Since they are not provided, we are creating some using a crude approximation. 
            # NOTE THAT BECAUSE THIS IS A RECTANGULAR GRID, ARRAYS ARE 1D HERE AS OPPOSED TO 2D IN THE GRID FILE ABOVE
            lon_corners = 0.25 * (
                lon_centers[:-1]
                + lon_centers[1:]
                + lon_centers[:-1]
                + lon_centers[1:]
            )
        
            bath['lon_corners'] = xr.DataArray(data=lon_corners, dims=("nxp"))
            
            lat_corners = 0.25 * (
                lat_centers[:-1]
                + lat_centers[1:]
                + lat_centers[:-1]
                + lat_centers[1:]
            )
            
            # assign these as coordinates in future?
            bath['lat_corners'] = xr.DataArray(data=lat_corners, dims=("nyp"))
            
            # drop elevation and bring it back - have to do this because for some reason adding dimensinos 
            bath = bath.drop_vars(bathVarName)
            bath[bathVarName] = (('lon_centers', 'lat_centers'), elev)
        
        if bathGeoLoc == "corner":
            if 'lat_corners' not in bath.variables:
                if bathLatName != None:
                    bath = bath.rename({bathLatName: 'lat_corners'})
                elif 'y' in bath.variables:
                    bath = bath.rename({'y': 'lat_corners'})
                elif 'lat' in bath.variables:
                    bath = bath.rename({'latitude': 'lat_corners'})
                elif 'latitude' in bath.variables:
                    bath = bath.rename({'latitude': 'lat_corners'})
                else:
                    print('Error: please define gridlatname')
            
                    
            if 'lon_corners' not in bath.variables:
                if bathLonName != None:
                    bath = bath.rename({bathLonName : 'lon_corners'})
                elif 'x' in bath.variables:
                    bath = bath.rename({'x': 'lon_corners'})
                elif 'lon' in bath.variables:
                    bath = bath.rename({'lon': 'lon_corners'})
                elif 'longitude' in grid.variables:
                    bath = bath.rename({'longitude': 'lon_corners'})
                else:
                    print('Error: Please define gridLonName')
        
            bath = bath.sel(lon_centers=slice(np.min(grid.lon_centers.values), np.max(grid.lon_centers.values)), 
                          lat_centers=slice(np.min(grid.lat_centers.values), np.max(grid.lat_centers.values)))
            bath = bath.coarsen(lon_centers=coarsenInt,lat_centers=coarsenInt, boundary='pad').mean()    
            lon_corners = bath['lon_corners'].values
            lat_corners = bath['lat_corners'].values
            
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
            
            bath['lat_centers'] = xr.DataArray(data=lat_centers, dims=("ny", "nx"))
            bath['lon_centers'] = xr.DataArray(data=lon_centers, dims=("ny", "nx"))
        
        
        ### CONVERT RECTANGULAR GRID TO CURVILINEAR HERE
        # I am not sure if we need this curvilinear portion or not - test this
        # get 2D versions of the lat and lon variables
        #lon2d, lat2d = np.meshgrid(bath.lon_centers.values, bath.lat_centers.values)
        
        # assign 2d coordinates as lat/lon 
        #bath = bath.assign_coords({"lon" : (("ny", "nx"), lon2d)})
        #bath = bath.assign_coords({"lat" : (("ny", "nx"), lat2d)})
                
        # rename for xesmf
        grid["lon"] = grid["lon_centers"]
        grid["lat"] = grid["lat_centers"]
        grid["lon_b"] = grid["lon_corners"]
        grid["lat_b"] = grid["lat_corners"]
        
        bath["lon"] = bath["lon_centers"]
        bath["lat"] = bath["lat_centers"]
        bath["lon_b"] = bath["lon_corners"]
        bath["lat_b"] = bath["lat_corners"]
        
        # create xarray data array of the bathymetry variable name containing the topographic elevation data
        dr = bath[bathVarName]
        
        # create ocean fraction array where ocean cells are 1 and land cells are 0
        lm_ds = bath[bathVarName].where(bath[bathVarName] < 0)
        lm_ds = lm_ds.fillna(0)
        lm_ds = lm_ds.where(lm_ds > -0.000001)
        lm_ds = lm_ds.fillna(1)
        
        
        # regrid our bathymetry and land/ocean mask
        regridder = xe.Regridder(bath, grid, method="conservative")
        dr_out = regridder(dr)
        lm_ds_out = regridder(lm_ds)
        
        # coarsen the bathymetry and landmask fraction from supergrid to regular grid supergrid is True
        dr_out = dr_out.coarsen(nx=2,ny=2, boundary='pad').mean()
        lm_ds_out = lm_ds_out.coarsen(nx=2,ny=2, boundary='pad').mean()
        
        # save our netCDF files
        opath = os.path.dirname(gridFile)
        dr_out.to_netcdf(opath + "ocean_topog.nc")
        lm_ds_out.to_netcdf(opath + "ocean_mask.nc")
        
        return
    
    
    

