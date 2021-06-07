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

    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    # Topography functions
    def regridTopo(gridFile, topoFile, gridGeoLoc = "corner", topoVarName="elevation", coarsenInt=10, superGrid=True, gridDimX=None, gridDimY=None, gridLatName=None, gridLonName=None, topoDimX=None, topoDimY=None, topoLatName=None, topoLonName=None):
        """Regrid topography file to the grid of a given grid file. It is assumed that the topography file is on a rectangular grid
        and has a finer resolution than the grid file. It is also assumed that the topography values are defined at the cell centers.
        We recommend using GEBCO2020 topography dataset for the topography file. It can be found here: https://www.gebco.net/data_and_products/gridded_bathymetry_data/
        
        gridFile: Path to netCDF gridFile
        topoFile: Path to netCDF topography file
        gridGeoLoc: "center" or "corner" - cell location of geographic placement of grid file - for example, x/y are located at 'q' points, which are located at the corner of each cell
        topoVarName: varialbe name of the topography variable within the topography file
        coarsenInt: Integer value used to decrease resolution of a given topography file - see `xarray.coarsen`
        superGrid: When true, this assumes the gridFile is a supergrid and the resulting topography is coarsened to a regular grid.
        gridDimX: The name of the dimension along the X axis of the grid file
        gridDimY: The name of the dimension along the Y axis of the grid file
        gridLatName: The name of the latitude variable of the grid file
        gridLonName: The name of the longitude variable of the grid file
        topoDimX: The name of the dimension along the X axis of the topography file
        topoDimY: The name of the dimension along the Y axis of the topography file
        topoLatName: The name of the latitude variable within the topography file
        topoLonName: The name of the longitude variable within the topography file
        
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


        # TOPOGRAPHY XARRAY OBJECT ORGANIZATION
        # Note that we assume topography is on a rectangular grid
        topo = xr.open_dataset(topoFile)
        
        if 'nx' not in topo.dims:
            if topoDimX != None:
                topo = topo.rename_dims({topoDimX : "nx"})
            elif 'lon' in topo.dims:
                topo = topo.rename_dims({"lon" : "nx"})
            elif 'longitude' in topo.dims:
                topo = topo.rename_dims({"longitude" : "nx"})
            elif 'x' in topo.dims:
                topo = topo.rename_dims({"x" : "nx"})
            else:
                print ('Error: plase define topoDimX')
        if 'ny' not in topo.dims:
            if topoDimY != None:
                topo = topo.rename_dims({topoDimY : "ny"})
            elif 'lat' in topo.dims:
                topo = topo.rename_dims({"lat" : "ny"})
            elif 'latitude' in topo.dims:
                topo = topo.rename_dims({"latitude" : "ny"})
            elif 'y' in topo.dims:
                topo = topo.rename_dims({"y" : "ny"})
            else:
                print ('Error: plase define topoDimY')
                
        # coarsen topo file down based on coarsenInt
        topo = topo.coarsen(nx=coarsenInt,ny=coarsenInt, boundary='pad').mean()
        
        if 'lat_centers' not in topo.variables:
            if topoLatName != None:
                topo = topo.rename({topoLatName: 'lat_centers'})
            elif 'y' in topo.variables:
                topo = topo.rename({'y': 'lat_centers'})
            elif 'lat' in topo.variables:
                topo = topo.rename({'lat': 'lat_centers'})
            elif 'latitude' in topo.variables:
                topo = topo.rename({'latitude': 'lat_centers'})
            else:
                print('Error: please define gridlatname')
        if 'lon_centers' not in topo.variables:
            if topoLonName != None:
                topo = topo.rename({topoLonName : 'lon_centers'})
            elif 'x' in topo.variables:
                topo = topo.rename({'x': 'lon_centers'})
            elif 'lon' in topo.variables:
                topo = topo.rename({'lon': 'lon_centers'})
            elif 'longitude' in grid.variables:
                topo = topo.rename({'longitude': 'lon_centers'})
            else:
                print('Error: Please define gridLonName')

        latMinInd = find_nearest(array = topo.lat_centers.values, value = np.min(grid.lat_centers.values))
        latMaxInd = find_nearest(array = topo.lat_centers.values, value = np.max(grid.lat_centers.values))
        lonMinInd = find_nearest(array = topo.lon_centers.values, value = np.min(grid.lon_centers.values))
        lonMaxInd = find_nearest(array = topo.lon_centers.values, value = np.max(grid.lon_centers.values))
        
        if lonMinInd > lonMaxInd:
            temp = lonMinInd
            lonMinInd = lonMaxInd
            lonMaxInd = temp

        # slice the large topography file down to the extents of the grid file + 1 on either side because we will slice down 2 points 
        # after the corner points are calculated
        topo = topo.isel(nx=slice(lonMinInd - 1, lonMaxInd + 1), ny=slice(latMinInd - 1, latMaxInd + 1))
        
        lon_centers = topo['lon_centers'].values
        lat_centers = topo['lat_centers'].values
        
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
        topo = topo.isel(nx=slice(1,-1), ny=slice(1,-1))
        
        # extract the topo values and add them back later with proper dimensions
        elev = topo[topoVarName].values
        
        # add nxp and nyp dimensions for the lat/lon corners to latch onto
        topo = topo.expand_dims({'nyp':(len(topo.ny) + 1)})
        topo = topo.expand_dims({'nxp':(len(topo.nx) + 1)})
        
        # add the lat/lon corners as data variables
        topo['lat_corners'] = xr.DataArray(data=lat_corners, dims=("nyp"))
        topo['lon_corners'] = xr.DataArray(data=lon_corners, dims=("nxp"))
        
        # drop elevation and bring it back, this time constraining the dimensions to lat/lon centers
        topo = topo.drop_vars(topoVarName)
        topo[topoVarName] = (('ny', 'nx'), elev)
        

        lon2d, lat2d = np.meshgrid(topo.lon_centers.values, topo.lat_centers.values)
        lon2d_b, lat2d_b = np.meshgrid(topo.lon_corners.values, topo.lat_corners.values)
        
        # assign 2d coordinates as lat/lon 
        topo = topo.assign_coords({"lon" : (("ny", "nx"), lon2d)})
        topo = topo.assign_coords({"lat" : (("ny", "nx"), lat2d)})
        topo = topo.assign_coords({"lon_b" : (("nyp", "nxp"), lon2d_b)})
        topo = topo.assign_coords({"lat_b" : (("nyp", "nxp"), lat2d_b)})
        
        # rename for xesmf
        grid["lon"] = grid["lon_centers"]
        grid["lat"] = grid["lat_centers"]
        grid["lon_b"] = grid["lon_corners"]
        grid["lat_b"] = grid["lat_corners"]

        
        # create xarray data array of the topography variable name containing the topographic elevation data
        dr = topo[topoVarName]
        
        # create ocean fraction array where ocean cells are 1 and land cells are 0
        lm_ds = topo[topoVarName].where(topo[topoVarName] < 0)
        lm_ds = lm_ds.fillna(0)
        lm_ds = lm_ds.where(lm_ds > -0.000001)
        lm_ds = lm_ds.fillna(1)
        
        
        # regrid our topography and land/ocean mask based on method
        regridder = xe.Regridder(topo, grid, method="conservative")
        dr_out = regridder(dr)
        lm_ds_out = regridder(lm_ds)
        
        # coarsen the topography and landmask fraction from supergrid to regular grid supergrid is True
        dr_out = dr_out.coarsen(nx=2,ny=2, boundary='pad').mean()
        lm_ds_out = lm_ds_out.coarsen(nx=2,ny=2, boundary='pad').mean()
        
        # save our netCDF files
        opath = os.path.dirname(gridFile)
        dr_out.to_netcdf(opath + "/ocean_topog.nc")
        lm_ds_out.to_netcdf(opath + "/ocean_mask.nc")
        
        return
    
    
    

