# Lambert Conformal Conic grid generation provided by:
#  https://github.com/nikizadehgfdl/grid_generation/blob/dev/jupynotebooks/regional_grid_spherical.ipynb

# TODO:
#   [DONE] Revise for https://www.python.org/dev/peps/pep-0008/#package-and-module-names
#   [DONE] Refactor matplotlib code to use Matplotlib panes
#      [DONE] Requires refactoring subplot to subsplots. See: https://towardsdatascience.com/plt-subplot-or-plt-subplots-understanding-state-based-vs-object-oriented-programming-in-pyplot-4ba0c7283f5d
#      [X] Another good example: https://unidata.github.io/MetPy/latest/examples/Four_Panel_Map.html
#   [DONE] Allow display of important messages and warnings in panel application
#   [] Further consolidate matplotlib code
#   [] Do we have to declare everything in __init__ first or can be push all that to respective reset/clear functions?
#   [] refactor refineS and refineR options as Niki had them defined
#   [] readGrid() and makeGrid() make use of xarray structures instead so we can utilize to_netcdf() serialization
#   [] Pass back an error graphic instead of None for figures that do not render
#   [] Add a formal logging/message mechanism.
#      [] Move all this to options.  Interact with message buffer.
#      [] Maybe warnings are better? Try some out.
#      [] Create a message buffer/system for information 
#   [] Bring in code that converts ROMS grids to MOM6 grids
#      [] Allow conversion of MOM6 grids to ROMS grids
#   [] For now, the gridParameters are always in reference to a center point in a grid
#     in the future, one may fix a side or point of the grid and grow out from that point
#     instead of the center.
#   [] Refactor plotting code.  It is mostly the same except for setting the projection.
#   [] makeGrid assumes degrees at this point.
#   [] grid reading and plot parameter defaults should be dynamic with grid type declaration and potentially
#      split out into separate library modules? lib/gridTools/grids/{MOM6,ROMS,WRF}

# Adopt:
#   https://www.python.org/dev/peps/pep-0008/#package-and-module-names
#   Function and variables are mixed type

# General imports and definitions
import os, sys
import cartopy
import numpy as np
import xarray as xr
import warnings
import pdb

# Needed for panel.pane                
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvas  # not needed for mpl >= 3.1

# Required for ROMS->MOM6 conversion
import Spherical

class GridUtils:

    def __init__(self, app={}):
        # Constants
        self.PI_180 = np.pi/180.
        self._default_Re = 6.378e6
        # File pointer
        self.ncfp = None
        # Internal parameters
        self.usePaneMatplotlib = False
        self.matplotlibImports = False
        # Debugging/logging
        self.debugLevel = 0
        self.verboseLevel = 0
        # Grid 
        self.grid = {}
        # Private variables begin with a _
        # Grid parameters
        self.grid['_meta'] = {}
        self.grid['_meta']['dimensions'] = {}
        self.grid['_meta']['gridParameters'] = {}
        self.grid['_meta']['gridParameterKeys'] = self.grid['_meta']['gridParameters'].keys()
        # Defaults
        self.plotParameterDefaults = {
            'figsize': (6.4, 4.8),
            'extent': [],
            'extentCRS': cartopy.crs.PlateCarree(),
            'projection': {
            },
            'showGrid': True,
            'showGridCells': False,
            'showSupergrid': False
        }          
        # Plot parameters
        self.grid['_meta']['plotParameters'] = self.plotParameterDefaults
        self.grid['_meta']['plotParameterKeys'] = self.grid['_meta']['plotParameters'].keys()
               
    # Utility functions
    
    def application(self, app={}):
        '''Attach application items to the GridUtil object.
        
            app = {
                'messages': panel.widget.TextBox     # Generally a pointer to a panel widget for display of text
                'defaultFigureSize': (8,6)           # Default figure size to return from matplotlib
                'usePaneMatplotlib': True/False      # Instructs GridUtils to use panel.pane.Matplotlib for plot objects 
            }
        
        '''
        # Setup links to panel, etc
        appKeys = app.keys()
        if 'messages' in appKeys:
            self.msgBox = app['messages']
            self.msgBox.value = "GridUtils attached to application."
        if 'defaultFigureSize' in appKeys:
            self.plotParameterDefaults['figsize'] = app['defaultFigureSize']
        if 'usePaneMatplotlib' in appKeys:
            self.usePaneMatplotlib = app['usePaneMatplotlib']
        else:
            self.usePaneMatplotlib = False       
        
    def printVerbose(self, msg):
        if self.verboseLevel > 0:
            print(msg)
    
    def setDebugLevel(self, newLevel):
        '''Set a new debug level.

        :param newLevel: debug level to set or update
        :type newLevel: integer
        :return: none
        :rtype: none
        
        .. note::
            Areas of code that typically cause errors have try/except blocks.  Some of these
            have python debugging breakpoints that are active when the debug level is set
            to a positive number.        
        '''
        self.debugLevel = newLevel
    
    def setVerboseLevel(self, newLevel):
        '''Set a new verbose level.

        :param newLevel: verbose level to set or update
        :type newLevel: integer
        :return: none
        :rtype: none
        
        .. note::
            Setting this to a positive number will increase the feedback from this
            module.
        '''
        self.verboseLevel = newLevel
    
    # Grid operations
    
    def clearGrid(self):
        '''Call this when you want to erase the current grid and grid parameters.  This also
        clobbers any current plot parameters.
        Do not call this method between plots of the same grid.'''
        self.grid = {}
        self.grid['_meta'] = {}
        self.grid['_meta']['dimensions'] = {}
        self.clearGridParameters()
        self.resetPlotParameters()
        
    def makeGrid(self):
        '''Using supplied grid parameters, populate a grid in memory.'''
        if self.grid['_meta']['gridParameters']['projection']['name'] == 'LambertConformalConic':
            # Sometimes tilt may not be specified, so use a default of 0.0
            if 'tilt' in self.grid['_meta']['gridParameters'].keys():
                tilt = self.grid['_meta']['gridParameters']['tilt']
            else:
                tilt = 0.0
            lonGrid, latGrid = self.generate_regional_spherical(
                self.grid['_meta']['gridParameters']['projection']['lon_0'], self.grid['_meta']['gridParameters']['dx'],
                self.grid['_meta']['gridParameters']['projection']['lat_0'], self.grid['_meta']['gridParameters']['dy'],
                tilt,
                self.grid['_meta']['gridParameters']['gridResolution'] * self.grid['_meta']['gridParameters']['gridMode']
            )
            self.grid['x'] = lonGrid
            self.grid['y'] = latGrid
            self.grid['_meta']['shape'] = lonGrid.shape
            # This technique seems to return a Lambert Conformal Projection with the following properties
            # This only works if the grid does not overlap a polar point
            # (lat_0 - (dy/2), lat_0 + (dy/2))
            self.grid['_meta']['gridParameters']['projection']['lat_1'] =\
                self.grid['_meta']['gridParameters']['projection']['lat_0'] - (self.grid['_meta']['gridParameters']['dy'] / 2.0)
            self.grid['_meta']['gridParameters']['projection']['lat_2'] =\
                self.grid['_meta']['gridParameters']['projection']['lat_0'] + (self.grid['_meta']['gridParameters']['dy'] / 2.0)
    
    # Original functions provided by Niki Zadeh - Lambert Conformal Conic grids
    # Grid creation and rotation in spherical coordinates
    def mesh_plot(self, lon, lat, lon0=0., lat0=90.):
        """Plot a given mesh with a perspective centered at (lon0,lat0)"""
        f = plt.figure(figsize=(8,8))
        ax = plt.subplot(111, projection=cartopy.crs.NearsidePerspective(central_longitude=lon0, central_latitude=lat0))
        ax.set_global()
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()
        (nj,ni) = lon.shape 
        # plotting verticies
        for i in range(0,ni+1,2):
            ax.plot(lon[:,i], lat[:,i], 'k', transform=cartopy.crs.Geodetic())
        for j in range(0,nj+1,2):
            ax.plot(lon[j,:], lat[j,:], 'k', transform=cartopy.crs.Geodetic())
            
        return f, ax
    
    def rotate_x(self,x,y,z,theta):
        """Rotate vector (x,y,z) by angle theta around x axis."""
        """Returns the rotated components."""
        cost = np.cos(theta)
        sint = np.sin(theta)
        yp   = y*cost - z*sint
        zp   = y*sint + z*cost
        return x,yp,zp
    
    def rotate_y(self,x,y,z,theta):
        """Rotate vector (x,y,z) by angle theta around y axis."""
        """Returns the rotated components."""
        cost = np.cos(theta)
        sint = np.sin(theta)
        zp   = z*cost - x*sint
        xp   = z*sint + x*cost
        return xp,y,zp
    
    def rotate_z(self,x,y,z,theta):
        """Rotate vector (x,y,z) by angle theta around z axis."""
        """Returns the rotated components."""
        cost = np.cos(theta)
        sint = np.sin(theta)
        xp   = x*cost - y*sint
        yp   = x*sint + y*cost
        return xp,yp,z
    
    
    def cart2pol(self,x,y,z):
        """Transform a point on globe from Cartesian (x,y,z) to polar coordinates."""
        """Returns the polar coordinates"""
        lam=np.arctan2(y,x)/self.PI_180
        phi=np.arctan(z/np.sqrt(x**2+y**2))/self.PI_180
        return lam,phi
    
    def pol2cart(self,lam,phi):
        """Transform a point on globe from Polar (lam,phi) to Cartesian coordinates."""
        """Returns the Cartesian coordinates"""
        lam=lam*self.PI_180
        phi=phi*self.PI_180
        x=np.cos(phi)*np.cos(lam)
        y=np.cos(phi)*np.sin(lam)
        z=np.sin(phi)
        return x,y,z
        
    def rotate_z_mesh(self,lam,phi,theta):
        """Rotate the whole mesh on globe by angle theta around z axis (globe polar axis)."""
        """Returns the rotated mesh."""
        #Bring the angle to be in [-pi,pi] so that atan2 would work
        lam       = np.where(lam>180,lam-360,lam)
        #Change to Cartesian coord
        x,y,z     = self.pol2cart(lam,phi)
        #Rotate
        xp,yp,zp  = self.rotate_z(x,y,z,theta)
        #Change back to polar coords using atan2, in [-pi,pi]
        lamp,phip = self.cart2pol(xp,yp,zp)
        #Bring the angle back to be in [0,2*pi]
        lamp      = np.where(lamp<0,lamp+360,lamp)
        return lamp,phip
    
    def rotate_x_mesh(self,lam,phi,theta):
        """Rotate the whole mesh on globe by angle theta around x axis (passing through equator and prime meridian.)."""
        """Returns the rotated mesh."""
        #Bring the angle to be in [-pi,pi] so that atan2 would work
        lam       = np.where(lam>180,lam-360,lam)
        #Change to Cartesian coord
        x,y,z     = self.pol2cart(lam,phi)
        #Rotate
        xp,yp,zp  = self.rotate_x(x,y,z,theta)
        #Change back to polar coords using atan2, in [-pi,pi]
        lamp,phip = self.cart2pol(xp,yp,zp)
        #Bring the angle back to be in [0,2*pi]
        lamp      = np.where(lamp<0,lamp+360,lamp)
        return lamp,phip
    
    def rotate_y_mesh(self,lam,phi,theta):
        """Rotate the whole mesh on globe by angle theta around y axis (passing through equator and prime meridian+90.)."""
        """Returns the rotated mesh."""
        #Bring the angle to be in [-pi,pi] so that atan2 would work
        lam       = np.where(lam>180,lam-360,lam)
        #Change to Cartesian coord
        x,y,z     = self.pol2cart(lam,phi)
        #Rotate
        xp,yp,zp  = self.rotate_y(x,y,z,theta)
        #Change back to polar coords using atan2, in [-pi,pi]
        lamp,phip = self.cart2pol(xp,yp,zp)
        #Bring the angle back to be in [0,2*pi]
        lamp      = np.where(lamp<0,lamp+360,lamp)
        return lamp,phip
    
    def generate_latlon_mesh_centered(self, lni, lnj, llon0, llen_lon, llat0, llen_lat, ensure_nj_even=True):
        """Generate a regular lat-lon grid"""
        self.printVerbose('Generating regular lat-lon grid centered at %.2f %.2f on equator.' % (llon0, llat0))
        llonSP = llon0 - llen_lon/2 + np.arange(lni+1) * llen_lon/float(lni)
        llatSP = llat0 - llen_lat/2 + np.arange(lnj+1) * llen_lat/float(lnj)
        if(llatSP.shape[0]%2 == 0 and ensure_nj_even):
            self.printVerbose("   The number of j's is not even. Fixing this by cutting one row at south.")
            llatSP = np.delete(llatSP,0,0)
        llamSP = np.tile(llonSP,(llatSP.shape[0],1))
        lphiSP = np.tile(llatSP.reshape((llatSP.shape[0],1)),(1,llonSP.shape[0]))
        self.printVerbose('   Generated regular lat-lon grid between latitudes %.2f %.2f' % (lphiSP[0,0],lphiSP[-1,0]))
        self.printVerbose('   Number of js=%d' % (lphiSP.shape[0]))
        #h_i_inv=llen_lon*self.PI_180*np.cos(lphiSP*self.PI_180)/lni
        #h_j_inv=llen_lat*self.PI_180*np.ones(lphiSP.shape)/lnj
        #delsin_j = np.roll(np.sin(lphiSP*self.PI_180),shift=-1,axis=0) - np.sin(lphiSP*self.PI_180)
        #dx_h=h_i_inv[:,:-1]*self._default_Re
        #dy_h=h_j_inv[:-1,:]*self._default_Re
        #area=delsin_j[:-1,:-1]*self._default_Re*self._default_Re*llen_lon*self.self.PI_180/lni
        return llamSP,lphiSP
    
    def generate_regional_spherical(self, lon0, lon_span, lat0, lat_span, tilt, refine):
        """Generate a regional grid centered at (lon0,lat0) with spans of (lon_span,lat_span) and tilted by angle tilt"""
        Ni = int(lon_span*refine)
        Nj = int(lat_span*refine)
       
        #Generate a mesh at equator centered at (lon0, 0)
        lam_,phi_ = self.generate_latlon_mesh_centered(Ni,Nj,lon0,lon_span,0.0,lat_span)
        lam_,phi_ = self.rotate_z_mesh(lam_,phi_, (90.-lon0)*self.PI_180)  #rotate around z to bring it centered at y axis
        lam_,phi_ = self.rotate_y_mesh(lam_,phi_,tilt*self.PI_180)         #rotate around y axis to tilt it as desired
        lam_,phi_ = self.rotate_x_mesh(lam_,phi_,lat0*self.PI_180)         #rotate around x to bring it centered at (lon0,lat0)
        lam_,phi_ = self.rotate_z_mesh(lam_,phi_,-(90.-lon0)*self.PI_180)  #rotate around z to bring it back
                
        return lam_,phi_

    # Grid generation functions
    
    # netCDF operations
    
    def closeDataset(self):
        '''Closes and open dataset file pointer.'''
        if self.ncfp:
            self.ncfp.close()
            self.ncfp = None
            
    def openDataset(self, inputFilename):
        '''Open a netCDF4 grid file.  The file pointer is internal to the object.
        To access it, use: obj.ncfp'''
        # check if we have a vailid inputFilename
        if not(os.path.isfile(inputFilename)):
            self.printVerbose("Dataset not found: %s" % (inputFilename))
            return
                
        # If we have a file pointer and it is open, close it and re-open the new file
        if self.ncfp:
            if self.ncfp.isopen():
                self.ncfp.close()
        
        self.ncfp = nc.Dataset(inputFilename, mode='r')
            
    def readGrid(self, opts={'type': 'MOM6'}):
        '''Read a grid.'''
        '''This can be generalized to work with "other" grids if we desired? (ROMS, HyCOM, etc)'''
        if self.ncfp:
            if opts['type'] == 'MOM6':
                # Load dims
                # Get the grid shape from the dimensions instead of the shape of a variable
                for dimKey in self.ncfp.dimensions:
                    self.grid['_meta']['dimensions'][dimKey] = self.ncfp.dimensions[dimKey].size
                
                self.grid['_meta']['shape'] = (
                    self.grid['_meta']['dimensions']['nyp'], self.grid['_meta']['dimensions']['nxp']
                )
                
                # Load variables ['lons', 'lats']
                read_variables = ['x','y']
                for var in read_variables:
                    self.grid[var] = self.ncfp.variables[var][:][:]

                # Save grid metadata
                self.grid["_meta"]['type'] = opts['type']
                # This method of computing the extent is subject to problems
                self.grid["_meta"]["extent"] = [
                    self.grid['x'].min(), self.grid['x'].max(), self.grid['y'].min(), self.grid['y'].max()
                ]
    
    # Plotting specific functions
    # These functions should not care what grid is loaded. 
    # Plotting is affected by plotParameters and gridParameters.

    def newFigure(self, figsize=(8,6)):
        '''Establish a new matplotlib figure.'''
                     
        fig = Figure(figsize=figsize)
        
        return fig
    
    def plotGrid(self):
        '''Perform a plot operation.
        
        :return: Returns a tuple of matplotlib objects (figure, axes)
        :rtype: tuple
        
        To plot a grid, you first must have the projection set.
        
        :Example:
        
        >>> grd = gridUtils()
        >>> grd.setPlotParameters(
                {
                    ...other grid options...,
                    'projection': {
                        'name': 'Mercator',
                        ...other projection options...,
                    },
        >>> grd.plotGrid()
        '''
        
        if not('shape' in self.grid['_meta'].keys()):
            warnings.warn("Unable to plot the grid.  Missing its 'shape'.")
            return (None, None)
        
        plotProjection = self.getPlotParameter('name', subKey='projection', default=None)
        
        if not(plotProjection):
            warnings.warn("Please set the plot 'projection' parameter 'name'")
            help(self.plotGrid)
            return (None, None)
        
        if plotProjection == 'LambertConformalConic':
            return (self.plotGridLambertConformalConic())
        if plotProjection == 'Mercator':
            return (self.plotGridMercator())
        if plotProjection == 'NearsidePerspective':
            return (self.plotGridNearsidePerspective())
        if plotProjection == 'NorthPolarStereo':
            return (self.plotGridNorthPolarStereo())
        
        warnings.warn("Unable to plot this projection: %s" % (plotProjection))
        return (None, None)

    def plotGridLambertConformalConic(self):
        '''Plot a given mesh using Lambert Conformal Conic projection.'''
        '''Requires: central_latitude, central_longitude and two standard parallels (latitude).'''
        figsize = self.getPlotParameter('figsize', default=(8,8))
        #f = plt.figure(figsize=figsize)
        f = self.newFigure(figsize=figsize)
        central_longitude = self.getPlotParameter('lon_0', subKey='projection', default=-96.0)
        central_latitude = self.getPlotParameter('lat_0', subKey='projection', default=39.0)
        lat_1 = self.getPlotParameter('lat_1', subKey='projection', default=33.0)
        lat_2 = self.getPlotParameter('lat_2', subKey='projection', default=45.0)
        standard_parallels = (lat_1, lat_2)
        crs = cartopy.crs.LambertConformal(
                central_longitude=central_longitude, central_latitude=central_latitude,
                standard_parallels=standard_parallels)
        #ax = plt.subplot(111, 
        #    projection=cartopy.crs.LambertConformal(
        #        central_longitude=central_longitude, central_latitude=central_latitude,
        #        standard_parallels=standard_parallels))
        ax = f.subplots(subplot_kw={'projection': crs})
        mapExtent = self.getPlotParameter('extent', default=[])
        mapCRS = self.getPlotParameter('extentCRS', default=cartopy.crs.PlateCarree())
        if len(mapExtent) == 0:
            ax.set_global()
        else:
            ax.set_extent(mapExtent, crs=mapCRS)
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()
        title = self.getPlotParameter('title', default=None)
        if title:
            ax.set_title(title)
        (nj,ni) = self.grid['_meta']['shape']
        plotAllVertices = self.getPlotParameter('showGridCells', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plot vertices
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.grid['x'][:,i], self.grid['y'][:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.grid['x'][j,:], self.grid['y'][j,:], jColor, linewidth=jLinewidth, transform=transform)
        
        return f, ax

    def plotGridMercator(self):
        '''Plot a given mesh using Mercator projection.'''
        figsize = self.getPlotParameter('figsize', default=(8,8))
        #f = plt.figure(figsize=figsize)
        f = self.newFigure(figsize=figsize)
        central_longitude = self.getPlotParameter('lon_0', subKey='projection', default=0.0)
        crs = cartopy.crs.Mercator(
                central_longitude=central_longitude)
        #ax = plt.subplot(111, 
        #    projection=cartopy.crs.Mercator(
        #        central_longitude=central_longitude))
        ax = f.subplots(subplot_kw={'projection': crs})
        mapExtent = self.getPlotParameter('extent', default=[])
        mapCRS = self.getPlotParameter('extentCRS', default=cartopy.crs.PlateCarree())
        if len(mapExtent) == 0:
            ax.set_global()
        else:
            ax.set_extent(mapExtent, crs=mapCRS)
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()
        title = self.getPlotParameter('title', default=None)
        if title:
            ax.set_title(title)
        (nj,ni) = self.grid['_meta']['shape']
        plotAllVertices = self.getPlotParameter('showGridCells', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plotting vertices
        # For a non conforming projection, we have to plot every line between the points of each grid box
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.grid['x'][:,i], self.grid['y'][:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.grid['x'][j,:], self.grid['y'][j,:], jColor, linewidth=jLinewidth, transform=transform)
                
        return f, ax
    
    def plotGridNearsidePerspective(self):
        """Plot a given mesh using the nearside perspective centered at (central_longitude,central_latitude)"""
        figsize = self.getPlotParameter('figsize', default=(8,8))
        #f = plt.figure(figsize=figsize)
        f = self.newFigure(figsize=figsize)
        central_longitude = self.getPlotParameter('lon_0', subKey='projection', default=0.0)
        central_latitude = self.getPlotParameter('lat_0', subKey='projection', default=90.0)
        satellite_height = self.getPlotParameter('satellite_height', default=35785831)
        crs = cartopy.crs.NearsidePerspective(central_longitude=central_longitude, central_latitude=central_latitude, satellite_height=satellite_height)
        #ax = f.subplots(111, 
        #    projection=cartopy.crs.NearsidePerspective(
        #        central_longitude=central_longitude, central_latitude=central_latitude, satellite_height=satellite_height))
        ax = f.subplots(subplot_kw={'projection': crs})
        if self.usePaneMatplotlib:
            FigureCanvas(f)
        mapExtent = self.getPlotParameter('extent', default=[])
        mapCRS = self.getPlotParameter('extentCRS', default=cartopy.crs.PlateCarree())
        if len(mapExtent) == 0:
            ax.set_global()
        else:
            ax.set_extent(mapExtent, crs=mapCRS)
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()
        title = self.getPlotParameter('title', default=None)
        if title:
            ax.set_title(title)
        (nj,ni) = self.grid['_meta']['shape']
        plotAllVertices = self.getPlotParameter('showGridCells', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plotting vertices
        # For a non conforming projection, we have to plot every line between the points of each grid box
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.grid['x'][:,i], self.grid['y'][:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.grid['x'][j,:], self.grid['y'][j,:], jColor, linewidth=jLinewidth, transform=transform)
                
        return f, ax
        
    def plotGridNorthPolarStereo(self):
        '''Generic plotting function for North Polar Stereo maps'''
        figsize = self.getPlotParameter('figsize', default=(8,8))
        #f = plt.figure(figsize=figsize)
        f = self.newFigure(figsize=figsize)
        central_longitude = self.getPlotParameter('lon_0', subKey='projection', default=0.0)
        true_scale_latitude = self.getPlotParameter('lat_ts', subKey='projection', default=75.0)
        crs = cartopy.crs.NorthPolarStereo(central_longitude=central_longitude, true_scale_latitude=true_scale_latitude)
        #ax = plt.subplot(111, projection=
        #    cartopy.crs.NorthPolarStereo(central_longitude=central_longitude, true_scale_latitude=true_scale_latitude))
        ax = f.subplots(subplot_kw={'projection': crs})
        mapExtent = self.getPlotParameter('extent', default=[])
        mapCRS = self.getPlotParameter('extentCRS', default=cartopy.crs.PlateCarree())
        if len(mapExtent) == 0:
            ax.set_global()
        else:
            ax.set_extent(mapExtent, crs=mapCRS)
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()
        title = self.getPlotParameter('title', default=None)
        if title:
            ax.set_title(title)
        (nj,ni) = self.grid['_meta']['shape']
        plotAllVertices = self.getPlotParameter('showGridCells', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plot vertices
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.grid['x'][:,i], self.grid['y'][:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.grid['x'][j,:], self.grid['y'][j,:], jColor, linewidth=jLinewidth, transform=transform) 
                
        return f, ax
        
    # Grid parameter operations

    def clearGridParameters(self):
        '''Clear grid parameters.  This does not erase any grid data.'''
        self.grid['_meta']['gridParameters'] = {}
        self.grid['_meta']['gridParameterKeys'] = self.grid['_meta']['gridParameters'].keys()
        
    def deleteGridParameters(self, gList, subKey=None):
        """This deletes a given list of grid parameters."""
        
        # Top level subkeys
        if subKey:
            if subKey in self.grid['_meta']['gridParameterKeys']:
                subKeys = self.grid['_meta'][subKey].keys()
                for k in gList:
                    if k in subKeys:
                        self.grid['_meta'][subKey].pop(k, None)
            return
        
        # Top level keys
        for k in gList:
            if k in self.grid['_meta']['gridParameterKeys']:
                self.self.grid['_meta']['gridParameters'].pop(k, None)
                            
        self.grid['_meta']['gridParameterKeys'] = self.grid['_meta']['gridParameters'].keys()

    def getGridParameter(self, gkey, subKey=None, default=None):
        '''Return the requested grid parameter or the default if none is available.'''
        if subKey:
            if subKey in self.grid['_meta']['gridParameterKeys']:
                if gkey in self.grid['_meta']['gridParameters'][subKey].keys():
                    return self.grid['_meta']['gridParameters'][subKey][gkey]
            return default
        
        if gkey in self.grid['_meta']['gridParameterKeys']:
            return self.grid['_meta']['gridParameters'][gkey]
        
        return default
        
    def setGridParameters(self, gridParameters, subKey=None):
        """Generic method for setting gridding parameters using dictionary arguments.
    
        :param gridParameters: grid parameters to set or update
        :type gridParameters: dictionary
        :param subkey: an entry in gridParameters that contains a dictionary of information to set or update
        :type subKey: string
        :return: none
        :rtype: none
        
        .. note::
            Core gridParameter list.  See other grid functions for other potential options.  
            Defaults are marked with an asterisk(*) below.
            
            The gridParameter has a 'projection' subkey that allows 
            
            In general, coordinates are consistent between degrees or meters.  There may
            be some obscure cases where options may be mixed.
            
                'centerUnits': Grid center point units ['degrees'(*), 'meters']
                'east0': Meters east of grid center 
                'north0': Meters north of grid center
                'lon0': Longitude of grid center (may not be the same as the projection center)
                'lat0': Latitude of grid center (may not be the same as the projection center)
                'dx': grid length along x or i axis (generally EW)
                'dy': grid length along y or j axis (generally NS)
                'dxUnits': grid cell units ['degrees'(*), 'meters']
                'dyUnits': grid cell units ['degrees'(*), 'meters']
                'nx': number of grid points along the x or i axis [integer]
                'ny': number of grid points along the y or i axis [integer]
                'tilt': degrees to rotate the grid [float, only available in LambertConformalConic]
                
                SUBKEY: 'projection' (mostly follows proj.org terminology)
                    'name': Grid projection ['LambertConformalConic','Mercator','NorthPolarStereo']
                    'lat_0': Latitude of projection center [degrees, 0.0(*)]
                    'lat_1': First standard parallel (latitude) [degrees, 0.0(*)]
                    'lat_2': Second standard parallel (latitude) [degrees, 0.0(*)]
                    'lat_ts': Latitude of true scale. Defines the latitude where scale is not distorted.
                              Takes precedence over k_0 if both options are used together.
                              For stereographic, if not set, will default to lat_0.
                    'lon_0': Longitude of projection center [degrees, 0.0(*)]
                    'ellps': See proj -le for a list of available ellipsoids [GRS80(*)]
                    'R': Radius of the sphere given in meters.  If both R and ellps are given, R takes precedence.
                    'x_0': False easting (meters, 0.0(*))
                    'y_0': False northing (meters, 0.0(*))
                    'k_0': Depending on projection, this value determines the scale factor for natural origin or the ellipsoid (1.0(*))
                
                MOM6 specific options:
                
                'gridMode': 2 = supergrid(*); 1 = actual grid [integer, 1 or 2(*)]
                'gridResolution': Inverse grid resolution scale factor [float, 1.0(*)]        

            Not to be confused with plotParameters which control how this grid or other
            information is plotted.  For instance, the grid projection and the requested plot
            can be in another projection.
            
        """
        
        # For now pass all keys into the plot parameter dictionary.  Sanity checking is done
        # by the respective makeGrid functions.
        for k in gridParameters.keys():
            if subKey:
                self.grid['_meta']['gridParameters'][subKey][k] = gridParameters[k]
            else:
                self.grid['_meta']['gridParameters'][k] = gridParameters[k]
        
        if not(subKey):
            self.grid['_meta']['gridParameterKeys'] = self.grid['_meta']['gridParameters'].keys()

    def showGridMetadata(self):
        """Show current grid metadata."""
        print(self.grid['_meta'])
            
    def showGridParameters(self):
        """Show current grid parameters."""
        if len(self.grid['_meta']['gridParameterKeys']) > 0:
            print("Current grid parameters:")
            for k in self.grid['_meta']['gridParameterKeys']:
                print("%20s: %s" % (k,self.grid['_meta']['gridParameters'][k]))
        else:
            print("No grid parameters found.")
    
    # Plot parameter operations
        
    def deletePlotParameters(self, pList, subKey=None):
        """This deletes a given list of plot parameters."""
        
        # Top level subkeys
        if subKey:
            if subKey in self.grid['_meta']['plotParameterKeys']:
                subKeys = self.grid['_meta'][subKey].keys()
                for k in pList:
                    if k in subKeys:
                        self.grid['_meta'][subKey].pop(k, None)
            return

        # Top level keys
        for k in pList:
            if k in self.grid['_meta']['plotParameterKeys']:
                self.self.grid['_meta']['plotParameters'].pop(k, None)
                
        self.grid['_meta']['plotParameterKeys'] = self.grid['_meta']['plotParameters'].keys()

    def getPlotParameter(self, pkey, subKey=None, default=None):
        '''Return the requested plot parameter or the default if none is available.
        
           To access dictionary values in projection, use the subKey argument.
        '''
        
        # Top level subkey access
        if subKey:
            if subKey in self.grid['_meta']['plotParameterKeys']:
                try:
                    if pkey in self.grid['_meta']['plotParameters'][subKey].keys():
                        return self.grid['_meta']['plotParameters'][subKey][pkey]
                except:
                    warnings.warn("Attempt to use a subkey(%s) which is not really a subkey? or maybe it should be?" % (subKey))
            return default
        
        # Top level key access
        if pkey in self.grid['_meta']['plotParameterKeys']:
            return self.grid['_meta']['plotParameters'][pkey]
        
        return default

    def resetPlotParameters(self):
        '''Resets plot parameters for a grid.'''
        # Need to use .copy on plotParameterDefaults or we get odd results
        self.grid['_meta']['plotParameters'] = self.plotParameterDefaults.copy()
        self.grid['_meta']['plotParameterKeys'] = self.grid['_meta']['plotParameters'].keys()
    
    def setPlotParameters(self, plotParameters, subKey=None):
        """A generic method for setting plotting parameters using dictionary arguments.

        :param plotParameters: plot parameters to set or update
        :type plotParameters: dictionary
        :param subkey: an entry in plotParameters that contains a dictionary of information to set or update
        :type subKey: string
        :return: none
        :rtype: none
        
        .. note::
            Plot parameters persist for as long as the object exists.
            
            Here is a core list of plot parameters.  Some parameters may be
            grid type specific.
            
                'figsize': tells matplotlib the figure size [width, height in inches (6.4, 4.8)]
                'extent': [x0, x1, y0, y1] map extent of given coordinate system (see extentCRS) [default is []]
                    If no extent is given, [], then set_global() is used. 
                    REF: https://scitools.org.uk/cartopy/docs/latest/matplotlib/geoaxes.html
                'extentCRS': cartopy crs [cartopy.crs.PlateCarree()] 
                    You must have the cartopy.crs module loaded to change the setting.
                'showGrid': show the grid outline [True(*)/False]
                'showGridCells': show the grid cells [True/False(*)]
                'showSupergrid': show the MOM6 supergrid cells [True/False(*)]
                'title': add a title to the plot [None(*)]
                'iColor': matplotlib color for i vertices ['k'(*) black]
                'jColor': matplotlib color for j vertices ['k'(*) black]
                'iLinewidth': matplotlib linewidth for i vertices [points: 1.0(*)]
                'jLinewidth': matplotlib linewidth for j vertices [points: 1.0(*)]
                    For dense gridcells, you can try a very thin linewidth of 0.1.

                SUBKEY: 'projection' (mostly follows proj.org terminology)
                    'name': Grid projection ['LambertConformalConic','Mercator','NorthPolarStereo']
                    'lat_0': Latitude of projection center [degrees, 0.0(*)]
                    'lat_1': First standard parallel (latitude) [degrees, 0.0(*)]
                    'lat_2': Second standard parallel (latitude) [degrees, 0.0(*)]
                    'lat_ts': Latitude of true scale. Defines the latitude where scale is not distorted.
                              Takes precedence over k_0 if both options are used together.
                              For stereographic, if not set, will default to lat_0.
                    'lon_0': Longitude of projection center [degrees, 0.0(*)]
                    'ellps': See proj -le for a list of available ellipsoids [GRS80(*)]
                    'R': Radius of the sphere given in meters.  If both R and ellps are given, R takes precedence.
                    'x_0': False easting (meters, 0.0(*))
                    'y_0': False northing (meters, 0.0(*))
                    'k_0': Depending on projection, this value determines the scale factor for natural origin or the ellipsoid (1.0(*))
                
        """
        
        # For now pass all keys into the plot parameter dictionary.  Sanity checking is done
        # by the respective plotGrid* fuctions.
        for k in plotParameters.keys():
            if subKey:
                self.grid['_meta']['plotParameters'][subKey][k] = plotParameters[k]
            else:
                try:
                    self.grid['_meta']['plotParameters'][k] = plotParameters[k]
                except:
                    if self.debugLevel > 0:
                        pdb.set_trace()
                    raise

        if not(subKey):
            self.grid['_meta']['plotParameterKeys'] = self.grid['_meta']['plotParameters'].keys()

    # Functions from pyroms/examples/grid_MOM6/convert_ROMS_grid_to_MOM6.py
    # Attribution: Mehmet Ilicak via Alistair Adcroft
    # Requires Spherical.py (copied to local lib)
    # Based on code written by Alistair Adcroft and Matthew Harrison of GFDL
    
