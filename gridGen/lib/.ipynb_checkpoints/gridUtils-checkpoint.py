# From: 
#  https://github.com/nikizadehgfdl/grid_generation/blob/dev/jupynotebooks/regional_grid_spherical.ipynb

# TODO:
#   * We can place default settings into __init__ and spread through the code.
#     this would allow users to change defaults globally.
#   * Add a formal logging mechanism.
#   * Bring in code that converts ROMS grids to MOM6 grids
#     * Allow conversion of MOM6 grids to ROMS grids
#   * writeGrid()

#General imports and definitions
import os, sys
import cartopy
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc

# For use with python notebooks
#%matplotlib inline

class gridUtils:

    def __init__(self):
        # Constants
        self.PI_180 = np.pi/180.
        self._default_Re = 6.378e6
        # File pointer
        self.ncfp = None
        # Debugging/logging
        self.debugLevel = 0
        self.verboseLevel = 0
        # Grid variables
        self.gridX = None
        self.gridY = None
        # matplotlib plot extent
        # degrees: [minLon, maxLon, minLat, maxLat]
        self.gridExtent = []
        self.gridParameters = {}
        self.gridShape = None
        self.resetGridParameters()
        # Plotting controls
        self.plotParameters = {}
        self.plotParameterKeys = []

    # Utility functions
        
    def printVerbose(self, msg):
        if self.verboseLevel > 0:
            print(msg)        
    
    # Grid operations
    
    def resetGrid(self):
        self.gridX = None
        self.gridY = None
        self.gridShape = None
        self.gridExtent = []

    def makeGrid(self):
        '''Using supplied grid parameters, populate a grid in memory.'''
        if self.gridParameters['gridProjection'] == 'Mercator':
            #lamc, phic = grd.generate_regional_spherical(lon0, lon_span, lat0, lat_span, tilt, refineR*refineS)
            lonGrid, latGrid = self.generate_regional_spherical(
                self.gridParameters['lonGridCenter'], self.gridParameters['lonSpan'],
                self.gridParameters['latGridCenter'], self.gridParameters['latSpan'],
                self.gridParameters['gridTilt'],
                self.gridParameters['nominalResolution'] * self.gridParameters['nominalSpacing']
            )
            #grd.mesh_plot(lamc, phic, lon0, lat0)
            self.gridX = lonGrid
            self.gridY = latGrid
            self.gridShape = lonGrid.shape
    
    # Original functions provided by Niki Zadeh
    # Grid creation and rotation in spherical coordinates
    def mesh_plot(self, lon, lat, lon0=0., lat0=90.):
        """Plot a given mesh with a perspective centered at (lon0,lat0)"""
        plt.figure(figsize=(8,8))
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
        self.printVerbose('Generating regular lat-lon grid between centered at %.2f %.2f' % (llon0, llat0))
        llonSP = llon0 - llen_lon/2 + np.arange(lni+1) * llen_lon/float(lni)
        llatSP = llat0 - llen_lat/2 + np.arange(lnj+1) * llen_lat/float(lnj)
        if(llatSP.shape[0]%2 == 0 and ensure_nj_even):
            self.printVerbose("   The number of j's is not even. Fixing this by cutting one row at south.")
            llatSP = np.delete(llatSP,0,0)
        llamSP = np.tile(llonSP,(llatSP.shape[0],1))
        lphiSP = np.tile(llatSP.reshape((llatSP.shape[0],1)),(1,llonSP.shape[0]))
        self.printVerbose('   generated regular lat-lon grid between latitudes %.2f %.2f' % (lphiSP[0,0],lphiSP[-1,0]))
        self.printVerbose('   number of js=%d' % (lphiSP.shape[0]))
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
    
    # MOM6 notes
    # nx, ny : grid centers
    # nxp, nxp : grid verticies
    
    # ROMS
    # has an extra grid box around the regular grid
    
    # netCDF operations
    
    def closeDataset(self):
        if self.ncfp:
            self.ncfp.close()
            self.ncfp = None
            
    def openDataset(self, inputFilename):
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
                # lons
                self.gridX = self.ncfp.variables['x'][:][:]
                # lats
                self.gridY = self.ncfp.variables['y'][:][:]
                # nj, ni (nyp, nxp)
                self.gridShape = self.gridX.shape
                self.gridExtent = [self.gridX.min(), self.gridX.max(), self.gridY.min(), self.gridY.max()]
    
    # Plotting specific functions
    # These functions should not care what grid is loaded.
    
    def plotGrid(self):
        '''Perform a plot grid operation.'''
        if not('view' in self.plotParameterKeys):
            self.printVerbose("Please set the 'view' plot parameter before plotting.")
            self.printVerbose("view = ['NearsidePerspective', 'NorthPolarStereo', 'LambertConformalConic', 'Mercator']")
            return
        
        if self.plotParameters['view'] == 'LambertConformalConic':
            self.plotGridLambertConformalConic()
        if self.plotParameters['view'] == 'Mercator':
            self.plotGridMercator()
        if self.plotParameters['view'] == 'NearsidePerspective':
            self.plotGridNearsidePerspective()
        if self.plotParameters['view'] == 'NorthPolarStereo':
            self.plotGridNorthPolarStereo()

    def plotGridLambertConformalConic(self):
        '''Plot a given mesh using Lambert Conformal Conic projection.'''
        '''Requires: central_latitude, central_longitude and two standard parallels (latitude).'''
        figsize = self.getPlotParameter('figsize', default=(8,8))
        plt.figure(figsize=figsize)
        central_longitude = self.getPlotParameter('central_longitude', default=-96.0)
        central_latitude = self.getPlotParameter('central_latitude', default=39.0)
        standard_parallels = self.getPlotParameter('standard_parallels', default=(33.0, 45.0))
        ax = plt.subplot(111, 
            projection=cartopy.crs.LambertConformal(
                central_longitude=central_longitude, central_latitude=central_latitude,
                standard_parallels=standard_parallels))
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
        (nj,ni) = self.gridShape
        plotAllVertices = self.getPlotParameter('showGrid', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plot vertices
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.gridX[:,i], self.gridY[:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.gridX[j,:], self.gridY[j,:], jColor, linewidth=jLinewidth, transform=transform) 

    def plotGridMercator(self):
        '''Plot a given mesh using Mercator projection.'''
        figsize = self.getPlotParameter('figsize', default=(8,8))
        plt.figure(figsize=figsize)
        central_longitude = self.getPlotParameter('central_longitude', default=0.0)
        ax = plt.subplot(111, 
            projection=cartopy.crs.Mercator(
                central_longitude=central_longitude))
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
        (nj,ni) = self.gridShape
        plotAllVertices = self.getPlotParameter('showGrid', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plotting vertices
        # For a non conforming projection, we have to plot every line between the points of each grid box
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.gridX[:,i], self.gridY[:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.gridX[j,:], self.gridY[j,:], jColor, linewidth=jLinewidth, transform=transform)
    
    def plotGridNearsidePerspective(self):
        """Plot a given mesh using the nearside perspective centered at (central_longitude,central_latitude)"""
        figsize = self.getPlotParameter('figsize', default=(8,8))
        plt.figure(figsize=figsize)
        central_longitude = self.getPlotParameter('central_longitude', default=0.0)
        central_latitude = self.getPlotParameter('central_latitude', default=90.0)
        satellite_height = self.getPlotParameter('satellite_height', default=35785831)
        ax = plt.subplot(111, 
            projection=cartopy.crs.NearsidePerspective(
                central_longitude=central_longitude, central_latitude=central_latitude, satellite_height=satellite_height))
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
        (nj,ni) = self.gridShape
        plotAllVertices = self.getPlotParameter('showGrid', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plotting vertices
        # For a non conforming projection, we have to plot every line between the points of each grid box
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.gridX[:,i], self.gridY[:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.gridX[j,:], self.gridY[j,:], jColor, linewidth=jLinewidth, transform=transform)
        
    def plotGridNorthPolarStereo(self):
        '''Generic plotting function for North Polar Stereo maps'''
        figsize = self.getPlotParameter('figsize', default=(8,8))
        plt.figure(figsize=figsize)
        central_longitude = self.getPlotParameter('central_longitude', default=0.0)
        true_scale_latitude = self.getPlotParameter('true_scale_latitude', default=75.0)
        ax = plt.subplot(111, projection=cartopy.crs.NorthPolarStereo(central_longitude=central_longitude, true_scale_latitude=true_scale_latitude))
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
        (nj,ni) = self.gridShape
        plotAllVertices = self.getPlotParameter('showGrid', default=False)
        iColor = self.getPlotParameter('iColor', default='k')
        jColor = self.getPlotParameter('jColor', default='k')
        transform = self.getPlotParameter('transform', default=cartopy.crs.Geodetic())
        iLinewidth = self.getPlotParameter('iLinewidth', default=1.0)
        jLinewidth = self.getPlotParameter('jLinewidth', default=1.0)
        
        # plot vertices
        for i in range(0,ni+1,2):
            if (i == 0 or i == (ni-1)) or plotAllVertices:
                ax.plot(self.gridX[:,i], self.gridY[:,i], iColor, linewidth=iLinewidth, transform=transform)
        for j in range(0,nj+1,2):
            if (j == 0 or j == (nj-1)) or plotAllVertices:
                ax.plot(self.gridX[j,:], self.gridY[j,:], jColor, linewidth=jLinewidth, transform=transform)                
        
    # Grid parameter operations
    
    def deleteGridParameters(self, pList):
        """This deletes a given list of grid parameters."""
        
        for k in pList:
            if k in self.gridParameterKeys:
                self.gridParameters.pop(k, None)
                
        self.gridParameterKeys = self.gridParameters.keys()

    def getGridParameter(self, pkey, default=None):
        '''Return the requested grid parameter or the default if none is available.'''
        if pkey in self.gridParameterKeys:
            return self.gridParameters[pkey]
        
        return default

    def resetGridParameters(self):
        self.gridParameters = {
            'lonGridCenter': None,
            'latGridCenter': None,
            'lonSpan': None,
            'latSpan': None,
            'gridTilt': 0.0,
            
        }
        
    def setGridParameters(self, pDict):
        """A generic method for setting gridding parameters using dictionary arguments."""
        """The information is passed into a dictionary that is persistent within the object."""
        
        # For now pass all keys into the plot parameter dictionary.  Sanity checking is done
        # by the respective makeGrid functions.
        for k in pDict.keys():
            self.gridParameters[k] = pDict[k]
            
        self.gridParameterKeys = self.gridParameters.keys()    
    
    # Plot parameter operations
        
    def deletePlotParameters(self, pList):
        """This deletes a given list of plot parameters."""
        
        for k in pList:
            if k in self.plotParameterKeys:
                self.plotParameters.pop(k, None)
                
        self.plotParameterKeys = self.plotParameters.keys()

    def getPlotParameter(self, pkey, default=None):
        '''Return the requested plot parameter or the default if none is available.'''
        if pkey in self.plotParameterKeys:
            return self.plotParameters[pkey]
        
        return default

    def resetPlotParameters(self):
        self.plotParameters = {}
        self.plotParameterKeys = []
    
    def setPlotParameters(self, pDict):
        """A generic method for setting plotting parameters using dictionary arguments."""
        """The information is passed into a dictionary that is persistent within the object."""
        
        # For now pass all keys into the plot parameter dictionary.  Sanity checking is done
        # by the respective plotGrid* fuctions.
        for k in pDict.keys():
            self.plotParameters[k] = pDict[k]
            
        self.plotParameterKeys = self.plotParameters.keys()