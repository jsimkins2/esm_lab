# General imports and definitions
import os, sys
import cartopy
import numpy as np
import xarray as xr
import warnings
import pdb
import logging

# Needed for panel.pane                
from matplotlib.figure import Figure
# not needed for mpl >= 3.1
# Does not cause any problems to continue to use it
from matplotlib.backends.backend_agg import FigureCanvas  

# Required for:
#  * ROMS to MOM6 grid conversion
#  * Computation of MOM6 grid metrics
import spherical

# GridUtils() application
from app import App

class GridUtils:

    def __init__(self, app={}):
        # Constants
        self.PI_180 = np.pi/180.
        self._default_Re = 6.378e6
        
        # File pointer
        self.xrOpen = False
        self.xrFilename = None
        self.xrDS = xr.Dataset()
        self.grid = self.xrDS
        # Internal parameters
        self.usePaneMatplotlib = False
        self.msgBox = None
        # Private variables begin with a _
        # Grid parameters
        self.gridInfo = {}
        self.gridInfo['dimensions'] = {}
        self.gridInfo['gridParameters'] = {}
        self.gridInfo['gridParameterKeys'] = self.gridInfo['gridParameters'].keys()
        # Defaults
        self.plotParameterDefaults = {
            'figsize': (8, 6),
            'extent': [],
            'extentCRS': cartopy.crs.PlateCarree(),
            'projection': {
            },
            'showGrid': True,
            'showGridCells': False,
            'showSupergrid': False
        }          
        # Plot parameters
        self.gridInfo['plotParameters'] = self.plotParameterDefaults
        self.gridInfo['plotParameterKeys'] = self.gridInfo['plotParameters'].keys()
        
        # Messages
        # Logging and Verbosity Levels
        # CRITICAL:50; ERROR:40; WARNING:30; INFO:20, DEBUG:10, NOTSET:0
        self.debugLevel = 0
        self.verboseLevel = logging.INFO
        # Logging options
        self.msgBuffer = []
        self.msgLogger = False
        self.msgLogFile = None
        self.msgLogLevel = logging.WARNING
        self.stdLogValues = {
                0: 'NOTSET',
                10: 'DEBUG',
                20: 'INFO',
                30: 'WARNING',
                40: 'ERROR',
                50: 'CRITICAL'
        }

    # Utility functions

    def addMessage(self, msg):
        '''Append new message to message buffer.'''
        self.msgBuffer.append(msg)
        return

    def app(self):
        '''By calling this function, the user is requesting the application functionality of GridUtils().
           return the dashboard, but GridUtils() also has an internal pointer to the application.'''
        appObj = App(grd=self)
        self.app = appObj
        return appObj.dashboard

    def application(self, app={}):
        '''Convienence function to attach application items to GridUtil so it can update certain portions of the application.
        
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
            msg = "GridUtils application initialized."
            self.printMsg(msg, level=logging.INFO)
        if 'defaultFigureSize' in appKeys:
            self.plotParameterDefaults['figsize'] = app['defaultFigureSize']
        if 'usePaneMatplotlib' in appKeys:
            self.usePaneMatplotlib = app['usePaneMatplotlib']
        else:
            self.usePaneMatplotlib = False

    def adjustExternalLoggers(self):
        '''This adjusts some noisy loggers from other python modules.'''

        # Readjust logger levels to specific noisy modules
        noisyModules = {
                'PIL.PngImagePlugin': logging.ERROR,
                'fiona._env': logging.ERROR,
                'fiona.collection': logging.ERROR,
                'fiona.env': logging.ERROR,
                'fiona.ogrext': logging.ERROR,
                'matplotlib.backends.backend_pdf': logging.ERROR,
                'matplotlib.font_manager': logging.ERROR,
        }

        for moduleName in noisyModules.keys():
            lh = logging.getLogger(moduleName)
            lh.setLevel(noisyModules[moduleName])

        return

    def clearMessage(self):
        '''This clears the message buffer of messages.'''
        self.msgBuffer = []
        return

    def debugMsg(self, msg, level = -1):
        '''This function has a specific purpose to aid in debugging and
        activating pdb breakpoints.  NOTE: pdb breakpoints tend not to
        work very well when running under the application.  It tends to
        terminate the bokeh/tornado server.

        The debug level can be zero(0) and you can forcibly add a break
        point in the code by using `debugMsg(msg, level=3)` anywhere in
        the code.

        .. note::

            Currently defined debug levels:
                0=off 
                1=log debugging messages
                2=raise an exception
                3=stop with a pdb breakpoint after logging message
        '''

        # If level is not set (-1), set to self.debugLevel
        if level == -1:
            level = self.debugLevel

        if level < 1:
            return

        if level >= 1 and msg != "":
            # Send the message at the DEBUG level if 
            # the message is not empty
            self.printMsg(msg, level=logging.DEBUG)

        if level == 2:
            raise

        if level == 3:
            pdb.set_trace()

        return

    def deleteLogfile(self, logFilename):
        '''Delete a log file.  Logging must be off.'''
        if self.msgLogger:
            self.printMsg("Logging active: Unable to delete log file.", level=logging.ERROR)
            return

        if not(os.path.isfile(logFilename)):
            self.printMsg("Logfile (%s) does not exist." % (logFilename), level=logging.ERROR)
            return

        os.unlink(logFilename)
        self.printMsg("Logfile (%s) removed." % (logFilename), level=logging.INFO)

    def disableLogging(self):
        '''Disable logging of messages to a file'''
        self.logHandle.disable = True
        self.logHandle = None
        self.msgLogger = False
        self.printMsg("Logging disabled", level=logging.INFO)

    def enableLogging(self, logFilename):
        '''Enable logging of messages to a file'''

        # Do not permit re-enabling a logging event if logging is already been activated
        if self.msgLogger:
            self.printMsg("Logging already active.  Ignoring request.", level=logging.ERROR)
            return

        # Test to see if we can write to the log file
        success = False
        try:
            fn = open(logFilename, 'a')
            success = True
        except:
            self.printMsg("Failed to open logfile (%s)" % (logFilename), level=logging.CRITICAL)
            return

        if success:
            fn.close()

        logging.basicConfig(filename=logFilename, level=self.msgLogLevel)
        self.adjustExternalLoggers()
        self.logHandle = logging.getLogger(__name__)
        self.logHandle.disabled = False
        self.msgLogger = True
        self.printMsg("Logging enabled", level=logging.INFO)

    def filterLogMessages(self, record):
        '''This may not be needed after all.'''
        print(">>",record.name,__name__)
        if record.name == __name__:
            return True
        return False

    def getDebugLevel(self):
        '''Get the current debug level for GridUtils().  See setDebugLevel() for available
        levels.'''
        return self.debugLevel

    def getLogLevel(self):
        '''Get the current debug level for GridUtils().  See setDebugLevel() for available
        levels.'''
        return self.msgLogLevel

    def getVerboseLevel(self):
        '''Get the current verbose level for GridUtils()'''
        return self.verboseLevel

    def printMsg(self, msg, level = logging.INFO):
        '''
        The verboseLevel and msgLogLevel can be set separately to any level.
        If this is attached to a panel application with a message
        box, the output is sent to that object.  Messages omitting the level argument
        will default to INFO.
        '''

        # Debugging messaging system: may be removed later
        #if self.debugLevel >= 1:
        #    print(">>(%s)(%d)(%d)" % (msg, level, self.verboseLevel))

        # If logging is enabled, send it to the logger
        if self.msgLogger:
            self.logHandle.log(level, msg)

        if level >= self.verboseLevel:
            self.addMessage(msg)

            # Always update the application message box
            # If we don't have a msgBox, then print to STDOUT.
            if hasattr(self, 'msgBox'):
                if self.msgBox:
                    self.msgBox.value = self.showMessages()
                else:
                    print(msg)
            else:
                print(msg)

        return

    def showLoggers(self):
        '''Display an alphabetical list of loggers.  Messages sent at the
        INFO level.'''

        # List of all known loggers: logging.Logger.manager.loggerDict
        loggerNames = list(logging.Logger.manager.loggerDict)
        loggerNames.sort()
        for loggerName in loggerNames:
            logger = logging.getLogger(loggerName)
            msg = "%-40s : %s" % (logger.name, logging.getLevelName(logger.getEffectiveLevel()))
            self.printMsg(msg, level=logging.INFO)

    def showMessages(self):
        '''This converts the message buffer to text with linefeeds.'''
        return '\n'.join(self.msgBuffer) 

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

            Currently defined debug levels:
                0=off 
                1=extra messages
                2=raise an exception
                3=stop at breakpoints
        '''
        self.printMsg("New DEBUG level (%d)" % (newLevel))
        self.debugLevel = newLevel
    
    def setLogLevel(self, newLevel):
        '''Set a new verbose level.

        :param newLevel: verbose level to set or update
        :type newLevel: integer
        :return: none
        :rtype: none
        
        .. note::
            Setting this to a positive number will increase the feedback from this
            module.
        '''
        self.msgLogLevel = newLevel
        # Also update the logger, if active
        if self.msgLogger:
            self.logHandle.setLevel(newLevel)
    
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
        
        # If there are file resources open, close them first.
        if self.xrOpen:
            self.closeDataset()
        
        self.xrFilename = None
        self.xrDS = xr.Dataset()
        self.grid = self.xrDS
        self.gridInfo = {}
        self.gridInfo['dimensions'] = {}
        self.clearGridParameters()
        self.resetPlotParameters()
        
    def computeGridMetrics(self):
        '''Compute MOM6 grid metrics: angle_dx, dx, dy and area.'''

        self.grid.attrs['grid_version'] = "0.2"
        self.grid.attrs['code_version'] = "GridTools: beta"
        self.grid.attrs['history'] = "sometime: GridTools"
        self.grid.attrs['projection'] = self.gridInfo['gridParameters']['projection']['name']
        self.grid.attrs['proj'] = self.gridInfo['gridParameters']['projection']['proj']
        
        R = 6370.e3 # Radius of sphere        

        # Make a copy of the lon grid as values are changed for computation
        lon = self.grid.x.copy()
        lat = self.grid.y
        
        # Approximate edge lengths as great arcs
        self.grid['dx'] = (('nyp', 'nx'),  R * spherical.angle_through_center( (lat[ :,1:],lon[ :,1:]), (lat[:  ,:-1],lon[:  ,:-1]) ))
        self.grid.dx.attrs['units'] = 'meters'
        self.grid['dy'] = (('ny' , 'nxp'), R * spherical.angle_through_center( (lat[1:, :],lon[1:, :]), (lat[:-1,:  ],lon[:-1,:  ]) ))
        self.grid.dy.attrs['units'] = 'meters'
        
        # Scaling by latitude?
        cos_lat = np.cos(np.radians(lat))
        
        # Presize the angle_dx array
        angle_dx = np.zeros(lat.shape)
        # Fix lon so they are 0 to 360 for computation of angle_dx
        lon = np.where(lon < 0., lon+360, lon)
        angle_dx[:,1:-1] = np.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
        angle_dx[:, 0  ] = np.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
        angle_dx[:,-1  ] = np.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
        self.grid['angle_dx'] = (('nyp', 'nxp'), angle_dx)
        self.grid.angle_dx.attrs['units'] = 'radians'
        
        self.grid['area'] = (('ny','nx'), R * R * spherical.quad_area(lat, lon))
        self.grid.area.attrs['units'] = 'meters^2'

        return
   



    def makeGrid(self):
        '''Using supplied grid parameters, populate a grid in memory.'''

        # Make a grid in the Mercator projection
        if self.gridInfo['gridParameters']['projection']['name'] == "Mercator":      
            if 'tilt' in self.gridInfo['gridParameters'].keys():
                tilt = self.gridInfo['gridParameters']['tilt']
            else:
                tilt = 0.0

            lonGrid, latGrid = self.generate_regional_spherical(
                self.gridInfo['gridParameters']['projection']['lon_0'], self.gridInfo['gridParameters']['dx'],
                0, # lat0 is 0 in mercator projection?
                self.gridInfo['gridParameters']['dy'],
                tilt,
                self.gridInfo['gridParameters']['gridResolution'] * self.gridInfo['gridParameters']['gridMode']
            )

            (nxp, nyp) = lonGrid.shape

            self.grid['x'] = (('nyp','nxp'), lonGrid)
            self.grid.x.attrs['units'] = 'degrees_east'
            self.grid['y'] = (('nyp','nxp'), latGrid)
            self.grid.y.attrs['units'] = 'degrees_north'

            # This technique seems to return a Lambert Conformal Projection with the following properties
            # This only works if the grid does not overlap a polar point
            # (lat_0 - (dy/2), lat_0 + (dy/2))
            #self.gridInfo['gridParameters']['projection']['lat_1'] =\
               # self.gridInfo['gridParameters']['projection']['lat_0'] - (self.gridInfo['gridParameters']['dy'] / 2.0)
            #self.gridInfo['gridParameters']['projection']['lat_2'] =\
              #  self.gridInfo['gridParameters']['projection']['lat_0'] + (self.gridInfo['gridParameters']['dy'] / 2.0)
            self.gridInfo['gridParameters']['projection']['proj'] =\
                    "+ellps=WGS84 +proj=merc +lon_0=%s +x_0=0.0 +y_0=0.0 +units=m +no_defs" %\
                        (self.gridInfo['gridParameters']['projection']['lon_0'])

            # Declare the xarray dataset open even though it is really only in memory at this point
            self.xrOpen = True

            # Compute grid metrics
            self.computeGridMetrics()


        # Make a grid in the North Polar Stereo projection
        if self.gridInfo['gridParameters']['projection']['name'] == "NorthPolarStereo":
            if 'tilt' in self.gridInfo['gridParameters'].keys():
                tilt = self.gridInfo['gridParameters']['tilt']
            else:
                tilt = 0.0

            lonGrid, latGrid = self.generate_regional_spherical(
                self.gridInfo['gridParameters']['projection']['lon_0'], self.gridInfo['gridParameters']['dx'],
                self.gridInfo['gridParameters']['projection']['lat_0'], self.gridInfo['gridParameters']['dy'],
                tilt,
                self.gridInfo['gridParameters']['gridResolution'] * self.gridInfo['gridParameters']['gridMode']
            )

            (nxp, nyp) = lonGrid.shape

            self.grid['x'] = (('nyp','nxp'), lonGrid)
            self.grid.x.attrs['units'] = 'degrees_east'
            self.grid['y'] = (('nyp','nxp'), latGrid)
            self.grid.y.attrs['units'] = 'degrees_north'

            # This technique seems to return a Lambert Conformal Projection with the following properties
            # This only works if the grid does not overlap a polar point
            # (lat_0 - (dy/2), lat_0 + (dy/2))
            #self.gridInfo['gridParameters']['projection']['lat_1'] =\
                #self.gridInfo['gridParameters']['projection']['lat_0'] - (self.gridInfo['gridParameters']['dy'] / 2.0)
            #self.gridInfo['gridParameters']['projection']['lat_2'] =\
               # self.gridInfo['gridParameters']['projection']['lat_0'] + (self.gridInfo['gridParameters']['dy'] / 2.0)
            self.gridInfo['gridParameters']['projection']['proj'] =\
                    "+ellps=WGS84 +proj=stere +lat_0=%s +lon_0=%s  +x_0=0.0 +y_0=0.0 +no_defs" %\
                        (self.gridInfo['gridParameters']['projection']['lat_0'],
                        self.gridInfo['gridParameters']['projection']['lon_0'])

            # Declare the xarray dataset open even though it is really only in memory at this point
            self.xrOpen = True

            # Compute grid metrics
            self.computeGridMetrics()


        # Make a grid in the South Polar Stereo projection
        if self.gridInfo['gridParameters']['projection']['name'] == "SouthPolarStereo":
            if 'tilt' in self.gridInfo['gridParameters'].keys():
                tilt = self.gridInfo['gridParameters']['tilt']
            else:
                tilt = 0.0

            lonGrid, latGrid = self.generate_regional_spherical(
                self.gridInfo['gridParameters']['projection']['lon_0'], self.gridInfo['gridParameters']['dx'],
                self.gridInfo['gridParameters']['projection']['lat_0'], self.gridInfo['gridParameters']['dy'],
                tilt,
                self.gridInfo['gridParameters']['gridResolution'] * self.gridInfo['gridParameters']['gridMode']
            )

            (nxp, nyp) = lonGrid.shape

            self.grid['x'] = (('nyp','nxp'), lonGrid)
            self.grid.x.attrs['units'] = 'degrees_east'
            self.grid['y'] = (('nyp','nxp'), latGrid)
            self.grid.y.attrs['units'] = 'degrees_north'

            # This technique seems to return a Lambert Conformal Projection with the following properties
            # This only works if the grid does not overlap a polar point
            # (lat_0 - (dy/2), lat_0 + (dy/2))
            #self.gridInfo['gridParameters']['projection']['lat_1'] =\
                #self.gridInfo['gridParameters']['projection']['lat_0'] - (self.gridInfo['gridParameters']['dy'] / 2.0)
            #self.gridInfo['gridParameters']['projection']['lat_2'] =\
                #self.gridInfo['gridParameters']['projection']['lat_0'] + (self.gridInfo['gridParameters']['dy'] / 2.0)
            self.gridInfo['gridParameters']['projection']['proj'] =\
                    "+ellps=WGS84 +proj=stere +lat_0=%s +lon_0=%s  +x_0=0.0 +y_0=0.0 +no_defs" %\
                        (self.gridInfo['gridParameters']['projection']['lat_0'],
                        self.gridInfo['gridParameters']['projection']['lon_0'])

            # Declare the xarray dataset open even though it is really only in memory at this point
            self.xrOpen = True

            # Compute grid metrics
            self.computeGridMetrics()


        # Make a grid in the Lambert Conformal Conic projection
        if self.gridInfo['gridParameters']['projection']['name'] == 'LambertConformalConic':
            # Sometimes tilt may not be specified, so use a default of 0.0
            if 'tilt' in self.gridInfo['gridParameters'].keys():
                tilt = self.gridInfo['gridParameters']['tilt']
            else:
                tilt = 0.0

            lonGrid, latGrid = self.generate_regional_spherical(
                self.gridInfo['gridParameters']['projection']['lon_0'], self.gridInfo['gridParameters']['dx'],
                self.gridInfo['gridParameters']['projection']['lat_0'], self.gridInfo['gridParameters']['dy'],
                tilt,
                self.gridInfo['gridParameters']['gridResolution'] * self.gridInfo['gridParameters']['gridMode']
            )

            (nxp, nyp) = lonGrid.shape

            self.grid['x'] = (('nyp','nxp'), lonGrid)
            self.grid.x.attrs['units'] = 'degrees_east'
            self.grid['y'] = (('nyp','nxp'), latGrid)
            self.grid.y.attrs['units'] = 'degrees_north'

            # This technique seems to return a Lambert Conformal Projection with the following properties
            # This only works if the grid does not overlap a polar point
            # (lat_0 - (dy/2), lat_0 + (dy/2))
            self.gridInfo['gridParameters']['projection']['lat_1'] =\
                self.gridInfo['gridParameters']['projection']['lat_0'] - (self.gridInfo['gridParameters']['dy'] / 2.0)
            self.gridInfo['gridParameters']['projection']['lat_2'] =\
                self.gridInfo['gridParameters']['projection']['lat_0'] + (self.gridInfo['gridParameters']['dy'] / 2.0)
            self.gridInfo['gridParameters']['projection']['proj'] =\
                    "+ellps=WGS84 +proj=lcc +lon_0=%s +lat_0=%s +x_0=0.0 +y_0=0.0 +lat_1=%s +lat_2=%s +no_defs" %\
                        (self.gridInfo['gridParameters']['projection']['lat_0'],
                        self.gridInfo['gridParameters']['projection']['lon_0'],
                        self.gridInfo['gridParameters']['projection']['lat_1'],
                        self.gridInfo['gridParameters']['projection']['lat_2'])

            # Declare the xarray dataset open even though it is really only in memory at this point
            self.xrOpen = True

            # Compute grid metrics
            self.computeGridMetrics()
       

    
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
        msg = 'Generating regular lat-lon grid centered at %.2f %.2f on equator.' % (llon0, llat0)
        self.printMsg(msg, level=logging.INFO)
        llonSP = llon0 - llen_lon/2 + np.arange(lni+1) * llen_lon/float(lni)
        llatSP = llat0 - llen_lat/2 + np.arange(lnj+1) * llen_lat/float(lnj)
        if(llatSP.shape[0]%2 == 0 and ensure_nj_even):
            msg = "   The number of j's is not even. Fixing this by cutting one row at south."
            self.printMsg(msg, level=logging.INFO)
            llatSP = np.delete(llatSP,0,0)
        llamSP = np.tile(llonSP,(llatSP.shape[0],1))
        lphiSP = np.tile(llatSP.reshape((llatSP.shape[0],1)),(1,llonSP.shape[0]))
        msg = '   Generated regular lat-lon grid between latitudes %.2f %.2f' % (lphiSP[0,0],lphiSP[-1,0])
        self.printMsg(msg, level=logging.INFO)
        msg = '   Number of js=%d' % (lphiSP.shape[0])
        self.printMsg(msg, level=logging.INFO)
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
    
    # xarray Dataset operations
    
    def closeDataset(self):
        '''Closes and open dataset file pointer.'''
        if self.xrOpen:
            self.xrDS.close()
            self.xrOpen = False
            
    def openDataset(self, inputFilename):
        '''Open a grid file.  The file pointer is internal to the object.
        To access it, use: obj.xrDS or obj.grid'''
        # check if we have a vailid inputFilename
        if not(os.path.isfile(inputFilename)):
            self.printMsg("Dataset not found: %s" % (inputFilename), level=logging.INFO)
            return
                
        # If we have a file pointer and it is open, close it and re-open the new file
        if self.xrOpen:
            self.closeDataset()
            
        try:
            self.xrDS = xr.open_dataset(inputFilename)
            self.xrOpen = True
            self.xrFilename = inputFilename
        except:
            msg = "ERROR: Unable to load dataset: %s" % (inputFilename)
            self.printMsg(msg, level=logging.ERROR)
            self.xrDS = None
            self.xrOpen = False
            # Stop on error to load a file
            self.debugMsg("")
            
    def readGrid(self, opts={'type': 'MOM6'}, local=None, localFilename=None):
        '''Read a grid.
        
        This can be generalized to work with "other" grids if we desired? (ROMS, HyCOM, etc)
        '''
        # if a dataset is being loaded via readGrid(local=), close any existing dataset
        if local:
            if self.xrOpen:
                self.closeDataset()
            self.xrOpen = True
            self.xrDS = local
            self.grid = local
        else:
            if self.xrOpen:
                if opts['type'] == 'MOM6':
                    # Save grid metadata
                    self.gridInfo['type'] = opts['type']
                    self.grid = self.xrDS
        
        if localFilename:
            self.xrFilename = localFilename

    def removeFillValueAttributes(self):

        ncEncoding = {}
        ncVars = list(self.grid.variables)
        for ncVar in ncVars:
            ncEncoding[ncVar] = {'_FillValue': None}

        return ncEncoding

    
    def saveGrid(self, filename=None):
        '''
        This operation is destructive using the last known filename which can be overridden.
        '''
        if filename:
            self.xrFilename = filename
            
        try:
            self.grid.to_netcdf(self.xrFilename, encoding=self.removeFillValueAttributes())
            msg = "Successfully wrote netCDF file to %s" % (self.xrFilename)
            self.printMsg(msg, level=logging.INFO)
        except:
            msg = "Failed to write netCDF file to %s" % (self.xrFilename)
            self.printMsg(msg, level=logging.INFO)
    
    # Plotting specific functions
    # These functions should not care what grid is loaded. 
    # Plotting is affected by plotParameters and gridParameters.

    def newFigure(self, figsize=None):
        '''Establish a new matplotlib figure.'''
        
        if figsize:
            figsize = self.getPlotParameter('figsize', default=figsize)             
        else:
            figsize = self.getPlotParameter('figsize', default=self.plotParameterDefaults['figsize'])
            
        fig = Figure(figsize=figsize)
        
        return fig
    
    # insert plotGrid.ipynb here
    
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

        #if not('shape' in self.gridInfo.keys()):
        #    warnings.warn("Unable to plot the grid.  Missing its 'shape'.")
        #    return (None, None)

        plotProjection = self.getPlotParameter('name', subKey='projection', default=None)

        if not(plotProjection):
            msg = "Please set the plot 'projection' parameter 'name'"
            self.printMsg(msg, level=logging.ERROR)
            #warnings.warn("Please set the plot 'projection' parameter 'name'")
            return (None, None)

        # initiate new plot, infer projection within the plotting procedure
        f = self.newFigure()

        # declare projection options - note that each projection uses a different
        # combination of these parameters and rarely all are used for one projection
        central_longitude = self.getPlotParameter('lon_0', subKey='projection', default=0.0)
        central_latitude = self.getPlotParameter('lat_0', subKey='projection', default=39.0)
        lat_1 = self.getPlotParameter('lat_1', subKey='projection', default=33.0)
        lat_2 = self.getPlotParameter('lat_2', subKey='projection', default=45.0)
        standard_parallels = (lat_1, lat_2)
        satellite_height = self.getPlotParameter('satellite_height', default=35785831)
        true_scale_latitude = self.getPlotParameter('lat_ts', subKey='projection', default=75.0)

        # declare varying crs based on plotProjection
        crs = None
        if plotProjection == 'LambertConformalConic':
            crs = cartopy.crs.LambertConformal(
            central_longitude=central_longitude, central_latitude=central_latitude,
            standard_parallels=standard_parallels)
        if plotProjection == 'Mercator':
            crs = cartopy.crs.Mercator(central_longitude=central_longitude)
        if plotProjection == 'NearsidePerspective':
            crs = cartopy.crs.NearsidePerspective(central_longitude=central_longitude,
                central_latitude=central_latitude, satellite_height=satellite_height)
        if plotProjection == 'NorthPolarStereo':
            crs = cartopy.crs.NorthPolarStereo(central_longitude=central_longitude,
                true_scale_latitude=true_scale_latitude)
        if plotProjection == 'SouthPolarStereo':
            crs = cartopy.crs.SouthPolarStereo(central_longitude=central_longitude,
                true_scale_latitude=true_scale_latitude)

        if crs == None:
            #warnings.warn("Unable to plot this projection: %s" % (plotProjection))
            msg = "Unable to plot this projection: %s" % (plotProjection)
            self.printMsg(msg, level=logging.ERROR)
            return (None, None)

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
        nj = self.grid.dims['nyp']
        ni = self.grid.dims['nxp']
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

    # Grid parameter operations

    def clearGridParameters(self):
        '''Clear grid parameters.  This does not erase any grid data.'''
        self.gridInfo['gridParameters'] = {}
        self.gridInfo['gridParameterKeys'] = self.gridInfo['gridParameters'].keys()

    def deleteGridParameters(self, gList, subKey=None):
        """This deletes a given list of grid parameters."""

        # Top level subkeys
        if subKey:
            if subKey in self.gridInfo['gridParameterKeys']:
                subKeys = self.gridInfo[subKey].keys()
                for k in gList:
                    if k in subKeys:
                        self.gridInfo[subKey].pop(k, None)
            return

        # Top level keys
        for k in gList:
            if k in self.gridInfo['gridParameterKeys']:
                self.self.gridInfo['gridParameters'].pop(k, None)
                            
        self.gridInfo['gridParameterKeys'] = self.gridInfo['gridParameters'].keys()

    def getGridParameter(self, gkey, subKey=None, default=None):
        '''Return the requested grid parameter or the default if none is available.'''
        if subKey:
            if subKey in self.gridInfo['gridParameterKeys']:
                if gkey in self.gridInfo['gridParameters'][subKey].keys():
                    return self.gridInfo['gridParameters'][subKey][gkey]
            return default
        
        if gkey in self.gridInfo['gridParameterKeys']:
            return self.gridInfo['gridParameters'][gkey]
        
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
                self.gridInfo['gridParameters'][subKey][k] = gridParameters[k]
            else:
                self.gridInfo['gridParameters'][k] = gridParameters[k]
        
        if not(subKey):
            self.gridInfo['gridParameterKeys'] = self.gridInfo['gridParameters'].keys()

    def showGridMetadata(self):
        """Show current grid metadata."""
        print(self.gridInfo)
            
    def showGridParameters(self):
        """Show current grid parameters."""
        if len(self.gridInfo['gridParameterKeys']) > 0:
            self.printMsg("Current grid parameters:", level=logging.INFO)
            for k in self.gridInfo['gridParameterKeys']:
                self.printMsg("%20s: %s" % (k,self.gridInfo['gridParameters'][k]), level=logging.INFO)
        else:
            self.printMsg("No grid parameters found.", level=logging.ERROR)
    
    # Plot parameter operations
        
    def deletePlotParameters(self, pList, subKey=None):
        """This deletes a given list of plot parameters."""
        
        # Top level subkeys
        if subKey:
            if subKey in self.gridInfo['plotParameterKeys']:
                subKeys = self.gridInfo[subKey].keys()
                for k in pList:
                    if k in subKeys:
                        self.gridInfo[subKey].pop(k, None)
            return

        # Top level keys
        for k in pList:
            if k in self.gridInfo['plotParameterKeys']:
                self.self.gridInfo['plotParameters'].pop(k, None)
                
        self.gridInfo['plotParameterKeys'] = self.gridInfo['plotParameters'].keys()

    def getPlotParameter(self, pkey, subKey=None, default=None):
        '''Return the requested plot parameter or the default if none is available.
        
           To access dictionary values in projection, use the subKey argument.
        '''
        
        # Top level subkey access
        if subKey:
            if subKey in self.gridInfo['plotParameterKeys']:
                try:
                    if pkey in self.gridInfo['plotParameters'][subKey].keys():
                        return self.gridInfo['plotParameters'][subKey][pkey]
                except:
                    msg = "Attempt to use a subkey(%s) which is not really a subkey? or maybe it should be?" % (subKey)
                    self.printMsg(msg, level=logging.WARNING)
            return default
        
        # Top level key access
        if pkey in self.gridInfo['plotParameterKeys']:
            return self.gridInfo['plotParameters'][pkey]
        
        return default

    def resetPlotParameters(self):
        '''Resets plot parameters for a grid.'''
        # Need to use .copy on plotParameterDefaults or we get odd results
        self.gridInfo['plotParameters'] = self.plotParameterDefaults.copy()
        self.gridInfo['plotParameterKeys'] = self.gridInfo['plotParameters'].keys()
    
    def showPlotParameters(self):
        """Show current plot parameters."""
        if len(self.gridInfo['plotParameterKeys']) > 0:
            self.printMsg("Current plot parameters:", level=logging.INFO)
            for k in self.gridInfo['plotParameterKeys']:
                self.printMsg("%20s: %s" % (k,self.gridInfo['plotParameters'][k]), level=logging.INFO)
        else:
            self.printMsg("No plot parameters found.", level=logging.INFO)
    
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
            try:
                if subKey:
                    self.gridInfo['plotParameters'][subKey][k] = plotParameters[k]
                else:
                    self.gridInfo['plotParameters'][k] = plotParameters[k]
            except:
                msg = 'Failed to assign a plotParameter(%s) with "%s"' %\
                        (k, plotParameters[k])
                if subKey:
                    msg = '%s with subkey "%s"' % (msg, subKey)
                self.debugMsg(msg)

        if not(subKey):
            self.gridInfo['plotParameterKeys'] = self.gridInfo['plotParameters'].keys()

    # Functions from pyroms/examples/grid_MOM6/convert_ROMS_grid_to_MOM6.py
    # Attribution: Mehmet Ilicak via Alistair Adcroft
    # Requires spherical.py (copied to local lib)
    # Based on code written by Alistair Adcroft and Matthew Harrison of GFDL
    
