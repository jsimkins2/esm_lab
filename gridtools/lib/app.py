# GridUtils.App()

# Modules

import os, sys, io, logging
import numpy as np
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import netCDF4 as nc
import warnings
import xarray as xr
import xgcm
from io import BytesIO
import panel as pn
pn.extension()

# This is called by GridTools() and can't be
# called by itself.

class App:

    def __init__(self, grd=None):
        # Globals
        
        # This application has its own copy of GridTools() object
        self.grd = grd

        # Default filenames
        self.defaultGridFilename = 'gridFile.nc'
        self.defaultLogFilename = 'logFile.log'

        # How we grow the grid from the specified latitude (lat0) or longitude (lon0)
        # TODO: If 'Center' is not chosen then the controls for Central Latitude and Central Longitude no longer make sense.
        # For now: we assume Center/Center for making grids.
        self.xGridModes = ['Left', 'Center', 'Right']
        self.yGridModes = ['Lower', 'Center', 'Upper']

        # For now only MOM6 works.  It could work for other grids!
        #gridTypes = ["MOM6", "ROMS", "WRF"]
        self.gridTypes = ["MOM6"]

        # Generic true/false indicators
        self.trueFalseNames = ["False","True"]
        self.trueFalseValues = [False, True]
        self.trueFalseDict = dict(zip(self.trueFalseNames, self.trueFalseValues))

        # Log levels
        self.logLevelNames = ["NOTSET", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        self.logLevelValues = [0, 10, 20, 30, 40, 50]
        self.logLevelDict = dict(zip(self.logLevelNames, self.logLevelValues))

        # Debug levels
        self.debugLevelNames = ["OFF", "MESSAGE", "RAISE", "BREAKPOINT"]
        self.debugLevelValues = [0, 1, 2, 3]
        self.debugLevelDict = dict(zip(self.debugLevelNames, self.debugLevelValues))

        # Available plot projections
        self.plotProjections = ["Nearside Perspective", "Mercator", "Lambert Conformal Conic", "North Polar Stereographic", "South Polar Stereographic"]
        self.plotProjectionsGridTools = ["NearsidePerspective", "Mercator", "LambertConformalConic", "NorthPolarStereo", "SouthPolarStereo"]
        self.plotProjectionsDict = dict(zip(self.plotProjections, self.plotProjectionsGridTools))

        # Supported grid projections
        self.projNames = ["Mercator", "Lambert Conformal Conic", "North Polar Stereographic", "South Polar Stereographic"]
        self.gridToolNames = ["Mercator", "LambertConformalConic", "NorthPolarStereo", "SouthPolarStereo"]
        self.projNamesGridTools = dict(zip(self.projNames, self.gridToolNames))
        self.projCarto = [ccrs.Mercator(), ccrs.LambertConformal(), ccrs.NorthPolarStereo(), ccrs.SouthPolarStereo()]
        self.projDict = dict(zip(self.projNames, self.projCarto))

        # Plot grid modes
        # gridExtent, gridCells, superGrid
        self.plotGridModes = ['gridExtent', 'gridCells', 'superGrid']
        self.plotGridModesDescriptions = ['Grid Extent', 'Grid Cells', 'Supergrid Cells']
        self.plotGridModeDict = dict(zip(self.plotGridModesDescriptions, self.plotGridModes))

        self.plotColors = ['b', 'c', 'g', 'k', 'm', 'r', 'y']
        self.plotColorsDescriptions = ['Blue', 'Cyan', 'Green', 'Black', 'Magenta', 'Red', 'Yellow']
        self.plotColorDict = dict(zip(self.plotColorsDescriptions, self.plotColors))

        self.plotLineStyles = ['solid', 'dotted', 'dashed', 'dashdot']
        self.plotLineStylesDescriptions = ['Solid', 'Dotted', 'Dashed', 'DashDot']
        self.plotLineStyleDict = dict(zip(self.plotLineStylesDescriptions, self.plotLineStyles))

        # This controls the default figure size of the plot in the panel application
        # TODO: Improve integration
        # aspect 4:3, default dpi=144
        self.widthIn = 5.0
        self.heightIn = (self.widthIn * 3.0) / 4.0
        self.defaultPlotFigureSize = (self.widthIn, self.heightIn)

        # plotWidgetWidth and plotWidgetHeight
        # NOTE: These values are used by other controls
        self.plotWidgetWidth = 800
        self.plotWidgetHeight = 600

        # Setup for other internals of GridTools()
        self.useNumpyPi = False
        self.enableLogging = False
        self.loggingFilename = None
        self.verboseLevel = 0
        self.debugLevel = 0
        
        self.initializeWidgets()
        self.initializeTabs()
        self.initializeDashboard()

    # Panel application functions

    def clearInformationWindow(self, event):
        self.grd.clearMessage()
        self.statusWidget.value = self.grd.showMessages()
        return

    def updateDataView(self):
        self.dataView[0] = self.grd.grid
        return

    def updateFilename(self, newFilename):
        self.gridFilenameLocal.value = newFilename
        self.gridFilenameRemote.value = newFilename
        return

    def downloadNetCDF(self):
        # See if we have a race condition
        self.saveLocalGridButton.filename = self.gridFilenameLocal.value
        bout = self.grd.grid.to_netcdf(encoding=self.grd.removeFillValueAttributes())
        bio = BytesIO()
        bio.write(bout)
        bio.seek(0)
        return bio

    def loadLocalGrid(self, event):
        if self.localFileSelection.value == None:
            msg = "A grid file has not been selected in the Local File tab."
            self.grd.printMsg(msg, logging.INFO)
            return

        # Test to see if xarray can load the selected file
        try:
            ncTest = xr.load_dataset(localFileSelection.value)
            msg = "The grid file %s was loaded." % (self.localFileSelection.filename)
            self.grd.printMsg(msg, logging.INFO)
            self.grd.clearGrid()
            self.grd.readGrid(local=ncTest, localFilename=self.localFileSelection.filename)
            self.updateDataView()
            self.updateFilename(self.localFileSelection.filename)
        except:
            msg = "The grid file %s was not loadable." % (self.localFileSelection.filename)
            self.grd.printMsg(msg, logging.ERROR)

        return
    
    def loadRemoteGrid(self, event):
        ct = len(self.remoteFileSelection.value)

        if ct == 0:
            msg = "A grid file has not been selected in the Remote File tab."
            self.grd.printMsg(msg, logging.INFO)
            return

        try:
            fileToOpen = self.remoteFileSelection.value[0]
            self.grd.openDataset(self.remoteFileSelection.value[0])
            self.grd.readGrid()
            msg = "The grid file %s was loaded." % (fileToOpen)
            self.grd.printMsg(msg, logging.INFO)
            self.updateDataView()
            self.updateFilename(fileToOpen)
        except:
            msg = "Failed to load grid file: %s" % (fileToOpen)
            self.grd.printMsg(msg, logging.ERROR)

        return

    def saveRemoteGrid(self, event):
        '''Attempt to save grid to remote filesystem using last known grid filename.'''
        self.grd.saveGrid(filename=self.gridFilenameRemote.value)

    def make_grid(self, event):
        updateMessage = "Nothing happened."
        msg = "Running make_grid()"
        self.grd.printMsg(msg, logging.INFO)
        self.grd.clearGrid()
        self.grd.setGridParameters({
            'projection': {
                'name': self.projNamesGridTools[self.gridProjection.value],
                'lon_0': float(self.glon0.value),
                'lat_0': float(self.glat0.value)
            },
            'dx': int(self.dx.value),
            'dy': int(self.dy.value),
            'dxUnits': 'degrees',
            'dyUnits': 'degrees',    
            'gridResolution': float(self.gridResolution.value),
            'gridMode': float(self.gridMode.value),
            'tilt': float(self.gtilt.value)
        })
        self.grd.makeGrid()

        # Update the plot if we updated the grid
        self.plotWindow.object = self.make_plot()

        # Update grid info
        self.updateDataView()

        if self.projNamesGridTools[self.gridProjection.value] == 'LambertConformalConic':
            # For this projection LCC sets lat_1 and lat_2 based on grid inputs.
            updateMessage = "NOTICE: Grid first and second parallels (lat_1, lat_2) have been changed to (%s, %s)." %\
                (self.grd.gridInfo['gridParameters']['projection']['lat_1'], self.grd.gridInfo['gridParameters']['projection']['lat_2'])
            self.glat1.value = self.grd.gridInfo['gridParameters']['projection']['lat_1']
            self.glat2.value = self.grd.gridInfo['gridParameters']['projection']['lat_2']

        msg = "Make grid succeeded: %s" % (updateMessage)
        self.grd.printMsg(msg, logging.INFO)

        return

    def make_plot(self):
        msg = "Running make_plot()"
        self.grd.printMsg(msg, logging.INFO)

        if self.plotTitle.value != "":
            mp_title = plotTitle.value
        else:
            selectedProjection = self.plotProjection.value
            mp_title = "%s: " % (selectedProjection) + str(self.dx.value) + "x" + str(self.dy.value) + " with " + str(self.gtilt.value) + " degree tilt"

        # Check plotGridMode.value to set plot parameter showGridCells
        showGridCellsState = False
        pGridMode = self.plotGridModeDict[self.plotGridMode.value]
        if pGridMode == 'gridCells':
            showGridCellsState = True

        # Determine plot extent (this may vary depending on selected projection)
        plotExtentState = []
        # If we are not using the global projection, use the user supplied extents
        if not(self.plotUseGlobal.value):
            x0pt = self.plotExtentX0.value
            x1pt = self.plotExtentX1.value
            y0pt = self.plotExtentY0.value
            y1pt = self.plotExtentY1.value

            plotExtentState = [x0pt, x1pt, y0pt, y1pt]

        # These inputs will have to change based on selected projection
        self.grd.setPlotParameters(
            {
                'figsize': self.defaultPlotFigureSize,
                'projection' : {
                    'name': self.plotProjectionsDict[self.plotProjection.value],
                    'lat_0': float(self.plat0.value),
                    'lon_0': float(self.plon0.value)
                },
                'extent': plotExtentState,
                'iLinewidth': self.plotXLineWidth.value,
                'jLinewidth': self.plotYLineWidth.value,
                'showGridCells': showGridCellsState,
                'title': mp_title,
                'iColor': self.plotColorDict[self.plotXColor.value],
                'jColor': self.plotColorDict[self.plotYColor.value]
            }
        )
        if self.grd.xrOpen:
            (figure, axes) = self.grd.plotGrid()
            msg = "Running make_plot(): done"     
            self.grd.printMsg(msg, logging.INFO)
        else:
            (figure, axes) = self.errorFigure()
            msg = "Running make_plot(): plotting failure"
            self.grd.printMsg(msg, logging.ERROR)

        return figure

    
    def errorFigure(self):
        '''Create Blank Plot to signal plotting failure. This signals a problem within the user specifications in relation to code capabilities.  '''
        
        f = self.grd.newFigure()
        central_longitude = self.grd.getPlotParameter('lon_0', subKey='projection', default=0.0)
        central_latitude = self.grd.getPlotParameter('lat_0', subKey='projection', default=90.0)
        satellite_height = self.grd.getPlotParameter('satellite_height', default=35785831)
        crs = cartopy.crs.NearsidePerspective(central_longitude=central_longitude, central_latitude=central_latitude, satellite_height=satellite_height)
        ax = f.subplots(subplot_kw={'projection': crs})
        if self.grd.usePaneMatplotlib:
            FigureCanvas(f)
        mapExtent = self.grd.getPlotParameter('extent', default=[])
        mapCRS = self.grd.getPlotParameter('extentCRS', default=cartopy.crs.PlateCarree())
        ax.set_global()
        ax.coastlines()
        ax.gridlines()
        ax.set_title("Plot Failure", color='red')
        ax.text(0.5, 0.4, 'please check plot/grid parameters and retry', transform=ax.transAxes,
                fontsize=10, color='red', alpha=0.2,
                ha='center', va='center', rotation='0')
        return f, ax
    
    def initializePlot(self):
        ''' Plot the initial image upon loading up the application. This is developed to differentiate between plot failure and the first plot.'''
        f = self.grd.newFigure()
        satellite_height = self.grd.getPlotParameter('satellite_height', default=35785831)
        crs = cartopy.crs.NearsidePerspective(central_longitude=290, central_latitude=30, satellite_height=satellite_height)
        ax = f.subplots(subplot_kw={'projection': crs})
        if self.grd.usePaneMatplotlib:
            FigureCanvas(f)
        mapExtent = self.grd.getPlotParameter('extent', default=[])
        mapCRS = self.grd.getPlotParameter('extentCRS', default=cartopy.crs.PlateCarree())
        ax.set_global()
        ax.stock_img()
        ax.coastlines()
        ax.gridlines()
        ax.set_title("Welcome! Please specify your grid and plot parameters")

        return f
    
    def plotRefresh(self, event):
        self.plotWindow.object = self.make_plot()
        return

    def showManual(self):
        manualTabs = pn.Tabs()

        pageMain = pn.WidgetBox('''
        # Instructions
        This will be the eventual location for the instruction manual for this application.  Information
        found here will mostly pertain to the operation of this application.  Additional details about
        the MOM6 model can be found in the [MOM6 User Manual](https://mom6.readthedocs.io/){target="_blank"}.
        
        ''', width=self.plotWidgetWidth)

        pageGrids = pn.WidgetBox('''
        # Grids
        ***We love grids!***
        ## Grid Reference
        This controls how the grid is grown from the selected latitude (lat0) and longitude (lon0) using
        degrees or x (x0) or y (y0) using meters.  By default, the grid is grown from the center point
        in both directions based on the size (dy, dx) and grid resolution.  For now, dy and dx can only
        operate in degrees.  In the future, grids may be build with other fixed points of reference.

        ## Grid Type
        For now, only MOM6 is supported.  Other grid types may be possible in the future.

        ## Grid Mode
        Internally, this mode is 2 which really means computations
        are done to compute vertices for the grid cells and vertices through the center points of the
        grid cells.  At present, this mode should not be anything other than 2 for MOM6 grids.

        ## Grid Representation
        Here is a representation of a (2, 3) MOM6 grid adapted from convert_ROMS_grid_to_MOM6.py
        by Mehmet Ilicak and Alistair Adcroft.  NOTE: The MOM6 supergrid is (5, 7) in shape.

        ```text
          G SG
             5 + | + | + | +
          2  4 - p - p - p -
             3 + | + | + | +
          1  2 - p - p - p -
             1 + | + | + | +
                 1   2   3    G
               1 2 3 4 5 6 7  SG

        KEY: p = ROMS rho (center) points; also MOM6 h (center) points
             + = ROMS psi (corner) points
             - = ROMS u points
             | = ROMS v points
             G = grid points
            SG = supergrid points
        ```

        A MOM6 grid of (ny, nx) will have (ny\*2+1, nx\*2+1) points on the supergrid.
        NOTE: In python, array storage is zero based.  In the above example, supergrid 
        point (1, 1) is stored in memory location (0, 0).
        ''', width=self.plotWidgetWidth)

        pageLogging = pn.WidgetBox('''
        # Logging
        The logging mechanism in this application and GridUtils() is slightly complex.  For messages
        emitted to the "Information" panel or using iterative means, you can control the amount of 
        detail presented to you or logged in a file.  The logging levels from low to high are: NOTSET, 
        DEBUG, INFO, WARNING, ERROR and CRITICAL.  The level set means only messages of that level or higher
        will be shown or logged.  If you want to see all available detail, use NOTSET.  NOTE: The detail sent
        to the "Information" window by default is INFO or higher.  The detail sent to a log file, if enabled,
        is WARNING or higher.  The function for emitting messages is `GridUtils.printMsg()`.
        # Log file
        You can only change the log file name or delete the log file when the logging system is disabled.
        To assist the software developers, we request that you provide a log file of activity to help us
        discover problems with the code.  The log file will continue to grow over time.  It is a good idea
        to periodically erase the log file.
        # Debug level
        This is a special feature mainly for developers.  If you are planning to "hack" this code, you can
        utilize this feature to assist with debugging existing or new code.  The available debug levels
        do not operate like the logging levels.  For operational use, the debug level is usually OFF.  You
        can use the MESSAGE level to simply emit messsages for debugging.  The debug level RAISE, can emit a
        message and then raise a python exception.  This can normally be done in a try/except block where you
        can try a bit of code and in the except block raise the exception after emitting a debugging message.
        The last level is BREAKPOINT.  This is similar to RAISE except that after the message is emitted, the
        program will attempt to start the python debugger (pdb) using `pdb.set_trace()`.  All messages sent
        via `GridUtils.debugMsg()` are shown at the DEBUG level.

        **NOTE**: Setting breakpoints do not work very well in the application.  The application
        is running a server.  When a breakpoint is triggered, it will crash the server running the application.
        ''', width=self.plotWidgetWidth)

        pageNumpyPi = pn.WidgetBox('''
        # numpypi
        Activating numpypi will replace some mathematical routines in numpy with slower
        routines that will produce bitwise identical results.  For more information on
        this package, please [consult this webpage](https://github.com/adcroft/numpypi){target="_blank"}.
        NOTE: The numpypi module provides portable intrinsic functions that return the
        same bitwise floating-point values on different platforms.  Not all numpy 
        routines are replaced.
        ''', width=self.plotWidgetWidth)

        manualTabs.extend([
            ('Main', pageMain),
            ('Grids', pageGrids),
            ('Logging', pageLogging),
            ('Numpypi', pageNumpyPi),
        ])

        return manualTabs
    
    def initializeWidgets(self):
        # Widgets
        
        # The text area input box can show many lines and automatically adds a scroll bar for a long message
        self.statusWidget = pn.widgets.TextAreaInput(name='Information', value="", background="skyblue", height=100,
                width=self.plotWidgetWidth+100)
        self.clearLogButton = pn.widgets.Button(name='Clear Information', button_type='primary', height=50, width=125)
        self.clearLogButton.on_click(self.clearInformationWindow)

        # Grid Controls
        # Use: Niki's defaults for rapid testing
        # 30x20 tilt 30 deg lat_0 40.0 lon_0 230.0 Res 1.0
        self.gridProjection = pn.widgets.Select(name='Projection', options=self.projNames, value=self.projNames[1])
        self.gridType = pn.widgets.Select(name="Grid Type", options=self.gridTypes, value=self.gridTypes[0])
        self.gridType.disabled = True
        self.gridResolution = pn.widgets.Spinner(name="Grid Resolution", value=1.0, step=0.1, start=0.0, end=10.0, width=80)
        self.gridMode = pn.widgets.Spinner(name="Grid Mode", value=2, step=1, start=1, end=2, width=80)
        self.gridMode.disabled = True
        self.unitNames = ['degrees','meters']
        self.dxdyUnits = pn.widgets.Select(name='Units', options=self.unitNames, value=self.unitNames[0])
        self.dx = pn.widgets.Spinner(name="dx", value=20, step=1, start=0, end=100, width=100)
        self.dy = pn.widgets.Spinner(name="dy", value=30, step=1, start=0, end=100, width=100)
        self.glon0 = pn.widgets.Spinner(name="Central Longitude(lon_0) (0 to 360)", value=230.0, step=1.0, start=0.0, end=360.0, width=100)
        self.glat0 = pn.widgets.Spinner(name="Central Latitude(lat_0) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.glat1 = pn.widgets.Spinner(name="First Parallel(lat_1) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.glat2 = pn.widgets.Spinner(name="Second Parallel(lat_2) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.glatts = pn.widgets.Spinner(name="Latitude of True Scale(lat_ts) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.gtilt = pn.widgets.Spinner(name="Tilt (-90 to 90)", value=30.0, step=0.1, start=-90.0, end=90.0, width=100)
        self.gridControlUpdateButton = pn.widgets.Button(name='Make Grid', button_type='primary')
        self.gridControlUpdateButton.on_click(self.make_grid)
        self.xGridControl = pn.widgets.Select(name='X grid mode', options=self.xGridModes, value=self.xGridModes[1])
        self.yGridControl = pn.widgets.Select(name='Y grid mode', options=self.yGridModes, value=self.yGridModes[1])

        # Plot Controls
        # Use Niki's defaults for rapid testing
        # extent: -160, -100, 60, 20

        # Projection
        self.plotProjection = pn.widgets.Select(name='Projection', options=self.plotProjections, value=self.plotProjections[0])
        self.plon0 = pn.widgets.Spinner(name="Central Longitude(lon_0) (0 to 360)", value=230.0, step=1.0, start=0.0, end=360.0, width=100)
        self.plat0 = pn.widgets.Spinner(name="Central Latitude(lat_0) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.plat1 = pn.widgets.Spinner(name="First Parallel(lat_1) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.plat2 = pn.widgets.Spinner(name="Second Parallel(lat_2) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.platts = pn.widgets.Spinner(name="Latitude of True Scale(lat_ts) (-90 to 90)", value=40.0, step=1.0, start=-90.0, end=90.0, width=100)

        # Extent
        # CARTOPY: (x0, x1, y0, y1)
        #  https://scitools.org.uk/cartopy/docs/latest/matplotlib/geoaxes.html
        self.plotExtentX0 = pn.widgets.Spinner(name="Longitude(x0) (-180 to 180)", value=-160.0, step=1.0, start=-180.0, end=180.0, width=100)
        self.plotExtentX1 = pn.widgets.Spinner(name="Longitude(x1) (-180 to 180)", value=-100.0, step=1.0, start=-180.0, end=180.0, width=100)
        self.plotExtentY0 = pn.widgets.Spinner(name="Latitude(y0) (-90 to 90)", value=20.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.plotExtentY1 = pn.widgets.Spinner(name="Latitude(y1) (-90 to 90)", value=60.0, step=1.0, start=-90.0, end=90.0, width=100)
        self.plotUseGlobal = pn.widgets.Checkbox(name="Use global extent (disables custom extent)")

        # Style
        self.plotTitle = pn.widgets.TextInput(name='Plot title', value="", width=250)
        self.plotGridMode = pn.widgets.Select(name='Grid Style', options=self.plotGridModesDescriptions, value=self.plotGridModesDescriptions[1])
        self.plotXColor = pn.widgets.Select(name='x Color', options=self.plotColorsDescriptions, value=self.plotColorsDescriptions[3])
        self.plotYColor = pn.widgets.Select(name='y Color', options=self.plotColorsDescriptions, value=self.plotColorsDescriptions[3])
        self.plotXLineWidth = pn.widgets.Spinner(name="x Line Width", value=1.0, step=0.1, start=0.01, end=10.0, width=80)
        self.plotYLineWidth = pn.widgets.Spinner(name="y Line Width", value=1.0, step=0.1, start=0.01, end=10.0, width=80)

        # Grid Save/Load controls

        # Grid file name
        # NOTE: Sharing a text field is not recommended.  It will display
        # more than once, but only one will actually update.
        self.gridFilenameLocal = pn.widgets.TextInput(name='Grid filename', value=self.defaultGridFilename, width=200)
        self.gridFilenameRemote = pn.widgets.TextInput(name='Grid filename', value=self.defaultGridFilename, width=200)

        # File download button and call back function
        # We can't put a variable in the filename= argument below.  It isn't updated
        # when the assigned variable is updated.  Any updates need to be done to
        # saveLocalGridButton.filename when the local file name is changed.  We
        # also discovered if we updated in the callback it also works.  Might
        # encounter a race condition later.  Watch for it.
        self.saveLocalGridButton = pn.widgets.FileDownload(
            label="Download Grid",
            button_type='success',
            callback=self.downloadNetCDF,
            filename=self.defaultGridFilename)

        # Local file selection
        self.localFileSelection = pn.widgets.FileInput(accept='.nc')

        # Remote file selection
        self.remoteFileSelection = pn.widgets.FileSelector('~', file_pattern='*.nc')

        # Load grid buttons
        self.loadLocalGridButton = pn.widgets.Button(name='Load Local Grid', button_type='primary')
        self.loadLocalGridButton.on_click(self.loadLocalGrid)

        self.loadRemoteGridButton = pn.widgets.Button(name='Load Remote Grid', button_type='primary')
        self.loadRemoteGridButton.on_click(self.loadRemoteGrid)

        self.saveRemoteGridButton = pn.widgets.Button(name='Save Remote Grid', button_type='success')
        self.saveRemoteGridButton.on_click(self.saveRemoteGrid)

        # Plot controls
        self.plotControlUpdateButton = pn.widgets.Button(name='Plot', button_type='primary')
        self.plotControlUpdateButton.on_click(self.plotRefresh)

        # The plot itself wrapped in a widget
        # Use panel.pane.Matplotlib(matplotlib.figure)
        self.plotWindow = pn.pane.Matplotlib(self.initializePlot(), width=self.plotWidgetWidth, height=self.plotWidgetHeight)

        # This presents a data view summary of the xarray object
        self.dataView = pn.Column(self.grd.grid, width=self.plotWidgetWidth)

        # Setup controls
        self.logEnableControl = pn.widgets.Checkbox(name="Enable file logging")
        self.logEnableControl.param.watch(self.logEnableCallback, 'value')
        #self.grd.debugMsg('breakpoint',3)
        self.logFilenameControl = pn.widgets.TextInput(name='Log filename', value=self.defaultLogFilename, width=200)
        self.logLevelControl = pn.widgets.Select(name='Log level', options=self.logLevelNames, value=self.logLevelNames[3])
        self.logLevelControl.param.watch(self.logLevelCallback, 'value')
        self.logEraseButton = pn.widgets.Button(name='Erase log file', button_type='danger', height=50, width=200)
        self.logEraseButton.on_click(self.deleteLogfile)
        self.informationLevelControl = pn.widgets.Select(name='Information level', options=self.logLevelNames, value=self.logLevelNames[2])
        self.debugLevelControl = pn.widgets.Select(name='Debug level', options=self.debugLevelNames, value=self.debugLevelNames[0])
        self.debugLevelControl.param.watch(self.debugLevelCallback, 'value')
        self.enableNumpyPiControl = pn.widgets.Checkbox(name="Enable numpypi bitwise-the-same")

    def deleteLogfile(self, event):
        '''This function is called as a result of pushing the "Erase logfile" button in the application.
        This places a call into GridTools.deleteLogfile(filename).'''
        self.grd.deleteLogfile(self.logFilenameControl.value)
        return

    def logEnableCallback(self, event):
        msg = "logEnableCallback event"
        self.grd.printMsg(msg, logging.DEBUG)
        if hasattr(event, 'name') and hasattr(event, 'old') and hasattr(event, 'new') and hasattr(event, 'type'):
            if event.name == 'value' and event.type == 'changed':
                if event.new == True:
                    self.grd.enableLogging(self.logFilenameControl.value)
                if event.new == False:
                    self.grd.disableLogging()
        else:
            msg = "Illegal event passed to App.logEnableCallback()"
            self.grd.printMsg(msg, logging.ERROR)

        return

    def logLevelCallback(self, event):
        msg = "logLevelCallback event"
        self.grd.printMsg(msg, logging.DEBUG)
        if hasattr(event, 'name') and hasattr(event, 'old') and hasattr(event, 'new') and hasattr(event, 'type'):
            if event.name == 'value' and event.type == 'changed':
                self.grd.setLogLevel(self.logLevelDict[event.new])
        
    def debugLevelCallback(self, event):
        msg = "debugLevelCallback event"
        self.grd.printMsg(msg, logging.DEBUG)
        if hasattr(event, 'name') and hasattr(event, 'old') and hasattr(event, 'new') and hasattr(event, 'type'):
            if event.name == 'value' and event.type == 'changed':
                self.grd.setDebugLevel(self.debugLevelDict[event.new])
        
    def initializeTabs(self):
        # Tabs

        # Plot, Grid and Setup controls
        self.controlTabs = pn.Tabs()
        self.plotControlTabs = pn.Tabs()
        self.gridControlTabs = pn.Tabs()
        self.setupControlTabs = pn.Tabs()

        # If the Alt layout works, we can replace the existing.
        self.displayTabs = pn.Tabs()
        self.saveLoadTabs = pn.Tabs()

        # Pull controls together

        # Plot controls
        self.plotProjectionControls = pn.WidgetBox('# Plot Projection', self.plotProjection, self.plon0, self.plat0, self.plat1, self.plat2, self.platts, self.plotControlUpdateButton)
        self.plotExtentControls = pn.WidgetBox('# Plot Extent', self.plotExtentX0, self.plotExtentX1, self.plotExtentY0, self.plotExtentY1, self.plotUseGlobal)
        self.plotStyleControls = pn.WidgetBox('# Plot Style', self.plotTitle, self.plotGridMode, self.plotXColor, self.plotYColor, self.plotXLineWidth, self.plotYLineWidth)

        # Grid controls
        self.gridProjectionControls = pn.WidgetBox('# Grid Projection', self.gridProjection, self.glon0, self.glat0, self.glat1, self.glat2, self.glatts, self.gtilt, self.gridControlUpdateButton)
        self.gridSpacingControls = pn.WidgetBox('# Grid Spacing', self.dx, self.dy, self.gridResolution, self.dxdyUnits, self.gridControlUpdateButton)
        self.gridAdvancedControls = pn.WidgetBox(
            """
            See "Grids" Manual tab for details about these controls.
            ## Grid Reference
            """, self.xGridControl, self.yGridControl, """    
            ## Grid Type
            For now, only MOM6 grids are supported.
            """, self.gridType, """
            ## Grid Mode
            For now, MOM6 grids require grid mode 2.
            """, self.gridMode)

        # Setup controls
        self.loggingControls = pn.WidgetBox('# Logging','''
        These controls allow you to limit the level of output put into the Information window.  Logging to
        an external file is also available.  See the Manual tab called "Logging" for more information.
        ''',
                self.logFilenameControl,
                self.logEnableControl,
                self.logEraseButton,
                self.logLevelControl,
                self.informationLevelControl,
                self.debugLevelControl
        )

        self.numpyPiControls = pn.WidgetBox('''
        # numpypi
        See Manual tab "Numpypi" for more information.

        NOTE: This control does not do anything yet!
        ''',
                self.enableNumpyPiControl)

        # Place controls into respective tabs

        # Control hierarchy (left panel)
        # Plot
        #  Projection
        #  Extent
        #  Style
        # Grid
        #  Projection
        #  Spacing
        #  Advanced
        # Setup
        #  Logging

        # Top level
        self.controlTabs.extend([
            ('Plot', self.plotControlTabs),
            ('Grid', self.gridControlTabs),
            ('Setup', self.setupControlTabs)
        ])

        # Plot
        self.plotControlTabs.extend([
            ('Projection', self.plotProjectionControls),
            ('Extent', self.plotExtentControls),
            ('Style', self.plotStyleControls)
        ])

        # Grid
        self.gridControlTabs.extend([
            ('Projection', self.gridProjectionControls),
            ('Spacing', self.gridSpacingControls),
            ('Advanced', self.gridAdvancedControls)
        ])

        self.localFilesWindow = pn.WidgetBox(
            '''
            # Local Files
            If you are running this notebook on the same computer as your web browser, accessing files
            from the Local Files tab and the Remote Files tab should look the same.  If you are running
            this notebook on a remote system, you may need to use the Remote Files tab to load grids on
            the remote system.  There may be size limit for loading/downloading files via the web
            browser (Local Files). You may change the grid filename prior to saving the grid.  Do not
            use it for file selection.  At this time, we only accept NetCDF file formats. 
            ''', '''### Upload Grid''', self.localFileSelection, self.loadLocalGridButton,
            ''' ### Download Grid''', self.gridFilenameLocal, self.saveLocalGridButton)

        self.remoteFilesWindow = pn.WidgetBox(self.loadRemoteGridButton, self.gridFilenameRemote, self.saveRemoteGridButton,
            '''
            # Remote Files
            This tab loads and saves grids to the remote system.  If you are running this notebook
            on the same system, either file tab will work.  You may change the grid filename prior
            to saving the grid.  Do not use it for file selection.
            ''', self.remoteFileSelection)

        self.saveLoadTabs.extend([
            ('Local Files', self.localFilesWindow),
            ('Remote Files', self.remoteFilesWindow)
        ])

        # Plotting area, data view, local/remote file access and manual
        self.displayTabs.extend([
            ('Grid Plot', self.plotWindow),
            ('Grid Info', self.dataView),
            ('Local Files', self.localFilesWindow),
            ('Remote Files', self.remoteFilesWindow),
            ('Manual', self.showManual())
        ])

        # Setup
        self.setupControlTabs.extend([
            ('Logging', self.loggingControls),
            ('Numpypi', self.numpyPiControls)
        ])
        
    def initializeDashboard(self):

        # Pull all the final dashboard together in an application
        self.dashboard = pn.WidgetBox(
            pn.Column(pn.Row(self.clearLogButton, self.statusWidget), sizing_mode='stretch_width', width_policy='max'),
            pn.Row(self.controlTabs, self.displayTabs)
        )
        
        # Attach application to GridUtils for integration into panel, etc
        # Do this just before launching the application
        self.grd.application(
            app={
                'messages': self.statusWidget,
                'defaultFigureSize': self.defaultPlotFigureSize
            }
        )
