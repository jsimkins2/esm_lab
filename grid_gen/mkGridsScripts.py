#!/bin/env python3

# conda: xesmfTools

# Grid toolset demonstration in
#  * Command line script form
#  * ipython

import sys, os

#if (os.environ.get('LIBROOT')):
#    sys.path.append(os.environ.get('LIBROOT'))
sys.path.append('lib')

#from sysinfo import SysInfo
#info = SysInfo()
#info.show(vList=['platform','python','esmf','esmpy','xgcm','xesmf',
#                 'netcdf4','numpy','xarray',
#                 'cartopy','matplotlib',
#                 'jupyter_core','jupyterlab','notebook',
#                 'dask'])

from gridutils import GridUtils

# Initialize a grid object
grd = GridUtils()

# We can turn on extra output from the module
grd.verboseLevel = 1
grd.debugLevel = 1

# Make sure we erase any previous grid, grid parameters and plot parameters.
grd.clearGrid()

# Specify the grid parameters
# gridMode should be 2.0 for supergrid
grd.setGridParameters({
    'projection': {
        'name': 'LambertConformalConic',
        'lon_0': 230.0,
        'lat_0': 40.0
    },
    'dx': 20.0,
    'dxUnits': 'degrees',
    'dy': 30.0,
    'dyUnits': 'degrees',
    'tilt': 30,
    'gridResolution': 1.0,
    'gridMode': 2.0
})

# To set or update dictionary items in 'projection', you can use the dictionary format above with a direct assigment
# or use the subKey parameter as in below.
#grd.setGridParameters({
#    'name': 'LambertConformalConic',
#    'lon_0': 230.0,
#    'lat_0': 40.0
#}, subKey='projection')

# This forms a grid in memory using the specified grid parameters
grd.makeGrid()

# Save the new grid to a netCDF file
grd.saveGrid(filename="configs/test/nikiTest.nc")

# This prints out all the current grid parameters
# Note: for Lambert Conformal Conic grids, two additional projection parameters are computed.
#       First and second parallel for the grid (lat_1 and lat_2)
grd.showGridParameters()

# You can show the data summary from xarray for the grid
grd.grid

# Define plot parameters so we can see what the grid looks like
grd.setPlotParameters(
    {
        'figsize': (8,8),
        'projection': {
            'name': 'NearsidePerspective',
            'lat_0': 40.0,
            'lon_0': 230.0
        },
        'extent': [-160.0 ,-100.0, 20.0, 60.0],
        'iLinewidth': 1.0,
        'jLinewidth': 1.0,
        'showGridCells': True,
        'title': "Nearside Perspective: 20x30 with 30 degree tilt",
        'iColor': 'k',
        'jColor': 'k'
    }
)

# Projection may be specified separately
grd.setPlotParameters(
    {
        'name': 'NearsidePerspective',
        'lat_0': 40.0,
        'lon_0': 230.0        
    }, subKey='projection'
)

# When we call plotGrid() we have two python objects returned
# Figure object - you have control whether to show the
#   figure or save the contents to an output file
# Axes object - you can further fine tune plot parameters,
#   titles, axis, etc prior to the final plotting of the figure.
#   Some items may be configured via the figure object.
(figure, axes) = grd.plotGrid()

# You can save the figure using the savefig() method on the
# figure object.  Many formats are possible.
figure.savefig('configs/test/nikiTest.jpg', dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)

figure.savefig('configs/test/nikiTest.pdf', dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)

