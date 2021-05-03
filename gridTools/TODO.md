# Planned work

A milestone for version 1.0 has yet to be established.

# TASKS

 - [ ] grid creation/editor
   - [ ] grid metrics
     - [X] Spherical solution is complete via Niki's ROMS to MOM6 converter
     - [ ] Mercator (angle_dx might be 0 as it is lined up along latitude lines; except for tilt?)
     - [ ] Polar (might be the same as spherical?)
   - [ ] make Lambert Conformal Conic Grids; needs testing
     - [ ] LCC cannot take custom lat_1 and lat_2; it generates lat_1 and lat_2 based on grid inputs
     - [X] Update new lat_1 and lat_2 for application after makeGrid() is run
     - [ ] changing plot parameters lat_1 and lat_2 do not seem to impact the view
   - [ ] grid generation in other projections
   - [ ] on saveGrid() convert lon [+0,+360] to [-180,+180]
   - [ ] Unify ellipse radius (R) constants throughout code
     - [ ] Allow user control?
 - [ ] grid mask editor (land, etc)
 - [ ] integration of bathymetric sources and apply to grids
       Niki: https://github.com/nikizadehgfdl/ocean_model_topog_generator
 - [X] add nbserverproxy/xgcm to conda software stacks; copied to binder environment.yml
 - [ ] Add option to use Alistair's numpypi package as a configurable option in toolsets
 - [ ] turn numpypi into a loadable package via pip
 - [X] add datashader and numpypi from github sources; see postBuild script
   - [ ] document this need and steps for iterative methods
 - [X] xarray \_FillValue needs to be turned off somehow
 - [X] place display(dashboard) as a separate notebook cell
 - [ ] on load of a grid
   - [ ] calculate R
   - [ ] calculate tilt (may not be possible)
   - [ ] update any tool metadata that is appropriate for that grid
 - [X] Create an application method within the GridUtils() class; GridTools().app()

# TODO

 - [X] Further consolidate matplotlib plotting code
   - [X] Refactor plotting code.  It is mostly the same except for setting the projection.
 - [ ] Add "Refresh Plot" buttons to other Plot tabs or figure out how to squeeze a single plot button into the layout
 - [ ] Do we have to declare everything in __init__ first or can be push all that to respective reset/clear functions?
 - [ ] refactor messaging/logging out of GridUtils into its own package so we can import printMsg/debugMsg as standalone calls
 - [ ] refactor refineS and refineR options as Niki had them defined; allow working in meters too
 - [ ] makeGrid assumes degrees at this point.
 - [X] Pass back an error graphic instead of None for figures that do not render
 - [ ] Add a formal logging/message mechanism.
   - [X] Allow display of important messages and warnings in panel application: widget=TextAreaInput
   - [X] Create options in application and other tools for user configuration of logging and output.
   - [X] Create a message buffer/system for information.
   - [ ] Create a separate app to watch a log file? https://discourse.holoviz.org/t/scrollable-log-text-viewer/317
 - [ ] For now, the gridParameters are always in reference to a center point in a grid
   in the future, one may fix a side or point of the grid and grow out from that point
   instead of the center.
 - [ ] Add testing harnesses.
   - [ ] pytest: This will allow testing of core code via command line and iterative methods.
   - [ ] selenium: Testing interactive methods may be harder.

# WISH

 - [ ] tripolar grids: use FRE-NCtools via cython?
 - [ ] Bring in code that converts ROMS grids to MOM6 grids
   - [ ] Allow conversion of MOM6 grids to ROMS grids
 - [ ] grid reading and plot parameter defaults should be dynamic with grid type declaration and potentially
       split out into separate library modules? lib/gridTools/grids/{MOM6,ROMS,WRF}
 - [ ] Place additional projection metadata into MOM6 grid files
   - [X] Added proj string to netCDF file
   - [ ] Tri polar grid description
