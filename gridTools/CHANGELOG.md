# Change Log

# 2021-04-10

 - Shore up messaging and debugging code in GridUtils().  A lot of missing level= in 2nd arguments
   calls to printMsg and debugMsg.
 - Add a TODO to refactor messaging and debugging into its own package/module.
 - Add an example on how to work with logging levels and debug levels.
 - Add showPlotParameters function.
 - Add more explanation to example1.

# 2021-04-09

 - Fix warning in GridUtils.plotGrid()
 - Change warnings to logging.WARNING messages
 - Moved all print statements to printMsg calls
 - Fix printing to STDOUT when msgBox is not defined
 - Add a function that allows tweaking of noisy python modules that send information to the log.
 - Add Manual documentation
 - Add application Setup tab for other obscure toolset options; add setter and getter methods in GridUtils()
   - Use numpypi: True/False
   - Enable logging: True/False
   - Specify logfile: "filename"
   - Specify logging level
   - Specify verbose level
   - Specify debug level
   - Log erase button
 - Update some important references
 - We can add param.watch to any control to trigger events when certain things happen.
 - We also learned that a lot of python modules leverage the logging module and that some of
   those modules are very verbose.  We setup a function to reduce some of that noise.
 - Once a logger is created, it cannot be deleted.  It can be enabled and disabled.
 - Add a small program example to show all available loggers after a GridUtils object is created.
 - TODO: Creating more small program examples to demonstrate logging and debugging techniques.
 - Use the setter functions in mkGridScripts and examples instead of setting the object variables directly.
   - TODO: consider moving important variables to private/hidden variables.

# 2021-04-08

 - Experimentation with panel.pane.HTML did not work pan out.  No great control over width and height.
   Text updates did not automatically resize the window.  The TextAreaInput automatically adds a
   scrollbar to the box when enough lines are added to the window.
 - Fixed up documentation of Grid Representation in the app manual.
 - Panel markdown honors the usage of options after a link. 
   Ex: [MOM6 User Manual](https://mom6.readthedocs.io/){target="\_blank"}
 - Testing of new printMsg facility is working.
 - Added a clear information button to clear the inforamtion window.

# 2021-04-07

 - Move Spherical.py to spherical.py to match coding standards
 - Use R from GridUtils class in spherical

# 2021-04-03

 - BUG: add ccrs.SouthPolarStereo() to projCarto
 - Remove plotExtentX0,X1 checks for lon>180
 - Moving the boilerplate into its own module works
 - BUG(migration): makeGrid() contained a small bug, parallels were not updating; missing self
 - reworked the way we start the app via show() and display()
 - BUG(migration): Plot button stopped working; missing self 
 - BUG(migration): Saving local files are fixed
 - BUG(migration): Call the showManual method with ()

# 2021-04-02

 - combined xgcmTools with xecmfTools configuration
 - added nbserverproxy to xecmfTools
 - moved documentation around
 - merge code updates from James
   - proj; applied application changes into gridutils so it is available to all
   - applied widget clean up for local file selection
 - add todos: LCC limitations; consider numpypi code
 - initial move to hide application boilerplate; help diff/debugging
 - discovered how to suppress xarray \_FillValue attributes
 - created mkGridScripts.py to demonstrate command line/ipython use
 - more documentation needed
 - document how numpypi and datashader should be installed
 - update extent display to show -180 to +180
 - most spinner widgets updated to numeric values instead of integer
 - establish a global for the defaultGridFilename
 - still need to understand how libraries, modules and packages are handled in python
 - migration of mkMap to mkGrid filename; we started with maps but have moved onto grids

