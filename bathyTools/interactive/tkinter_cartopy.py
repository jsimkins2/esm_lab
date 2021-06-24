#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 09:14:24 2021

@author: james
"""

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import sys
import tkinter as Tk
import xarray as xr

root = Tk.Tk()
root.wm_title("Cartopy in TK")

gridFile = "/Users/james/Downloads/gridFile.nc"
grid = xr.open_dataset(gridFile)

fig = Figure(figsize=(8,4), dpi=100)
ax = fig.add_subplot(1,1,1, projection=ccrs.LambertConformal())
im = ax.pcolormesh(grid['x'], grid['y'], grid['area'], transform=ccrs.LambertConformal())
ax.set_extent([-150, -100, 21, 70], ccrs.Geodetic())

# a tk.DrawingArea
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.draw()
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

button = Tk.Button(master=root, text='Quit', command=sys.exit)
button.pack(side=Tk.BOTTOM)

Tk.mainloop()
