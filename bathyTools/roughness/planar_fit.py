#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 09:14:24 2021

@author: james
"""

import numpy as np
import scipy.optimize

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt



#determinant of matrix a
def det(a):
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[0][0]*a[1][2]*a[2][1]

#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#dot product of vectors a and b
def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

#cross product of vectors a and b
def cross(a, b):
    x = a[1] * b[2] - a[2] * b[1]
    y = a[2] * b[0] - a[0] * b[2]
    z = a[0] * b[1] - a[1] * b[0]
    return (x, y, z)

#area of polygon poly


#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)


fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.scatter(data[:, 0], data[:, 1], data[:, 2])
xx = np.array([[0,10],[0,10]])
yy  = np.array([[0,0],[3,3]])
zz = np.array([[0,1],[0,1]])
# compute needed points for plane plotting

# plot plane
ax.plot_surface(xx, yy, zz, vmin=0, vmax=10, alpha=0.5)
ax.set_zlim3d(0,4)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()







import numpy as np
import matplotlib.pyplot as plt
from   mpl_toolkits.mplot3d import Axes3D
from   math import pow, sqrt

pts = np.add.accumulate(np.random.random((10,3)))

x, y, z = pts.T

# plane parallel to the y-axis
A_xz = np.vstack((x, np.ones(len(x)))).T
m_xz, c_xz = np.linalg.lstsq(A_xz, z, rcond=None)[0]

# plane parallel to the x-axis
A_yz = np.vstack((y, np.ones(len(y)))).T
m_yz, c_yz = np.linalg.lstsq(A_yz, z, rcond=None)[0]

# the intersection of those two planes and
# the function for the line would be:
# z = m_yz * y + c_yz
# z = m_xz * x + c_xz
# or:
def lin(z):
    x = (z - c_xz)/m_xz
    y = (z - c_yz)/m_yz
    return x,y


# get 2 points on the intersection line 
za = z[0]
zb = z[len(z) - 1]
xa, ya = lin(za)
xb, yb = lin(zb)

# get distance between points
len = sqrt(pow(xb - xa, 2) + pow(yb - ya, 2) + pow(zb - za, 2))

# get slopes (projections onto x, y and z planes)
sx = (xb - xa) / len  # x slope
sy = (yb - ya) / len  # y slope
sz = (zb - za) / len  # z slope

# integrity check - the sum of squares of slopes should equal 1.0
# print (pow(sx, 2) + pow(sy, 2) + pow(sz, 2))

fig = plt.figure()
ax = Axes3D(fig)
ax.set_xlabel("x, slope: %.4f" %sx, color='blue')
ax.set_ylabel("y, slope: %.4f" %sy, color='blue')
ax.set_zlabel("z, slope: %.4f" %sz, color='blue')
ax.scatter(x, y, z)
ax.plot([xa], [ya], [za], markerfacecolor='k', markeredgecolor='k', marker = 'o')
ax.plot([xb], [yb], [zb], markerfacecolor='k', markeredgecolor='k', marker = 'o')
ax.plot([xa, xb], [ya, yb], [za, zb], color = 'r')

plt.show()



