#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:41:08 2017

@author: JoeRippke

This script calculates heat diffusion between the mantle and a vertically
subducting slab of oceanic lithosphere.

I tried to make this realistic by using realistic:
    - slab thickness
    - thermal diffusivity for oceanic lithosphere
    - subduction rate
    - slab and mantle temperatures

I know this is unrealistic in that it:
    - subducts vertically
    - is a homogenous slab (in temperature and composition)
    - does not include mantle convection
    - is a perfectly rectangular slab
    
This is written in Python 3, and tested using iPython on the command line

It gives a warning message about using the tight_layout command. Ignore the
warning. The plot looks much nicer when I use the tight_layout command.

"""

import numpy as np
import matplotlib.pyplot as plt

# domain size [km]
w = h = 100.

# grid resolution, [km]
dx = dy = 1

# Thermal diffusivity of oceanic lithosphere, [m^2/s]
K = 8.7e7

# subduction rate [m/yr]
v = 0.1

# plate thickness [km]
thickness = 50

# slab and mantle temperatures [°C]
Tslab, Tmantle = 0, 1300

# scaling for x and y axes
nx, ny = int(w/dx), int(h/dy)

# variables for use in the discretized diffusion equation
dx2, dy2 = dx*dx, dy*dy
dt = dx2 * dy2 / (2 * K * (dx2 + dy2))

# initial arrays
u0 = Tmantle * np.ones((nx, ny))
u = np.empty((nx, ny))

# initial slab placement at top center
u0[1,25:75] = Tslab

# discretized diffusion equation function
def diffuser(u0, u):
    
    # discretized diffusion equation
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + K * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )
    
    # boundary condition. maintains mantle temp at the boundaries
    for i in range(nx):
        for j in range(ny):
            if i == 0 or i == nx-1:
                u[j,i] = Tmantle
            if j == ny-1:
                u[j,i] = Tmantle
            if (j == 0 and i < 25) or (j == 0 and i > 75):
                u[j,i] = Tmantle

    u0 = u.copy()
    return u0, u

# Number of timesteps
nsteps = 100

# Output 4 figures at these timesteps [Ma]
mfig = [0, 10, 50, 99]
fignum = 0
fig = plt.figure()

# loop to run the diffusion equation function for each time step as the slab
# sinks into the mantle
for m in range(nsteps):
    
    # define the slab sinking
    u0[2:m+1,25:75] = u0[1:m,25:75]
    u0[1,25:75] = Tslab
    
    # run the diffusion eq function
    u0, u = diffuser(u0, u)
    
    # create images at certain time steps
    if m in mfig:
        fignum += 1
        ax = fig.add_subplot(220 + fignum)
        im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=Tslab, vmax=Tmantle)
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('Depth (km)')
        ax.set_title(str(m*10) + ' ka')

# finishing touches to the plot       
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
cbar_ax.set_xlabel('T (°C)', labelpad=20)
fig.colorbar(im, cax=cbar_ax)
plt.tight_layout()

# show the final plot
plt.show()