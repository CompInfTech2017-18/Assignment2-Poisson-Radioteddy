#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  3 01:20:16 2018

@author: radioteddy
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

#starting distribution of potential for right colorbar mapping
U = np.loadtxt('/home/radioteddy/scripts/poisson/distribution/potential_'+str(100)+'.dat')

x = np.linspace(0, 100, 100)
y = np.linspace(0, 100, 100)
xgrid, ygrid = np.meshgrid(x, y)

fig = plt.figure()
ax = fig.gca(projection='3d')

face = ax.plot_surface(xgrid, ygrid, U, cmap=cm.plasma, linewidth=0, antialiased=False) 
       
ax.set_zlim(0, 100.0)
ax.zaxis.set_major_locator(LinearLocator(6))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))
cb = fig.colorbar(face, shrink=0.5, aspect=10, ticks=[20, 40, 60, 80])


def data_gen(frame, U, plot):
    U = np.loadtxt('/home/radioteddy/scripts/poisson/distribution/potential_'+str(frame+1)+'.dat')
    ax.clear()
    plot = ax.plot_surface(xgrid, ygrid, U, cmap=cm.plasma, linewidth=0, antialiased=False)
    return plot,

face_ani = animation.FuncAnimation(fig, data_gen, frames=range(100), fargs=(U, face), interval=1, blit=False, repeat=False)

plt.show()
