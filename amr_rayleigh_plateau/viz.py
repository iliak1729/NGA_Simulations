from turtle import color

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special

# Colors
color_black = '#000000'
color_orange = '#E69F00'
color_skyblue = '#56B4E9'
color_teal = '#009E73'
color_yellow = '#F0E442'
color_blue = '#0072B2'
color_red = '#D55E00'
color_pink = '#CC79A7'
# Data
datafile = "amr_rayleigh_taylor/monitor/rayleigh_plateau"
dataLIVE = np.loadtxt(datafile, skiprows=2)


plt.figure()

width = 3
x = dataLIVE[:,1]
y = dataLIVE[:,3]
plt.plot(x,y,label = "Rmax",linewidth=width,color = color_red)

# x = dataLIVE[:,2]
# y = dataLIVE[:,4]
# plt.plot(x,y,label = "Rmin",linewidth=width,color = color_blue)

plt.legend()

plt.title("Rayleigh Plateau Instability")
plt.xlabel("t/tc")
plt.ylabel("R(t)")

plt.grid(True)
plt.show()