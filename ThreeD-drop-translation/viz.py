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
datafile = "ThreeD-drop-translation/SANITY_RUNS/NONAMR_PCST_N64/monitor/simulation"
data1 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-drop-translation/SANITY_RUNS/NONAMR_PUST_N64_NoP/monitor/simulation"
data2 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-drop-translation/SANITY_RUNS/NONAMR_PUST_N64_P/monitor/simulation"
data3 = np.loadtxt(datafile, skiprows=2)

x = data1[:,1]
y = data1[:,3]
Tnorm120 = x[-1]
Tnorm1200 = np.sqrt(10)*Tnorm120
Tnorm12000 = np.sqrt(10)*Tnorm1200

plt.semilogy(x/Tnorm120,y,label = "Paraboloid",linewidth=5,color = color_teal)
plt.legend()


x = data2[:,1]
y = data2[:,3]
plt.semilogy(x/Tnorm120,y,label = "Pu_NoP",linewidth=5,color = color_yellow)
plt.legend()

x = data3[:,1]
y = data3[:,3]
plt.semilogy(x/Tnorm120,y,label = "Pu_P",linewidth=5,color = color_red)
plt.legend()


plt.title("Laplace Eq")
plt.xlabel("t/Tv")
plt.ylabel("Vrms")
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.8, alpha=0.7)
plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.4)

plt.show()