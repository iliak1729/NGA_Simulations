from turtle import color

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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
datafile = "ThreeD-PU-convergence/monitor/Error Values 2D"
data1 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-PU-convergence/monitor/Error Values 3D"
data2 = np.loadtxt(datafile, skiprows=2)




plt.figure


x = data1[:,1]
y = data1[:,2]
refIndex = 1
refOrder = 2
plt.loglog(x,y,linewidth=5,color = color_red,label = '2D')
plt.loglog(x,0.9*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'Second Order',linestyle='--')

x = data2[:,1]
y = data2[:,2]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_blue,label = '3D')
plt.loglog(x,0.9*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='-')

ax = plt.gca()

ax.xaxis.set_minor_formatter(ticker.FormatStrFormatter('%g'))
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))

plt.xlabel('Cells Across Diameter',fontsize=20)
plt.ylabel('L2 Distance error',fontsize=20)

plt.legend()
plt.show()