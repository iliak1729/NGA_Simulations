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

# Get Data
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/WithSpeed/Error_Values_2D_EllipsoidStress_Reordered"
d1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/WithSpeed/Error_Values_2D_EllipsoidStress_Reordered_Jibben"
d2 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/WithSpeed/Error_Values_3D_EllipsoidStress_Reordered"
d3 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/WithSpeed/Error_Values_3D_EllipsoidStress_Reordered_Jibben"
d4 = np.loadtxt(datafile, skiprows=2)

# xIndex = 1
# yIndex = 5
# widthSize = 5
# markerForm = '.'
# markerSize = 15
# plt.figure()
# x = d2[:,xIndex]
# y = d2[:,yIndex]
# refIndex = 1
# refOrder = 1
# plt.loglog(x,y,linewidth=widthSize,color = color_blue,label = '2D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
# x = d1[:,xIndex]
# y = d1[:,yIndex]
# refIndex = 1
# refOrder = 2
# plt.loglog(x,y,linewidth=widthSize,color = color_red,label = '2D',marker = markerForm,markersize = markerSize)
# plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Second Order',linestyle='--')


# x = d4[:,xIndex]
# y = d4[:,yIndex]
# refIndex = 1
# refOrder = 3
# plt.loglog(x,y,linewidth=widthSize,color = color_teal,label = '3D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
# plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Third Order',linestyle='-.')
# x = d3[:,xIndex]
# y = d3[:,yIndex]
# refIndex = 1
# refOrder = 2
# plt.loglog(x,y,linewidth=widthSize,color = color_pink,label = '3D',marker = markerForm,markersize = markerSize)

# plt.xlabel('Cells Across Minor Axis',fontsize=15)
# plt.ylabel('L2 Radius error',fontsize=15)
# plt.legend()
# plt.title("Radius Error Convergence Plot",fontsize=15)

# plt.show()

# # ============================= Normal PLOT ========================
# xIndex = 1
# yIndex = 6
# widthSize = 5
# markerForm = '.'
# markerSize = 15
# plt.figure()
# x = d2[:,xIndex]
# y = d2[:,yIndex]
# refIndex = 3
# refOrder = 5
# plt.loglog(x,y,linewidth=widthSize,color = color_blue,label = '2D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
# # plt.loglog(x,1.5*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Fifth Order',linestyle='-.')
# refIndex = 0
# refOrder = 4
# plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Fourth Order',linestyle=':')

# x = d1[:,xIndex]
# y = d1[:,yIndex]
# refIndex = 1
# refOrder = 2
# plt.loglog(x,y,linewidth=widthSize,color = color_red,label = '2D',marker = markerForm,markersize = markerSize)
# plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Second Order',linestyle='--')

# x = d4[:,xIndex]
# y = d4[:,yIndex]
# refIndex = 1
# refOrder = 5
# plt.loglog(x,y,linewidth=widthSize,color = color_teal,label = '3D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')

# x = d3[:,xIndex]
# y = d3[:,yIndex]
# plt.loglog(x,y,linewidth=widthSize,color = color_pink,label = '3D',marker = markerForm,markersize = markerSize)

# plt.xlabel('Cells Across Minor Axis',fontsize=15)
# plt.ylabel('L2 Normal error',fontsize=15)
# plt.legend()
# plt.title("Normal Error Convergence Plot",fontsize=15)

# plt.show()

# # ============================= Curvature PLOT ========================
# xIndex = 1
# yIndex = 7
# widthSize = 5
# markerForm = '.'
# markerSize = 15
# plt.figure()
# x = d2[:,xIndex]
# y = d2[:,yIndex]

# plt.loglog(x,y,linewidth=widthSize,color = color_blue,label = '2D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
# refIndex = 0
# refOrder = 1
# plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'First Order',linestyle='-.')

# x = d1[:,xIndex]
# y = d1[:,yIndex]
# refIndex = 1
# refOrder = 0
# plt.loglog(x,y,linewidth=widthSize,color = color_red,label = '2D',marker = markerForm,markersize = markerSize)
# plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Zeroth Order',linestyle='--')


# x = d4[:,xIndex]
# y = d4[:,yIndex]

# plt.loglog(x,y,linewidth=widthSize,color = color_teal,label = '3D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')

# x = d3[:,xIndex]
# y = d3[:,yIndex]

# plt.loglog(x,y,linewidth=widthSize,color = color_pink,label = '3D',marker = markerForm,markersize = markerSize)

# plt.xlabel('Cells Across Minor Axis',fontsize=15)
# plt.ylabel('L2 Curvature error',fontsize=15)
# plt.legend()
# plt.title("Curvature Error Convergence Plot",fontsize=15)

# plt.show()

# # ============================= Total Stress PLOT ========================
# xIndex = 1
# yIndex = 17
# widthSize = 5
# markerForm = '.'
# markerSize = 15
# plt.figure()
# x = d2[:,xIndex]
# y = d2[:,yIndex]

# plt.loglog(x,y,linewidth=widthSize,color = color_blue,label = '2D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
# x = d1[:,xIndex]
# y = d1[:,yIndex]
# refIndex = 0
# refOrder = 0.5
# plt.loglog(x,y,linewidth=widthSize,color = color_red,label = '2D',marker = markerForm,markersize = markerSize)
# plt.loglog(x,0.9*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Order 0.5',linestyle='--')


# x = d4[:,xIndex]
# y = d4[:,yIndex]

# plt.loglog(x,y,linewidth=widthSize,color = color_teal,label = '3D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
# refIndex = 0
# refOrder = 2
# plt.loglog(x,0.9*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Second Order',linestyle='--')

# x = d3[:,xIndex]
# y = d3[:,yIndex]

# plt.loglog(x,y,linewidth=widthSize,color = color_pink,label = '3D',marker = markerForm,markersize = markerSize)

# plt.xlabel('Cells Across Minor Axis',fontsize=15)
# plt.ylabel('L2 Total Stress error',fontsize=15)
# plt.legend()
# plt.title("Total Stress Error Convergence Plot",fontsize=15)

# plt.show()


# ============================= Other Stress Plots ========================
xIndex = 1
yIndex = 20
widthSize = 5
markerForm = '.'
markerSize = 15
plt.figure()
x = d2[:,xIndex]
y = d2[:,yIndex]

plt.loglog(x,y,linewidth=widthSize,color = color_blue,label = '2D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
x = d1[:,xIndex]
y = d1[:,yIndex]
refIndex = 0
refOrder = 0.5
plt.loglog(x,y,linewidth=widthSize,color = color_red,label = '2D',marker = markerForm,markersize = markerSize)
plt.loglog(x,0.9*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'Order 0.5',linestyle='--')


x = d4[:,xIndex]
y = d4[:,yIndex]

plt.loglog(x,y,linewidth=widthSize,color = color_teal,label = '3D Jibben',marker = markerForm,markersize = markerSize,linestyle='--')
refIndex = 0
refOrder = 1
plt.loglog(x,0.9*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=widthSize/2,color = color_black,label = 'First Order',linestyle='--')

x = d3[:,xIndex]
y = d3[:,yIndex]

plt.loglog(x,y,linewidth=widthSize,color = color_pink,label = '3D',marker = markerForm,markersize = markerSize)

plt.xlabel('Cells Across Minor Axis',fontsize=15)
plt.ylabel('L2 Z Force error',fontsize=15)
plt.legend()
plt.title("Z Force Error Convergence Plot",fontsize=15)


plt.show()






