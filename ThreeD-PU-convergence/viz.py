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

# Radius Convergence Plot =====================================================================================
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_2D_ProjectedMetricConvergence"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_3D_ProjectedMetricConvergence"
data2 = np.loadtxt(datafile, skiprows=2)

plt.figure()
x = data2[:,2]
y = data2[:,3]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_blue,label = '3D')
x = data1[:,2]
y = data1[:,3]
refIndex = 1
refOrder = 2
plt.loglog(x,y,linewidth=5,color = color_red,label = '2D')
plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'Second Order',linestyle='--')


datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_2D_ProjectedMetricConvergence_Jibben"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_3D_ProjectedMetricConvergence_Jibben"
data2 = np.loadtxt(datafile, skiprows=2)
x = data2[:,2]
y = data2[:,3]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_teal,label = '3D Jibben')
x = data1[:,2]
y = data1[:,3]
refIndex = 1
refOrder = 2
plt.loglog(x,y,linewidth=5,color = color_pink,label = '2D Jibben')
# plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'Second Order',linestyle='--')

plt.xlabel('Cells Across Diameter',fontsize=15)
plt.ylabel('L2 Radius error',fontsize=15)
plt.legend()
plt.title("Radius Error Convergence Plot",fontsize=15)

plt.show()

# Tangent Convergence Plot =====================================================================================
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_2D_ProjectedMetricConvergence"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_3D_ProjectedMetricConvergence"
data2 = np.loadtxt(datafile, skiprows=2)

plt.figure()
x = data2[:,2]
y = data2[:,4]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_blue,label = '3D')
# plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='--')
x = data1[:,2]
y = data1[:,4]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_red,label = '2D')
plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='--')


datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_2D_ProjectedMetricConvergence_Jibben"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_3D_ProjectedMetricConvergence_Jibben"
data2 = np.loadtxt(datafile, skiprows=2)
x = data2[:,2]
y = data2[:,4]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_teal,label = '3D Jibben')
# plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='--')
x = data1[:,2]
y = data1[:,4]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_pink,label = '2D Jibben')
plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='--')
refIndex = 1
refOrder = 2
plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'Second Order',linestyle='-.')

plt.xlabel('Cells Across Diameter',fontsize=15)
plt.ylabel('L2 Tangent error',fontsize=15)
plt.legend()
plt.title("Tangent Error Convergence Plot",fontsize=15)

plt.show()


# Tangent Convergence Plot =====================================================================================
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_2D_ProjectedMetricConvergence"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_3D_ProjectedMetricConvergence"
data2 = np.loadtxt(datafile, skiprows=2)

plt.figure()
x = data2[:,2]
y = data2[:,5]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_blue,label = '3D')
# plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='--')
x = data1[:,2]
y = data1[:,5]
refIndex = 1
refOrder = 0
plt.loglog(x,y,linewidth=5,color = color_red,label = '2D')
plt.loglog(x,0.6*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='--')


datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_2D_ProjectedMetricConvergence_Jibben"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetricConvergence/Error_Values_3D_ProjectedMetricConvergence_Jibben"
data2 = np.loadtxt(datafile, skiprows=2)
x = data2[:,2]
y = data2[:,5]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_teal,label = '3D Jibben')
# plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'First Order',linestyle='--')
x = data1[:,2]
y = data1[:,5]
refIndex = 1
refOrder = 1
plt.loglog(x,y,linewidth=5,color = color_pink,label = '2D Jibben')
plt.loglog(x,1.3*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'Zeroth Order',linestyle='--')

plt.xlabel('Cells Across Diameter',fontsize=15)
plt.ylabel('L2 Curvature error',fontsize=15)
plt.legend()
plt.title("Curvature Error Convergence Plot",fontsize=15)

plt.show()



# Spread Variance 2D =====================================================================================
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetric_Sweep/Error_Values_2D_ProjectedMetric_Sweep"
data1 = np.loadtxt(datafile, skiprows=2)

plt.figure()
d1 = data1[data1[:,0]==16]
d2 = data1[data1[:,0]==32]
d3 = data1[data1[:,0]==64]
d4 = data1[data1[:,0]==128]
d5 = data1[data1[:,0]==256]
d6 = data1[data1[:,0]==512]

plt.semilogy(d1[:,1],d1[:,3],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,3],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,3],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,3],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,3],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,3],linewidth=2,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Radius error',fontsize=15)
plt.title("2D Sweep, Radius Error")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,4],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,4],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,4],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,4],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,4],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,4],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Tangent error',fontsize=15)
plt.title("2D Sweep, Tangent Error")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,5],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,5],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,5],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,5],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,5],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,5],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Curvature error',fontsize=15)
plt.title("2D Sweep, Curvature Error")
plt.legend()

plt.show()
# Spread Variance 2D Jibben =====================================================================================
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetric_Sweep/Error_Values_2D_ProjectedMetric_Sweep_Jibben"
data1 = np.loadtxt(datafile, skiprows=2)

plt.figure()
d1 = data1[data1[:,0]==16]
d2 = data1[data1[:,0]==32]
d3 = data1[data1[:,0]==64]
d4 = data1[data1[:,0]==128]
d5 = data1[data1[:,0]==256]
d6 = data1[data1[:,0]==512]

plt.semilogy(d1[:,1],d1[:,3],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,3],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,3],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,3],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,3],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,3],linewidth=2,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Radius error',fontsize=15)
plt.title("2D Sweep, Radius Error,Jibben")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,4],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,4],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,4],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,4],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,4],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,4],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Tangent error',fontsize=15)
plt.title("2D Sweep, Tangent Error,Jibben")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,5],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,5],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,5],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,5],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,5],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,5],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Curvature error',fontsize=15)
plt.title("2D Sweep, Curvature Error,Jibben")
plt.legend()

plt.show()


# Spread Variance 3D =====================================================================================
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetric_Sweep/Error_Values_3D_ProjectedMetric_Sweep"
data1 = np.loadtxt(datafile, skiprows=2)

plt.figure()
d1 = data1[data1[:,0]==16]
d2 = data1[data1[:,0]==32]
d3 = data1[data1[:,0]==64]
d4 = data1[data1[:,0]==128]
d5 = data1[data1[:,0]==256]
d6 = data1[data1[:,0]==512]

plt.semilogy(d1[:,1],d1[:,3],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,3],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,3],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,3],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,3],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,3],linewidth=2,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Radius error',fontsize=15)
plt.title("3D Sweep, Radius Error")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,4],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,4],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,4],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,4],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,4],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,4],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Tangent error',fontsize=15)
plt.title("3D Sweep, Tangent Error")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,5],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,5],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,5],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,5],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,5],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,5],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Curvature error',fontsize=15)
plt.title("3D Sweep, Curvature Error")
plt.legend()

plt.show()
# Spread Variance 3D Jibben =====================================================================================
datafile = "ThreeD-PU-convergence/OLD_SHARE/ProjectedMetric_Sweep/Error_Values_3D_ProjectedMetric_Sweep_Jibben"
data1 = np.loadtxt(datafile, skiprows=2)

plt.figure()
d1 = data1[data1[:,0]==16]
d2 = data1[data1[:,0]==32]
d3 = data1[data1[:,0]==64]
d4 = data1[data1[:,0]==128]
d5 = data1[data1[:,0]==256]
d6 = data1[data1[:,0]==512]

plt.semilogy(d1[:,1],d1[:,3],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,3],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,3],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,3],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,3],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,3],linewidth=2,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Radius error',fontsize=15)
plt.title("3D Sweep, Radius Error, Jibben")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,4],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,4],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,4],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,4],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,4],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,4],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Tangent error',fontsize=15)
plt.title("3D Sweep, Tangent Error, Jibben")
plt.legend()


plt.figure()
plt.semilogy(d1[:,1],d1[:,5],linewidth=5,color = color_black,label = 'N/D=8')
plt.semilogy(d2[:,1],d2[:,5],linewidth=5,color = color_orange,label = 'N/D=16')
plt.semilogy(d3[:,1],d3[:,5],linewidth=5,color = color_teal,label = 'N/D=32')
plt.semilogy(d4[:,1],d4[:,5],linewidth=5,color = color_yellow,label = 'N/D=64')
plt.semilogy(d5[:,1],d5[:,5],linewidth=5,color = color_blue,label = 'N/D=128')
plt.semilogy(d6[:,1],d6[:,5],linewidth=5,color = color_red,label = 'N/D=256')
plt.xlabel('PU Radius, normalized by dx',fontsize=15)
plt.ylabel('L2 Curvature error',fontsize=15)
plt.title("3D Sweep, Curvature Error, Jibben")
plt.legend()

plt.show()
