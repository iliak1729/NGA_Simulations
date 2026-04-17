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
datafile = "APS_2025_ARCHIVE/marangoni_drop_translation_true_runs/monitor/PU_Marangoni_128x192"
data0 = np.loadtxt(datafile, skiprows=2)

datafile = "pressure_term_marangoni/monitor/Marangoni_3D_64_CubicSpline_Shift"
data1 = np.loadtxt(datafile, skiprows=2)

datafile = "pressure_term_marangoni/monitor/Marangoni_3D_64_CubicSpline"
data2 = np.loadtxt(datafile, skiprows=2)
 
datafile = "APS_2025_ARCHIVE/marangoni_drop_translation_true_runs/monitor/Herrmann_Marangoni_256.txt"
data3 = np.loadtxt(datafile, skiprows=1,delimiter=',')


datafile = "pressure_term_marangoni/monitor/Marangoni_Shift_Peskin_64"
data4 = np.loadtxt(datafile, skiprows=2)

datafile = "pressure_term_marangoni/monitor/Marangoni_Seric_Update_128"
data5 = np.loadtxt(datafile, skiprows=2)

datafile = "pressure_term_marangoni/monitor/Marangoni_Shift_Peskin_128"
data6 = np.loadtxt(datafile, skiprows=2)

datafile = "pressure_term_marangoni/monitor/Marangoni_Peskin_128"
data7 = np.loadtxt(datafile, skiprows=2)

rho1 = .2
mu1 = 0.1
sigma0 = 0.1
sigmaT = -0.1
a = 0.5
dT = 2/15

kr = 1
mur = 1

U0 = -sigmaT*a*dT/mu1
tNorm = a/U0

Vygb = 2/((2+kr)*(2+3*mur))
print(Vygb)
Vygb = -2*sigmaT*dT*a/(6*mu1+9*mu1)

print("Re = ",rho1*U0*a/mu1)
print("Ca = ",mu1*U0/sigma0)
print("Ca Ur = ",mu1/sigma0)
print("Re Ur = ",rho1*a/mu1)
print("Tnorm = ",tNorm)
print("U0 = ",U0)
print("Vygb =",Vygb)

plt.figure()

# Cases from Before
# x = data3[:,0]/tNorm
# y = data3[:,1]
# plt.plot(x,y,label = "Herrmann",linewidth=3, color = color_black)

# x = data0[:,1]/tNorm
# y = data0[:,4]/Vygb
# plt.plot(x,y,label = "PCST,N=128",linewidth=5, color = color_skyblue)

# x = data2[:,1]/tNorm
# y = data2[:,4]/Vygb
# plt.plot(x,y,label = "Seric,N=64",linewidth=5,color = color_yellow)

x = data1[:,1]/tNorm
y = data1[:,4]/Vygb
plt.plot(x,y,label = "3D,Shift,N=32",linewidth=5,color = color_red)
x = data2[:,1]/tNorm
y = data2[:,4]/Vygb
plt.plot(x,y,label = "3D,Base",linewidth=5,color = color_blue)
# x = data4[:,1]/tNorm
# y = data4[:,4]/Vygb
# plt.plot(x,y,label = "ShiftPeskin,N=64",linewidth=5, color = color_blue)

# x = data5[:,1]/tNorm
# y = data5[:,4]/Vygb
# plt.plot(x,y,label = "Seric,N=128",linewidth=5,color = color_yellow)

# x = data7[:,1]/tNorm
# y = data7[:,4]/Vygb
# plt.plot(x,y,label = "Peskin,N=128",linewidth=5,color = color_red)

# x = data6[:,1]/tNorm
# y = data6[:,4]/Vygb
# plt.plot(x,y,label = "ShiftPeskin,N=128",linewidth=5, color = color_blue)

# x = data5[:,1]/tNorm
# y = data5[:,4]/Vygb
# plt.plot(x,y,label = "Seric,N=128",linewidth=5,color = color_pink)

# x = data6[:,1]/tNorm
# y = data6[:,4]/Vygb
# plt.plot(x,y,label = "ShiftPeskin64",linewidth=2,color = color_orange)


# x = data3[:,0]/tNorm
# y = data3[:,1]
# plt.plot(x,y,'k--',label = "Herrmann VOF",linewidth=3)

# x = data3B[:,0]/tNorm
# y = data3B[:,1]
# plt.plot(x,y,'k-',label = "Herrmann LS",linewidth=3)


plt.ylabel("$v/v_{ygb}$")
plt.xlabel("Time")
plt.legend()
plt.xlim([0,0.4])
# plt.savefig("marangoni_drop_translation_PUTESTINGb/graphs/MarangoniUpdated.png")
plt.show()