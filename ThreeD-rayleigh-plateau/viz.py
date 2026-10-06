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
datafile = "ThreeD-marangoni-rise/OLD/N64-64-96-3D/monitor/bubble"
data1 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-marangoni-rise/OLD/N64-64-96-2D/monitor/bubble"
data2 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-marangoni-rise/OLD/N128-128-192-2D/monitor/bubble"
data3 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-marangoni-rise/OLD/N128-128-192-3D/monitor/bubble"
data4 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-marangoni-rise/OLD/N64-64-96-3D-NOSHIFT/monitor/bubble"
data5 = np.loadtxt(datafile, skiprows=2)

rho1 = .2
mu1 = 0.1
k1 = 0.001
sigma0 = 0.1
sigmaT = -0.1
a = 0.5
dT = 2/15

kr = 1
mur = 1
cp = 0.1

U0 = -sigmaT*a*dT/mu1
tNorm = a/U0

Vygb = -2*sigmaT*dT*a/(6*mu1+9*mu1)

muL = mu1 
muG = mu1*kr

kL = k1
kG = k1*kr

gradT = dT*a
R = a 
dSigmadT = sigmaT 

VygbB = -2 * sigmaT * gradT * R/(2*muL + 3*muG) * ((kG+2*kL)/(2*kL+kG))
print("VygbB = ",VygbB)
print("Re = ",rho1*U0*a/mu1)
print("Ca = ",mu1*U0/sigma0)
print("Ma = ",rho1*cp*a*U0/k1)
print("Ca Ur = ",mu1/sigma0)
print("Re Ur = ",rho1*a/mu1)
print("Tnorm = ",tNorm)
print("U0 = ",U0)
print("Vygb =",Vygb)

plt.figure()


x = data1[:,1]/tNorm
y = data1[:,3]/Vygb
plt.plot(x,y,label = "3D,PCST-SHIFT,N=64",linewidth=5,color = color_red)


x = data2[:,1]/tNorm
y = data2[:,3]/Vygb
plt.plot(x,y,label = "2D,PCST-SHIFT,N=64",linewidth=5,color = color_blue)

# x = data3[:,1]/tNorm
# y = data3[:,3]/Vygb
# plt.plot(x,y,label = "2D,PCST,N=128",linewidth=5,color = color_black)

# x = data4[:,1]/tNorm
# y = data4[:,3]/Vygb
# plt.plot(x,y,label = "3D,PCST,N=128",linewidth=5,color = color_orange)

x = data5[:,1]/tNorm
y = data5[:,3]/Vygb
plt.plot(x,y,label = "3D,PCST,N=64",linewidth=5,color = color_teal)
plt.legend()

plt.title("Marangoni Rise Case")
plt.xlabel("t/tnorm")
plt.ylabel("V/Vygb")
plt.show()