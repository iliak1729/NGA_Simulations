

import numpy as np
import matplotlib.pyplot as plt
import math

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
datafile = "ThreeD-marangoni-rise/monitor/bubble"
data1 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-marangoni-rise/NewProperties/P_NoMatching_Jibben_64/monitor/bubble"
data2 = np.loadtxt(datafile, skiprows=2)


mug = 0.001
mul = 0.001

rhol = 1 
rhog = 100

sigma0 = 0.072
sigmaT =-0.1
dT = 22500

a = 0.00001/2.0
L = mug/mul
U0 = -sigmaT*a*dT/mul
taup = rhog*(2*a)*(2*a)*(3+3*L)/(18*mul*(2+2*L))
tNorm = a/U0
# tNorm = taup
# Vygb = -2*sigmaT*dT*a/(6*mu1+9*mu1)
Vygb = -2*a*dT*sigmaT/(6*mul + 9*mug)

gradT = dT*a
R = a 
dSigmadT = sigmaT 

# VygbB = -2 * sigmaT * gradT * R/(2*muL + 3*muG) * ((kG+2*kL)/(2*kL+kG))
# print("VygbB = ",VygbB)
# print("Re = ",rho1*U0*a/mu1)
# print("Ca = ",mu1*U0/sigma0)
# print("Ma = ",rho1*cp*a*U0/k1)
# print("Ca Ur = ",mu1/sigma0)
# print("Re Ur = ",rho1*a/mu1)
print("Tnorm = ",tNorm)
print("U0 = ",U0)
print("Vygb =",Vygb)
plt.figure()

# x = data5[:,1]
# y = data5[:,3]/Vygb
# plt.plot(x,y,label = "3D,OLD PU,N=64",linewidth=5,color = color_teal)


x = data2[:,1]/tNorm
y = data2[:,3]/Vygb
plt.plot(x,y,label = "3D,Face Paraboloid,12.8 Cells/Dia.",linewidth=5,color = color_teal)

x = data1[:,1]/tNorm
y = data1[:,3]/Vygb
plt.plot(x,y,label = "3D,Face Paraboloid,19.2 Cells/Dia.",linewidth=5,color = color_red)



# x = data2[:,1]/tNorm
# y = data2[:,3]/Vygb
# plt.plot(x,y,label = "2D,PCST-SHIFT,N=64",linewidth=5,color = color_blue)

# x = data3[:,1]/tNorm
# y = data3[:,3]/Vygb
# plt.plot(x,y,label = "2D,PCST,N=128",linewidth=5,color = color_black)

# x = data4[:,1]/tNorm
# y = data4[:,3]/Vygb
# plt.plot(x,y,label = "3D,PCST,N=128",linewidth=5,color = color_orange)


plt.legend()

plt.title("Marangoni Rise Case")
plt.xlabel("t/tnorm")
plt.ylabel("V/Vygb")
plt.show()