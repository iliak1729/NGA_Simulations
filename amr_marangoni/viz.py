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
datafile = "amr_marangoni/OLD/amr_3level_noP/monitor/statistics"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "amr_marangoni/OLD/amr_4level_noP/monitor/statistics"
data1b = np.loadtxt(datafile, skiprows=2)
datafile = "amr_marangoni/OLD/amr_3level_P/monitor/statistics"
data1c = np.loadtxt(datafile, skiprows=2)

datafile = "amr_marangoni/LargeDomain_results/amr_Level5NoP/monitor/statistics"
data_NoP = np.loadtxt(datafile, skiprows=2)
datafile = "amr_marangoni/LargeDomain_results/amr_Level5NoP_PPIC/monitor/statistics"
data_NoP_PPIC = np.loadtxt(datafile, skiprows=2)
datafile = "amr_marangoni/LargeDomain_results/amr_Level5P/monitor/statistics"
data2_P = np.loadtxt(datafile, skiprows=2)
datafile = "amr_marangoni/LargeDomain_results/amr_Level5P_PPIC/monitor/statistics"
data2_P_PPIC = np.loadtxt(datafile, skiprows=2)


datafile = "amr_marangoni/Domain_160_240_160_CFL01/amr_Level5NoP_PPIC/monitor/statistics"
data3_NoP_PPIC = np.loadtxt(datafile, skiprows=2)

datafile = "amr_marangoni/Domain_160_240_160/amr_Level5NoP_PPIC/monitor/statistics"
data4_NoP_PPIC = np.loadtxt(datafile, skiprows=2)

datafile = "amr_marangoni/Domain_160_240_160_Level9/amr_Level5NoP/monitor/statistics"
data5_NoP_Level9 = np.loadtxt(datafile, skiprows=2)
datafile = "amr_marangoni/Domain_160_240_160_Level9/amr_Level5P/monitor/statistics"
data5_P_Level9 = np.loadtxt(datafile, skiprows=2)

datafile = "amr_marangoni/Domain_160_240_160_Level9/amr_Level5NoP_PPIC/monitor/statistics"
data5_NoP_PPIC_Level9 = np.loadtxt(datafile, skiprows=2)
datafile = "amr_marangoni/Domain_160_240_160_Level9/amr_Level5P_PPIC/monitor/statistics"
data5_P_PPIC_Level9 = np.loadtxt(datafile, skiprows=2)


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

# Young et al from POpinet
term1 =2/((2+3*mur)*(2+kr))
term2 = sigmaT*a*dT/mu1
VygbY = -term1*term2
print("VygbY = ",VygbY)
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
plt.title("Marangoni Rise Case")

width = 3

x = data_NoP[:,1]/tNorm
y = data_NoP[:,7]/Vygb
plt.plot(x,y,label = "Large Domain,14 Cells/D,NoP,PLIC,OLD",linewidth=width,color = color_black,linestyle='dashed')

x = data2_P[:,1]/tNorm
y = data2_P[:,7]/Vygb
plt.plot(x,y,label = "Large Domain,14 Cells/D,P,PLIC,OLD",linewidth=width,color = color_blue,linestyle='dashed')

x = data5_NoP_Level9[:,1]/tNorm
y = data5_NoP_Level9[:,7]/Vygb
plt.plot(x,y,label = "Larger Domain,26 Cells/D,NoP,PLIC",linewidth=width,color = color_red,linestyle='solid')

x = data5_P_Level9[:,1]/tNorm
y = data5_P_Level9[:,7]/Vygb
plt.plot(x,y,label = "Larger Domain,26 Cells/D,P,PLIC",linewidth=width,color = color_teal,linestyle='solid')

x = data5_NoP_PPIC_Level9[:,1]/tNorm
y = data5_NoP_PPIC_Level9[:,7]/Vygb
plt.plot(x,y,label = "Larger Domain,26 Cells/D,NoP,PPIC",linewidth=width,color = color_orange,linestyle='solid')

x = data5_P_PPIC_Level9[:,1]/tNorm
y = data5_P_PPIC_Level9[:,7]/Vygb
plt.plot(x,y,label = "Larger Domain,26 Cells/D,P,PPIC",linewidth=width,color = color_black,linestyle='solid')

x = data_NoP_PPIC[:,1]/tNorm
y = data_NoP_PPIC[:,7]/Vygb
plt.plot(x,y,label = "Large Domain,14 Cells/D,NoP,PPIC,OLD",linewidth=width,color = color_yellow,linestyle='dashdot')

plt.legend()

plt.title("Marangoni Rise Case")
plt.xlabel("t/tnorm")
plt.ylabel("V/Vygb")
plt.grid(True)
plt.show()