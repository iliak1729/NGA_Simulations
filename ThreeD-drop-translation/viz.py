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

# Physics
sigma = 1
rho = 1
D = 0.4
La = 120
mu = np.sqrt(sigma*rho*D/La)
Usigma = 1.0
Tmu = rho*D**2/mu

print("Usigma = ", Usigma)
print("Tmu = ", Tmu)
print("Tmax = ", Tmu*0.013)
    
# Data
datafile = "ThreeD-drop-translation/LONG_RUNS/TranslatingDrop_PCST_N64_NoP_Matching_LA120/monitor/simulation"
data1 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-drop-translation/LONG_RUNS/TranslatingDrop_PCST_N64_P_Matching_LA120/monitor/simulation"
data2 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-drop-translation/LONG_RUNS/TranslatingDrop_PCST_N64_NoP_NoMatching_LA120/monitor/simulation"
data3 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-drop-translation/LONG_RUNS/TranslatingDrop_PCST_N64_P_NoMatching_LA120/monitor/simulation"
data4 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-drop-translation/LONG_RUNS/TranslatingDrop_CSF_N64_LA120/monitor/simulation"
data5 = np.loadtxt(datafile, skiprows=2)

datafile = "ThreeD-drop-translation/monitor/simulation"
data6 = np.loadtxt(datafile, skiprows=2)


x = data1[:,1]
y = data1[:,3]
Tnorm120 = Tmu
Tnorm1200 = np.sqrt(10)*Tnorm120
Tnorm12000 = np.sqrt(10)*Tnorm1200

plt.semilogy(x/Tnorm120,y/Usigma,label = "Matching, No P",linewidth=5,color = color_teal)


x = data2[:,1]
y = data2[:,3]
plt.semilogy(x/Tnorm120,y/Usigma,label = "Matching, P",linewidth=5,color = color_yellow)


x = data3[:,1]
y = data3[:,3]
plt.semilogy(x/Tnorm120,y/Usigma,label = "No Matching, No P",linewidth=5,color = color_red)


x = data4[:,1]
y = data4[:,3]
plt.semilogy(x/Tnorm120,y/Usigma,label = "No Matching, P",linewidth=5,color = color_blue)


x = data5[:,1]
y = data5[:,3]
plt.semilogy(x/Tnorm120,y/Usigma,label = "CSF",linewidth=5,color = color_pink)

x = data6[:,1]
y = data6[:,3]
plt.semilogy(x/Tnorm120,y/Usigma,label = "PUST,WithP",linewidth=5,color = color_black)

plt.legend(loc='best',bbox_to_anchor=(0.8, 0.5))


plt.title("Laplace Eq")
plt.xlabel("t/Tv")
plt.ylabel("Vrms/Ucap")
plt.minorticks_on()
plt.grid(which='major', linestyle='-', linewidth=0.8, alpha=0.7)
plt.grid(which='minor', linestyle=':', linewidth=0.5, alpha=0.4)

plt.show()