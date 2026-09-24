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
datafile = "amr_rayleigh_plateau/Baseline_Runs/amr_CSF/monitor/rayleigh_plateau"
dataCSF = np.loadtxt(datafile, skiprows=2)

datafile = "amr_rayleigh_plateau/Baseline_Runs/amr_CSF_PPIC/monitor/rayleigh_plateau"
dataCSF_PPIC = np.loadtxt(datafile, skiprows=2)

datafile = "amr_rayleigh_plateau/Baseline_Runs/amr_PCST_NoP/monitor/rayleigh_plateau"
dataPCST_NoP = np.loadtxt(datafile, skiprows=2)

datafile = "amr_rayleigh_plateau/Baseline_Runs/amr_PCST_P/monitor/rayleigh_plateau"
dataPCST_P = np.loadtxt(datafile, skiprows=2)

datafile = "amr_rayleigh_plateau/Baseline_Runs/amr_PCST_NoP_PPIC/monitor/rayleigh_plateau"
dataPCST_NoP_PPIC = np.loadtxt(datafile, skiprows=2)

datafile = "amr_rayleigh_plateau/Baseline_Runs/amr_PCST_P_PPIC/monitor/rayleigh_plateau"
dataPCST_P_PPIC = np.loadtxt(datafile, skiprows=2)

datafile = "amr_rayleigh_plateau/monitor/rayleigh_plateau"
dataLIVE = np.loadtxt(datafile, skiprows=2)

R0 = 0.1
rho = 1000
ST_Coeff = 0.072
# Given
ka = 0.697
k = ka/R0
a = R0
plt.figure()

CSF_color = color_black
PCST_NoP_color = color_red
PCST_P_color = color_skyblue

PLIC_style = '-'
PPIC_style = '--'

width = 3
# PLIC Versions
x = dataCSF[:,1]
y = dataCSF[:,3]
plt.plot(x,y,label = "CSF",linewidth=width,color = CSF_color,linestyle=PLIC_style)

x = dataPCST_NoP[:,1]
y = dataPCST_NoP[:,3]
plt.plot(x,y,label = "PCST_NoP",linewidth=width,color = PCST_NoP_color,linestyle=PLIC_style)

x = dataPCST_P[:,1]
y = dataPCST_P[:,3]
plt.plot(x,y,label = "PCST_P",linewidth=width,color = PCST_P_color,linestyle=PLIC_style)

# PPIC Versions
x = dataCSF_PPIC[:,1]
y = dataCSF_PPIC[:,3]
plt.plot(x,y,label = "CSF,PPIC",linewidth=width,color = CSF_color,linestyle=PPIC_style)

x = dataPCST_NoP_PPIC[:,1]
y = dataPCST_NoP_PPIC[:,3]
plt.plot(x,y,label = "PCST_NoP,PPIC",linewidth=width,color = PCST_NoP_color,linestyle=PPIC_style)

x = dataPCST_P_PPIC[:,1]
y = dataPCST_P_PPIC[:,3]
plt.plot(x,y,label = "PCST_P,PPIC",linewidth=width,color = PCST_P_color,linestyle=PPIC_style)



# Exact Value
# Compute Growth Rate
numer = ST_Coeff*k * special.iv(1,ka)*(1-(ka)**2)
denom = rho*a*a*a*special.iv(0,ka)
growth_rate = numer/denom
growth_rate = np.sqrt(growth_rate)

# Exponential Growth 
# plt.plot(x, R0 * np.exp(growth_rate * x), label="Exponential Growth", linewidth=width, color=color_black,linestyle='dashed')


plt.xlim([0,60])


plt.legend()
plt.title("Rayleigh Plateau Instability")
plt.xlabel("t/tc")
plt.ylabel("R(t)")

plt.grid(True)
plt.show()