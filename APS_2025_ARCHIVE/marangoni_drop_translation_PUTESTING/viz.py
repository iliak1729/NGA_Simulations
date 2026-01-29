import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special

datafile = "marangoni_drop_translation/monitor/PUTesting_Marangoni_NewSettings_WithPressure3"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation/monitor/PUTesting_Marangoni_NewSettings_2"
data2 = np.loadtxt(datafile, skiprows=2)

rho1 = 0.2
mu1 = 0.1
sigma0 = 0.1
sigmaT = 0.005
dT = 0.1
a = 0.5


kr = 0
mur = 2

U0 = sigmaT*a*dT/mu1
tNorm = a/U0

Vygb = 2/((2+kr)*(2+3*mur))
print("Re = ",rho1*U0*a/mu1)
print("Ca = ",mu1*U0/sigma0)
print("Ca Ur = ",mu1/sigma0)
print("Re Ur = ",rho1*a/mu1)
print("Tnorm = ",tNorm)
print("U0 = ",U0)
print("Vygb =",Vygb)

x = data1[:,1]
y = data1[:,6]
plt.figure()
plt.plot(x,y,label = 'dt = '+str(round(1000000*data1[1,2])/100) + "e-4")
plt.plot(data2[:,1],data2[:,6],label = 'dt = '+str(round(1000000*data2[1,2])/100) + "e-4")

plt.ylabel("COM Velocity")
plt.xlabel("Time")
plt.legend()

plt.savefig("marangoni_drop_translation/graphs/COM U Velocity New TimeB.png")
plt.show()