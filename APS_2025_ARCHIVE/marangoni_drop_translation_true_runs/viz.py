import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special


datafile = "marangoni_drop_translation_true_runs/monitor/PU_Marangoni_64x96"
data0 = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_true_runs/monitor/PU_Marangoni_64x96_Jibben"
data0P = np.loadtxt(datafile, skiprows=2)

datafile = "marangoni_drop_translation_true_runs/monitor/PU_Marangoni_128x192"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_true_runs/monitor/PU_Marangoni_128x192_Jibben"
data1P = np.loadtxt(datafile, skiprows=2)


datafile = "marangoni_drop_translation_true_runs/monitor/PU_Marangoni_256x384"
data2 = np.loadtxt(datafile, skiprows=2)
 
datafile = "marangoni_drop_translation_true_runs/monitor/PU_Marangoni_256x384_Jibben"
data2P = np.loadtxt(datafile, skiprows=2)

datafile = "marangoni_drop_translation_true_runs/monitor/Herrmann_Marangoni_256.txt"
data3 = np.loadtxt(datafile, skiprows=1,delimiter=',')
datafile = "marangoni_drop_translation_true_runs/monitor/Herrmann_MarangoniLevelset_256.txt"
data3B = np.loadtxt(datafile, skiprows=1,delimiter=',')


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


x = data0[:,1]/tNorm
y = data0[:,4]/Vygb
plt.plot(x,y,'r-',label = "$D/\Delta=12.8$",linewidth=2)
plt.xlim([0,10])
plt.savefig("marangoni_drop_translation_true_runs/graphs/PU_Marangoni_64x96.png")


x = data0P[:,1]/tNorm
y = data0P[:,4]/Vygb
plt.plot(x,y,'r-',label = "$D/\Delta=12.8$",linewidth=2)
plt.xlim([0,10])
plt.savefig("marangoni_drop_translation_true_runs/graphs/PU_Marangoni_64x96_Jibben.png")
plt.xlim([0,10])

x = data1[:,1]/tNorm
y = data1[:,4]/Vygb
plt.plot(x,y,'r-',label = "$D/\Delta=25.6$",linewidth=2)
plt.xlim([0,10])
plt.savefig("marangoni_drop_translation_true_runs/graphs/PU_Marangoni_128x192.png")
plt.xlim([0,10])

x = data1P[:,1]/tNorm
y = data1P[:,4]/Vygb
plt.plot(x,y,'r-',label = "$D/\Delta=25.6$",linewidth=2)
plt.xlim([0,10])
plt.savefig("marangoni_drop_translation_true_runs/graphs/PU_Marangoni_128x192_Jibben.png")
plt.xlim([0,10])

x = data2[:,1]/tNorm
y = data2[:,4]/Vygb
plt.plot(x,y,'r-',label = "N=256",linewidth=2)
plt.xlim([0,10])
plt.savefig("marangoni_drop_translation_true_runs/graphs/PU_Marangoni_256x384.png")
plt.xlim([0,10])

x = data2P[:,1]/tNorm
y = data2P[:,4]/Vygb
plt.plot(x,y,'r-',label = "N=256P",linewidth=2)
plt.xlim([0,10])
plt.savefig("marangoni_drop_translation_true_runs/graphs/PU_Marangoni_256x384_Jibben.png")
plt.xlim([0,10])

# x = data3[:,0]/tNorm
# y = data3[:,1]
# plt.plot(x,y,'k--',label = "Herrmann VOF",linewidth=3)

# x = data3B[:,0]/tNorm
# y = data3B[:,1]
# plt.plot(x,y,'k-',label = "Herrmann LS",linewidth=3)


plt.ylabel("$v/v_{ygb}$")
plt.xlabel("Time")
plt.legend()
plt.xlim([0,10])
# plt.savefig("marangoni_drop_translation_PUTESTINGb/graphs/MarangoniUpdated.png")
plt.show()