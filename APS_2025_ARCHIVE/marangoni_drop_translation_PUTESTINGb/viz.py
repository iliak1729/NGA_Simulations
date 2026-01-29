import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special


# datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_64x96"
# data0 = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_32x48_Parabolic"
data0P = np.loadtxt(datafile, skiprows=2)

datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_64x96"
data1 = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_64x96_Parabolic"
data1P = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_128x192"
data2 = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_128x192_Parabolic"
data2P = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_256x384"
data3 = np.loadtxt(datafile, skiprows=2)
datafile = "marangoni_drop_translation_PUTESTINGb/monitor/PUTesting_Marangoni_512x768"
data4 = np.loadtxt(datafile, skiprows=2)

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



# x = data0P[:,1]/tNorm
# y = data0P[:,7]/Vygb
# plt.plot(x,y,label = "N=32P",linewidth=3)
# x = data0[:,1]/tNorm
# y = data0[:,7]/Vygb
# plt.plot(x,y,label = "N=32",linewidth=3)


x = data1[:,1]/tNorm
y = data1[:,7]/Vygb
plt.plot(x,y,label = "N=64",linewidth=3)
# x = data1P[:,1]/tNorm
# y = data1P[:,7]/Vygb
# plt.plot(x,y,label = "N=64P",linewidth=3)
x = data2[:,1]/tNorm
y = data2[:,7]/Vygb
plt.plot(x,y,label = "N=128",linewidth=3)
# x = data2P[:,1]/tNorm
# y = data2P[:,7]/Vygb
# plt.plot(x,y,label = "N=128P",linewidth=3)
x = data3[:,1]/tNorm
y = data3[:,7]/Vygb
plt.plot(x,y,label = "N=256",linewidth=3)
# x = data4[:,1]/tNorm
# y = data4[:,7]/Vygb
# plt.plot(x,y,label = "N=512",linewidth=3)

plt.ylabel("$v/v_{ygb}$")
plt.xlabel("Time")
plt.legend()

# plt.savefig("marangoni_drop_translation_PUTESTINGb/graphs/MarangoniUpdated.png")
plt.show()