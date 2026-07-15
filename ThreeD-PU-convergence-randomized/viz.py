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

# Get Data
datafile = "ThreeD-PU-convergence-randomized/OLD_SHARE/Results/RandomizedCenters1000/Error_Values_2D_RandomizedCenters"
d1 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence-randomized/OLD_SHARE/Results/RandomizedCenters1000/Error_Values_2D_RandomizedCenters_Jibben"
d2 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence-randomized/OLD_SHARE/Results/RandomizedCenters1000/Error_Values_3D_RandomizedCenters"
d3 = np.loadtxt(datafile, skiprows=2)
datafile = "ThreeD-PU-convergence-randomized/OLD_SHARE/Results/RandomizedCenters1000/Error_Values_3D_RandomizedCenters_Jibben"
d4 = np.loadtxt(datafile, skiprows=2)


# Extract Specific Values
data2D16 = d1[d1[:, 0] == 16]
data2D32 = d1[d1[:, 0] == 32]
data2D64 = d1[d1[:, 0] == 64]
data2D128 = d1[d1[:, 0] == 128]
data2D256 = d1[d1[:, 0] == 256]

data2D16J = d2[d2[:, 0] == 16]
data2D32J = d2[d2[:, 0] == 32]
data2D64J = d2[d2[:, 0] == 64]
data2D128J = d2[d2[:, 0] == 128]
data2D256J = d2[d2[:, 0] == 256]

data3D16 = d3[d3[:, 0] == 16]
data3D32 = d3[d3[:, 0] == 32]
data3D64 = d3[d3[:, 0] == 64]

data3D16J = d4[d4[:, 0] == 16]
data3D32J = d4[d4[:, 0] == 32]
data3D64J = d4[d4[:, 0] == 64]

# Histogram
plt.figure()
plt.hist(data2D16[:, 5]*data2D16[:, 0]*data2D16[:, 0], bins=20, alpha=0.5, label='2D, N=16',weights=np.ones_like(data2D16[:, 5]) / len(data2D16))
plt.hist(data2D32[:, 5]*data2D32[:, 0]*data2D32[:, 0], bins=20, alpha=0.5, label='2D, N=32',weights=np.ones_like(data2D32[:, 5]) / len(data2D32))
plt.hist(data2D64[:, 5]*data2D64[:, 0]*data2D64[:, 0], bins=20, alpha=0.5, label='2D, N=64',weights=np.ones_like(data2D64[:, 5]) / len(data2D64))
plt.hist(data2D128[:, 5]*data2D128[:, 0]*data2D128[:, 0], bins=20, alpha=0.5, label='2D, N=128',weights=np.ones_like(data2D128[:, 5]) / len(data2D128))
plt.hist(data2D256[:, 5]*data2D256[:, 0]*data2D256[:, 0], bins=20, alpha=0.5, label='2D, N=256',weights=np.ones_like(data2D256[:, 5]) / len(data２D２５６))
plt.xlabel('Error')
plt.ylabel('Frequency')
plt.title('Histogram of Radius Errors')
plt.legend()

plt.figure()
plt.hist(data2D16J[:, 5]*data2D16J[:, 0]*data2D16J[:, 0], bins=20, alpha=0.5, label='2D, N=16',weights=np.ones_like(data2D16J[:, 5]) / len(data2D16J))
plt.hist(data2D32J[:, 5]*data2D32J[:, 0]*data2D32J[:, 0], bins=20, alpha=0.5, label='2D, N=32',weights=np.ones_like(data2D32J[:, 5]) / len(data2D32J))
plt.hist(data2D64J[:, 5]*data2D64J[:, 0]*data2D64J[:, 0], bins=20, alpha=0.5, label='2D, N=64',weights=np.ones_like(data2D64J[:, 5]) / len(data2D64J))
plt.hist(data2D128J[:, 5]*data2D128J[:, 0]*data2D128J[:, 0], bins=20, alpha=0.5, label='2D, N=128',weights=np.ones_like(data2D128J[:, 5]) / len(data2D128J))
plt.hist(data2D256J[:, 5]*data2D256J[:, 0]*data2D256J[:, 0], bins=20, alpha=0.5, label='2D, N=256',weights=np.ones_like(data2D256J[:, 5]) / len(data2D256J))
plt.xlabel('Error')
plt.ylabel('Frequency')
plt.title('Histogram of Radius Errors - Jibben')
plt.legend()

plt.figure()
plt.hist(data3D16[:, 5]*data3D16[:, 0]*data3D16[:, 0], bins=20, alpha=0.5, label='3D, N=16',weights=np.ones_like(data3D16[:, 5]) / len(data3D16))
plt.hist(data3D32[:, 5]*data3D32[:, 0]*data3D32[:, 0], bins=20, alpha=0.5, label='3D, N=32',weights=np.ones_like(data3D32[:, 5]) / len(data3D32))
plt.hist(data3D64[:, 5]*data3D64[:, 0]*data3D64[:, 0], bins=20, alpha=0.5, label='3D, N=64',weights=np.ones_like(data3D64[:, 5]) / len(data3D64))

plt.xlabel('Error')
plt.ylabel('Frequency')
plt.title('Histogram of Radius Errors')
plt.legend()

plt.figure()
plt.hist(data3D16J[:, 5]*data3D16J[:, 0]*data3D16J[:, 0], bins=20, alpha=0.5, label='3D, N=16',weights=np.ones_like(data3D16J[:, 5]) / len(data3D16J))
plt.hist(data3D32J[:, 5]*data3D32J[:, 0]*data3D32J[:, 0], bins=20, alpha=0.5, label='3D, N=32',weights=np.ones_like(data3D32J[:, 5]) / len(data3D32J))
plt.hist(data3D64J[:, 5]*data3D64J[:, 0]*data3D64J[:, 0], bins=20, alpha=0.5, label='3D, N=64',weights=np.ones_like(data3D64J[:, 5]) / len(data3D64J))

plt.xlabel('Error')
plt.ylabel('Frequency')
plt.title('Histogram of Radius Errors - Jibben')
plt.legend()

plt.show()

# Average Error vs N
error2D = []
error2DJ = []
error3D = []
error3DJ = []

err2D = []
err2DJ = []
err3D = []
err3DJ = []

datasets2D = [data2D16, data2D32, data2D64, data2D128, data2D256]
datasets2DJ = [data2D16J, data2D32J, data2D64J, data2D128J, data2D256J]
datasets3D = [data3D16, data3D32, data3D64]
datasets3DJ = [data3D16J, data3D32J, data3D64J]

# Mean and standard deviation
for data in datasets2D:
    error2D.append(np.mean(data[:, 5]))
    err2D.append(np.std(data[:, 5], ddof=1))

for data in datasets2DJ:
    error2DJ.append(np.mean(data[:, 5]))
    err2DJ.append(np.std(data[:, 5], ddof=1))

for data in datasets3D:
    error3D.append(np.mean(data[:, 5]))
    err3D.append(np.std(data[:, 5], ddof=1))

for data in datasets3DJ:
    error3DJ.append(np.mean(data[:, 5]))
    err3DJ.append(np.std(data[:, 5], ddof=1))

plt.figure()

N_values_2D = [16.0, 32.0, 64.0, 128.0, 256.0]
N_values_3D = [16.0, 32.0, 64.0]

# plt.errorbar(
#     N_values_2D, error2D, yerr=[abs(error2D[i] - err2D[i]) for i in range(len(error2D))],
#     fmt='o-', capsize=4, label='2D Randomized Centers'
# )

# plt.errorbar(
#     N_values_2D, error2DJ, yerr=[abs(error2DJ[i] - err2DJ[i]) for i in range(len(error2DJ))],
#     fmt='o-', capsize=4, label='2D Randomized Centers - Jibben'
# )

# plt.errorbar(
#     N_values_3D, error3D, yerr=[abs(error3D[i] - err3D[i]) for i in range(len(error3D))],
#     fmt='o-', capsize=4, label='3D Randomized Centers'
# )

# plt.errorbar(
#     N_values_3D, error3DJ, yerr=[abs(error3DJ[i] - err3DJ[i]) for i in range(len(error3DJ))],
#     fmt='o-', capsize=4, label='3D Randomized Centers - Jibben'
# )


plt.errorbar(
    N_values_2D, error2D, yerr=err2D,
    fmt='o-', capsize=4, label='2D Randomized Centers'
)

plt.errorbar(
    N_values_2D, error2DJ, yerr=err2DJ,
    fmt='o-', capsize=4, label='2D Randomized Centers - Jibben'
)

plt.errorbar(
    N_values_3D, error3D, yerr=err3D,
    fmt='o-', capsize=4, label='3D Randomized Centers'
)

plt.errorbar(
    N_values_3D, error3DJ, yerr=err3DJ,
    fmt='o-', capsize=4, label='3D Randomized Centers - Jibben'
)
refIndex = 1
refOrder = 2
x = np.array(N_values_2D)
y = np.array(error2D)
plt.loglog(x,0.8*y[refIndex]*x**(-refOrder)*x[refIndex]**(refOrder),linewidth=2,color = color_black,label = 'Second Order',linestyle='--')

plt.xscale('log')
plt.yscale('log')


plt.xlabel('N')
plt.ylabel('Average Error')
plt.title('Average Error vs N')
plt.legend()
plt.show()



plt.show()

