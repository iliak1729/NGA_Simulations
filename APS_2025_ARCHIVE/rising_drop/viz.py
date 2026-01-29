import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special

datafile = "rising_drop/vtk/DropTranslation/Isocontour050.csv"
data = np.loadtxt(datafile, delimiter=',',skiprows=1)

x = np.array(data[:,1])
y = np.array(data[:,2])
xOrdered = np.zeros(x.shape)
yOrdered = np.zeros(x.shape)

# Order Points by closes points
# First Point
xOrdered[0] = x[0]
yOrdered[0] = y[0]
np.delete(x,0)
np.delete(y,0)
count = 1
for j in range(len(xOrdered)):
    minLen = 10000000000000000000
    minIndex = -1
    for i in range(len(x)):
        currX = x[i]
        currY = y[i]

        dx = currX - xOrdered[j-1]
        dy = currY - yOrdered[j-1]
        dL = np.sqrt(dx**2+dy**2)

        if(dL < minLen):
            minLen = dL
            minIndex = i
    xOrdered[j] = x[minIndex]
    yOrdered[j] = y[minIndex]
    x = np.delete(x,minIndex)
    y = np.delete(y,minIndex)


xOrdered = np.append(xOrdered,xOrdered[0])
yOrdered = np.append(yOrdered,yOrdered[0])

plt.figure(figsize=(12,10))


plt.plot(xOrdered,yOrdered,label = "CSS128",linewidth = 3)


datafile = "rising_drop/vtk/DropTranslation/Isocontour050_256x512.csv"
data = np.loadtxt(datafile, delimiter=',',skiprows=1)

x = np.array(data[:,1])
y = np.array(data[:,2])
xOrdered = np.zeros(x.shape)
yOrdered = np.zeros(x.shape)

# Order Points by closes points
# First Point
xOrdered[0] = x[0]
yOrdered[0] = y[0]
np.delete(x,0)
np.delete(y,0)
count = 1
for j in range(len(xOrdered)):
    minLen = 10000000000000000000
    minIndex = -1
    for i in range(len(x)):
        currX = x[i]
        currY = y[i]

        dx = currX - xOrdered[j-1]
        dy = currY - yOrdered[j-1]
        dL = np.sqrt(dx**2+dy**2)

        if(dL < minLen):
            minLen = dL
            minIndex = i
    xOrdered[j] = x[minIndex]
    yOrdered[j] = y[minIndex]
    x = np.delete(x,minIndex)
    y = np.delete(y,minIndex)


xOrdered = np.append(xOrdered,xOrdered[0])
yOrdered = np.append(yOrdered,yOrdered[0])

plt.plot(xOrdered,yOrdered,label = "CSS256",linewidth = 3)
# Popinet 128
datafile = "rising_drop/graphs/Popinet128x256.csv"
data = np.loadtxt(datafile, delimiter=',',skiprows=1)
plt.scatter(data[:,0],data[:,1],label = "Popinet128", s=80, facecolors='none', edgecolors='k')

# Popinet256

datafile = "rising_drop/graphs/Popinet256x512.csv"
data = np.loadtxt(datafile, delimiter=',',skiprows=1)
plt.scatter(data[:,0],data[:,1],label = "Popinet256", s=80, facecolors='none', edgecolors='r',marker = "^")




plt.legend()
plt.axis("equal")

plt.savefig("rising_drop/graphs/Comparison0_50.png")
plt.show()



# Rising Drop Contour