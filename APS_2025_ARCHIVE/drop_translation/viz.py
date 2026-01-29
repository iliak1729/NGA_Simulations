import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special


datafile1 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La120_CSF"
data1 = np.loadtxt(datafile1, skiprows=2)

datafile2 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La1200_CSF"
data2 = np.loadtxt(datafile2, skiprows=2)

datafile3 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La12000_CSF"
data3 = np.loadtxt(datafile3, skiprows=2)

datafile4 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La120_CSS"
data4 = np.loadtxt(datafile4, skiprows=2)

datafile5 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La1200_CSS"
data5 = np.loadtxt(datafile5, skiprows=2)

datafile6 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La12000_CSS"
data6 = np.loadtxt(datafile6, skiprows=2)

datafile7 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La120_CSF_Jibben"
data7 = np.loadtxt(datafile7, skiprows=2)

datafile8 = "drop_translation/monitor/A_PUTesting_TranslatingDrop_La120_CSS_Jibben"
data8 = np.loadtxt(datafile8, skiprows=2)


plt.figure()
plt.title(datafile1)
x = data1[0::50,1]/data1[1,3]
y = data1[0::50,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSF120.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSF120")
# plt.savefig("drop_translation/graphs/La120CSF.png")

# plt.title(datafile2)
x = data2[0::110,1]/data2[1,3]
y = data2[0::110,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSF1200.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSF1200")
# plt.savefig("drop_translation/graphs/La1200CSF.png")

# plt.title(datafile3)
x = data3[0::800,1]/data3[1,3]
y = data3[0::800,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSF12000.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSF12000")
# plt.savefig("drop_translation/graphs/La12000CSF.png")

# plt.title(datafile4)
x = data4[0::50,1]/data4[1,3]
y = data4[0::50,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSS120.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSS120")
# plt.savefig("drop_translation/graphs/La120CSS.png")

# plt.title(datafile5)
x = data5[0::110,1]/data5[1,3]
y = data5[0::110,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSS1200.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSS1200")
# plt.savefig("drop_translation/graphs/La1200CSS.png")

# plt.figure()
# plt.title(datafile6)
x = data6[0::800,1]/data6[1,3]
y = data6[0::800,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSS12000.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSS12000")
# plt.savefig("drop_translation/graphs/La12000CSS.png")

# plt.figure()
# plt.title(datafile7)
x = data7[0::50,1]/data7[1,3]
y = data7[0::50,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSF120J.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSF120J")
# plt.savefig("drop_translation/graphs/La120CSF_Jibben.png")

# plt.figure()
# plt.title(datafile8)
x = data8[0::50,1]/data8[1,3]
y = data8[0::50,5]
output = np.array([x,y])
output = np.transpose(output)
np.savetxt("CSS120J.csv",output,delimiter=",",fmt="%g")
plt.semilogy(x,y,label = "CSS120J")
# plt.savefig("drop_translation/graphs/La120CSS_Jibben.png")


plt.ylabel("RMS Velocity")
plt.xlabel("Time")
plt.xlim([0,1])
plt.legend()

# plt.savefig("marangoni_drop_translation/graphs/COM U Velocity New TimeB.png")
plt.show()