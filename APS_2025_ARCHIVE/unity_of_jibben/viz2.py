import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special


def Prosperti1981(t):
    # Get Prospereti Data
    rho_l = 1
    rho_u = 1
    beta = (rho_l*rho_u)/((rho_l+rho_u)*(rho_l+rho_u))
    g = 0
    zeta = 1
    mu = 0.0182571749236
    a0 = 0.01
    u0 = 0
    nu = mu/rho_u

    lambdaX = 1
    kx = 2*np.pi/lambdaX
    ky = 0 
    k = np.sqrt(kx*kx+ky*ky)

    omega0 = (rho_l - rho_u)*g*k/(rho_l+rho_u) + zeta*k*k*k/(rho_l+rho_u)
    # print(omega0)
    omega0 = np.sqrt(omega0)

    a = 1
    b = -4*beta*np.sqrt(k*k*nu)
    c = 2*(1-6*beta)*k*k*nu
    d = 4*(1-3*beta)*np.pow(k*k*nu,3/2)
    e = (1-4*beta)*nu*nu*k*k*k*k + omega0*omega0

    coeffArray = [a,b,c,d,e]
    z = np.roots(coeffArray)
    Z1 = (z[1]-z[0])*(z[2]-z[0])*(z[3]-z[0])
    Z2 = (z[0]-z[1])*(z[2]-z[1])*(z[3]-z[1])
    Z3 = (z[1]-z[2])*(z[0]-z[2])*(z[3]-z[2])
    Z4 = (z[1]-z[3])*(z[2]-z[3])*(z[0]-z[3])
    Z = [Z1,Z2,Z3,Z4]
    # print(Z)
    numer = 4*(1-4*beta)*nu*nu*k*k*k*k 
    denom = 8*(1-4*beta)*nu*nu*k*k*k*k + omega0*omega0
    C1 = numer/denom
    # print(z)
    a = np.linspace(0,0,len(t))
    for j in range(len(t)):  
        a[j] = C1 *a0 * np.sqrt(special.erfc(nu*k*k*t[j]))
        for i in range(4):
            inside = (omega0*omega0*a0)/(z[i]*z[i]- nu*k*k)
            C = z[i]*inside/Z[i]
            # val = np.exp((z[i]*z[i]-nu*k*k)*t[j])*math.erfc(np.real(z[i])*math.sqrt(t[j]))
            val = np.exp((z[i]*z[i]-nu*k*k)*t[j])*special.erfc(z[i]*np.sqrt(t[j]))
            a[j] = a[j] + np.real(C*val)
    return(omega0,a)

omega0 = 1.11367E+01

# CSF PLIC
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF8_NewDt_3L"
dataCSF8 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF16_NewDt_3L"
dataCSF16 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF32_NewDt_3L"
dataCSF32 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF64_NewDt_3L"
dataCSF64 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF128_NewDt_3L"
dataCSF128 = np.loadtxt(datafile, skiprows=2)

# CSS PLIC
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS8_NewDt_3L"
dataCSS8 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS16_NewDt_3L"
dataCSS16 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS32_NewDt_3L"
dataCSS32 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS64_NewDt_3L"
dataCSS64 = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS128_NewDt_3L"
dataCSS128 = np.loadtxt(datafile, skiprows=2)

# CSF Jibben
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF8_Jibben_NewDt_3L"
dataCSF8_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF16_Jibben_NewDt_3L"
dataCSF16_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF32_Jibben_NewDt_3L"
dataCSF32_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF64_Jibben_NewDt_3L"
dataCSF64_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSF128_Jibben_NewDt_3L"
dataCSF128_Jibben = np.loadtxt(datafile, skiprows=2)

# CSS Jibben
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS8_Jibben_NewDt_3L"
dataCSS8_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS16_Jibben_NewDt_3L"
dataCSS16_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS32_Jibben_NewDt_3L"
dataCSS32_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS64_Jibben_NewDt_3L"
dataCSS64_Jibben = np.loadtxt(datafile, skiprows=2)
datafile = "unity_of_jibben/monitor/PU_CapWave_CSS128_Jibben_NewDt_3L"
dataCSS128_Jibben = np.loadtxt(datafile, skiprows=2)


# Exact Data
tExact = np.linspace(0,25/omega0,100)
(omega0,aExact) = Prosperti1981(tExact)
plt.figure()



plt.plot(dataCSS64[:,6],dataCSS64[:,4],label = "CSS64",linewidth = 8)
plt.plot(dataCSF64[:,6],dataCSF64[:,4],label = "CSF64",linewidth = 3) 
plt.plot(tExact*omega0,np.abs(aExact),'ko',label = 'Exact',linewidth=1)

plt.xlabel("$\omega_0 t$")
plt.ylabel("Amplitude")
plt.legend()


datas = [dataCSF8,
dataCSF16,
dataCSF32,
dataCSF64,
dataCSF128, 
dataCSS8,
dataCSS16,
dataCSS32,
dataCSS64,
dataCSS128,
dataCSF8_Jibben,
dataCSF16_Jibben,
dataCSF32_Jibben,
dataCSF64_Jibben,
dataCSF128_Jibben, 
dataCSS8_Jibben,
dataCSS16_Jibben,
dataCSS32_Jibben,
dataCSS64_Jibben,
dataCSS128_Jibben]


errors = np.linspace(0,0,len(datas))
NValues = np.array([8,16,32,64,128])
# Error Calculation
tTest = np.linspace(0,25/omega0,1000)
(omega0,aExactTest) = Prosperti1981(tTest)
aExactTest = np.abs(aExactTest)
for i in range(len(errors)):
    error = 0.0
    print(i)
    # Get Data
    data = datas[i]
    # Get Values
    omegaT = data[:,6]
    aList = data[:,4]
    t = omegaT/omega0
    # (o,aIn) = Prosperti1981(t)
    aList = np.array(aList)
    (omega0,aExact) = Prosperti1981(t)
    aExact = np.abs(aExact)
    # Next, we will interpolate values to new array
    aTest = np.interp(tTest,t,aList)

    # Calculate Error
    error = (aList-aExact)*(aList-aExact)
    # error = (aTest-aExactTest)*(aTest-aExactTest)
    # Get Value
    error = sum(error)/len(error)
    error = error
    errors[i] = np.sqrt(error)

fig,ax  = plt.subplots(1,1)
ax.plot(NValues,1.5*errors[0]*NValues[0]**2/NValues**2,label = 'O(N^2)',linewidth = 3,color = 'k')
# ax.plot(NValues,0.5*errors[0]*NValues[0]**1/NValues**1,label = 'O(N)',linewidth = 3,color = 'k')
print("ERROR DATA ================ ")
ax.plot(NValues,errors[0:5],'o-',label = 'CSF PLIC',linewidth = 3)
print("CSF PLIC ERRORS: ",errors[0:5])
ax.plot(NValues,errors[5:10],'o--',label = 'CSS PLIC',linewidth = 3)
print("CSS PLIC ERRORS: ",errors[5:10])
ax.plot(NValues,errors[10:15],'o-',label = 'CSF Jibben',linewidth = 3)
print("CSF Jibben ERRORS: ",errors[10:15])
ax.plot(NValues,errors[15:20],'o--',label = 'CSS Jibben',linewidth = 3)
print("CSS Jibben ERRORS: ",errors[15:20])

ax.set_xscale("log", base=2)
ax.set_yscale("log")
# plt.plot(64,errors[20],'kx',label = 'Larger Domain')
plt.xlabel("$N/\lambda$")
plt.ylabel("$\epsilon_{L2}$")
plt.grid(True)
plt.legend()
plt.show()