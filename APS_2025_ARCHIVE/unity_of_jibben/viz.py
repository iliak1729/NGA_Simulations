import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special

def plot_simulation_numpy(file_path, x_col, y_col,nor_x,nor_y,name):
    """
    Plot columns from the simulation data file using NumPy.

    Parameters
    ----------
    file_path : str
        Path to the simulation data file.
    x_col : str
        Column name for the x-axis.
    y_col : str
        Column name for the y-axis.
    """
    # Column names in the file (after skipping the 2 header rows)
    col_names = ["Timestep", "Time", "dt", "ViscousTS", "CapillaryTS",
                 "Vrms", "VOF_max", "VOF_min", "VOF_integral","COMX","COMY","COMXZ","Wave Amplitude","Wave Frequency","Wave Time"]

    # Load numeric data
    data = np.loadtxt(file_path, skiprows=2)

    # Map names to column indices
    name_to_idx = {name: i for i, name in enumerate(col_names)}

    if x_col not in name_to_idx or y_col not in name_to_idx:
        raise ValueError(f"Columns must be one of: {col_names}")

    x = data[:, name_to_idx[x_col]]
    y = data[:, name_to_idx[y_col]]
    x = x / nor_x
    y = y / nor_y
    # Plot
    
    plt.plot(x, y, linestyle="-", markersize=3,label = name,linewidth=5)
    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.grid(True)
    return(x,y)
    

def plot_COM_simulation_numpy(file_path, x_col, y_col,nor_x,nor_y,name):
    """
    Plot columns from the simulation data file using NumPy.

    Parameters
    ----------
    file_path : str
        Path to the simulation data file.
    x_col : str
        Column name for the x-axis.
    y_col : str
        Column name for the y-axis.
    """
    # Column names in the file (after skipping the 2 header rows)
    col_names = ["Timestep", "Time", "dt", "ViscousTS", "CapillaryTS",
                 "Vrms", "VOF_max", "VOF_min", "VOF_integral","COMX","COMY","COMXZ","Wave Amplitude","Wave Frequency","Wave Time"]

    # Load numeric data
    data = np.loadtxt(file_path, skiprows=2)

    # Map names to column indices
    name_to_idx = {name: i for i, name in enumerate(col_names)}

    if x_col not in name_to_idx or y_col not in name_to_idx:
        raise ValueError(f"Columns must be one of: {col_names}")

    x = data[:, name_to_idx[x_col]]
    y = data[:, name_to_idx[y_col]]
    x = x / nor_x
    y = y / nor_y
    
    # Plot
    return x,y

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
plt.figure()
(omegaTCss128,aCSS128)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSWithP_128","Wave Time","Wave Amplitude",1.0,1.0,"CSS128P")
(omegaTCss64,aCSS64)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSWithP_64","Wave Time","Wave Amplitude",1.0,1.0,"CSS64P")
(omegaTCss32,aCSS32)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSWithP_32","Wave Time","Wave Amplitude",1.0,1.0,"CSS32P")
(omegaTCss16,aCSS16)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSWithP_16","Wave Time","Wave Amplitude",1.0,1.0,"CSS16P")

(omegaTCss128b,aCSS128b)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSNoP_128","Wave Time","Wave Amplitude",1.0,1.0,"CSS128")
(omegaTCss64b,aCSS64b)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSNoP_64","Wave Time","Wave Amplitude",1.0,1.0,"CSS64")
(omegaTCss32b,aCSS32b)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSNoP_32","Wave Time","Wave Amplitude",1.0,1.0,"CSS32")
(omegaTCss16b,aCSS16b)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSNoP_16","Wave Time","Wave Amplitude",1.0,1.0,"CSS16")

(omegaTCsf128,aCSF128)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_128","Wave Time","Wave Amplitude",1.0,1.0,"CSF128")
(omegaTCsf64,aCSF64)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_64","Wave Time","Wave Amplitude",1.0,1.0,"CSF64")
(omegaTCsf32,aCSF32)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_32","Wave Time","Wave Amplitude",1.0,1.0,"CSF32")
(omegaTCsf16,aCSF16)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_16","Wave Time","Wave Amplitude",1.0,1.0,"CSF16")

(omegaTCsf128_HF,aCSF128_HF)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_128_HF","Wave Time","Wave Amplitude",1.0,1.0,"CSF128")
(omegaTCsf64_HF,aCSF64_HF)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_64_HF","Wave Time","Wave Amplitude",1.0,1.0,"CSF64")
(omegaTCsf32_HF,aCSF32_HF)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_32_HF","Wave Time","Wave Amplitude",1.0,1.0,"CSF32")
(omegaTCsf16_HF,aCSF16_HF)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_16_HF","Wave Time","Wave Amplitude",1.0,1.0,"CSF16")

(omegaTCss16c,aCSS16c)=plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSSNoP_16_CFL01","Wave Time","Wave Amplitude",1.0,1.0,"CSF16_01")

datafile = "unity_of_jibben/monitor/PUJibbenTesting_CapillaryWave_Testing"
data1 = np.loadtxt(datafile, skiprows=2)

ValuePairs = [(omegaTCss128,aCSS128),(omegaTCss64,aCSS64),(omegaTCss32,aCSS32),(omegaTCss16,aCSS16),
              (omegaTCsf128,aCSF128),(omegaTCsf64,aCSF64),(omegaTCsf32,aCSF32),(omegaTCsf16,aCSF16),
              (omegaTCss128b,aCSS128b),(omegaTCss64b,aCSS64b),(omegaTCss32b,aCSS32b),(omegaTCss16b,aCSS16b),
              (omegaTCsf128_HF,aCSF128_HF),(omegaTCsf64_HF,aCSF64_HF),(omegaTCsf32_HF,aCSF32_HF),(omegaTCsf16_HF,aCSF16_HF),
              (omegaTCss16c,aCSS16c)]
NValues = np.array([128,64,32,16])
plt.figure()
# Calculate Errors
errors = np.linspace(0,0,len(ValuePairs))
# Times to calculate values
length = min(len(aCSS128),len(aCSS64),len(aCSS32),len(aCSS16))
length = 738
omegaTTest = np.linspace(0,25,length)
tTest = omegaTTest/omega0
(o,aExact) = Prosperti1981(tTest)
aExact = np.abs(aExact)
for i in range(len(ValuePairs)):
    print("N = ",NValues[i%4])
    error = 0.0

    # Get Pair
    pair = ValuePairs[i]
    # Get Values
    omegaT = pair[0]
    aList = pair[1]
    print(np.shape(omegaT))
    t = omegaT/omega0
    (o,aIn) = Prosperti1981(t)
    aList = np.array(aList)

    # Next, we will interpolate values to new array
    aTest = np.interp(tTest,t,aList)
    aIn = np.abs(aIn)
    # Now, we have a constant length array, so we will use this to do the error
    error = (aTest-aExact)*(aTest-aExact)
    # error = (aList-aIn)*(aList-aIn)
    # if i < 4:
    #     plt.figure()
    #     plt.plot(omegaT,np.sqrt(error))
    #     plt.figure()
    #     plt.plot(omegaT,aIn)
    #     plt.plot(omegaT,aList)
    #     plt.title("N = "+str(NValues[i%4]))
    print(len(error))
    error = sum(error)/len(error)
    error = error*omega0/25
    errors[i] = np.sqrt(error)

plt.figure()
plt.loglog(NValues,errors[0:4],'-x',label = "CSS-P")
plt.loglog(NValues,errors[4:8],'-x',label = "CSF")
plt.loglog(NValues,errors[12:16],'-x',label = "CSF_HF")
plt.loglog(NValues[1:4],errors[9:12],'-x',label = "CSS-NoP")
O = 2
plt.loglog(NValues,errors[11]*(NValues[3]**O)/(NValues**O),label = "O(N^-"+str(O)+")")
plt.grid(True)
plt.legend()

# Visual Compare
plt.figure(figsize=(8,5))
omega0 = 1.11367E+01
t = np.linspace(0,25/omega0,10000)
(omega0,a) = Prosperti1981(t)



# plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave","Wave Time","Wave Amplitude",1.0,1.0,"CSS")
plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_128","Wave Time","Wave Amplitude",1.0,1.0,"CSF128")
plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_16","Wave Time","Wave Amplitude",1.0,1.0,"CSF16")
# plot_simulation_numpy("capillary_wave/monitor/PUTesting_CapillaryWave_CSF_128","Wave Time","Wave Amplitude",1.0,1.0,"CSF128")

plt.xlabel("omega*T")
plt.ylabel("a")

plt.title("Capillary Wave - CSS")
dataX = [1.8441308531109686,
         3.383579217447082,
5.339961513790892,
6.847338037203336,
8.723540731237973,
10.327132777421424,
12.235407312379731,
13.75881975625401,
15.68313021167415,
17.17447081462476,
19.146889031430405,
20.654265554842848,
22.530468248877487,
24.166132135984608]

dataY = [0.00006642066420664207,
0.0072619926199262,
0.0001033210332103321,
0.004996309963099631,
0.00016236162361623618,
0.0034686346863468634,
0.00005904059040590406,
0.0024206642066420666,
0.00005166051660516605,
0.0016383763837638377,
0.00005166051660516605,
0.0012029520295202953,
0.000007380073800738008,
0.0008413284132841328]
# plt.plot(dataX,dataY,'kx',label='Saini Plot Digitized')

plt.plot(t*omega0,np.abs(a),label = 'Prosperti (1981)',linewidth=3)

plt.xlim([0,25])
plt.legend()
# plt.savefig('capillary_wave/CapillaryWave_CompareToSaini.png')

plt.figure()
plt.plot(data1[:,14],data1[:,12],label = "Nz=3",linewidth=3)
# plt.plot(data1[:,14],data1[:,12]-0.09434,label = "Nz=3 Translated",linewidth=3)
plt.plot(t*omega0,np.abs(a),'--',label = 'Prosperti (1981)',linewidth=3)
plt.xlim([0,25])
plt.legend()
plt.show()

