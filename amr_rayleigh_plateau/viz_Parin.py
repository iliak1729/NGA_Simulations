import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv

# ============================================================
# Parameters
# ============================================================

# ligament parameters
r0 = 0.1        # initial radius of the jet
#eps = 0.1       # perturbation amplitude
eps = 0.01

# other parameters
Lz = 1.0        # axial domain length
mode = 1        # perturbation mode

# liquid parameters
rho = 1000.0    # density, kg/m^3
# sigma = 2e-7   # surface tension, N/m
sigma = 0.072   # surface tension, N/m
mu = 0.001      # dynamic viscosity, Pa.s

nu = mu / rho

# Capillary time scale
t_cap = np.sqrt(rho * r0**3 / sigma)
print(f"Capillary time scale: {t_cap:.6e} s")
# Laplace and Ohnesorge numbers
La = sigma * rho * r0 / mu**2
Oh = 1.0 / np.sqrt(La)

# Wavenumber
k = 2.0 * np.pi * mode / Lz
kR = k * r0

print(f"Capillary time scale: {t_cap:.6e} s")
print(f"Laplace number: {La:.6e}")
print(f"Ohnesorge number: {Oh:.6e}")
print(f"Dimensional wavenumber k: {k:.6e}")
print(f"Non-dimensional wavenumber kR: {kR:.6e}")

# ============================================================
# Growth rates
# ============================================================

omega_rayleigh = np.sqrt(sigma / (rho * r0**3)* kR* (1.0 - kR**2)* iv(1, kR) / iv(0, kR))

A = 3.0 * nu * k**2
B = sigma / (2.0 * rho * r0**3) * (1.0 - kR**2) * kR**2
omega_weber = -0.5 * A + np.sqrt((0.5 * A)**2 + B)


print(f"Rayleigh growth rate omega_R: {omega_rayleigh:.6e} 1/s")
print(f"Weber growth rate omega_W:    {omega_weber:.6e} 1/s")
print(f"Non-dimensional Rayleigh growth rate omega_R*t_cap: {omega_rayleigh*t_cap:.6e}")
print(f"Non-dimensional Weber growth rate omega_W*t_cap:    {omega_weber*t_cap:.6e}")


# time
tstar = np.linspace(0.0, 18, 500)   # nondimensional time t / t_cap
t = tstar * t_cap                     # dimensional time

# ============================================================
# Linear growth rate
# ============================================================

amp_rayleigh = np.exp(omega_rayleigh * t)
amp_weber = np.exp(omega_weber * t)

# Radius extrema from linear theory
rmax_rayleigh = r0 * (1.0 + eps * amp_rayleigh)
rmin_rayleigh = r0 * (1.0 - eps * amp_rayleigh)

rmax_weber = r0 * (1.0 + eps * amp_weber)
rmin_weber = r0 * (1.0 - eps * amp_weber)

# Normalized quantities used in paper
rmax_rayleigh_norm = (rmax_rayleigh - r0) / (r0 * eps)
rmin_rayleigh_norm = (r0 - rmin_rayleigh) / (r0 * eps)

rmax_weber_norm = (rmax_weber - r0) / (r0 * eps)
rmin_weber_norm = (r0 - rmin_weber) / (r0 * eps)

# ============================================================
# Sim data
# ============================================================

# List of simulation files and labels
sim_cases = [
    # (
    #     "amr_rayleigh_plateau/Baseline_Runs_Level4_PCellForm_SameFluid/amr_PCST_NoP/monitor/rayleigh_plateau",
    #     "AMR,PCentered,NoP,PLIC"
    # ),
    # (
    #     "amr_rayleigh_plateau/Baseline_Runs_Level4_PCellForm_SameFluid/amr_PCST_NoP_PPIC/monitor/rayleigh_plateau",
    #     "AMR,PCentered,NoP,PPIC"
    # ),
    # (
    #     "amr_rayleigh_plateau/Baseline_Runs_Level4_PCellForm_SameFluid/amr_PCST_NoP_Laplacian/monitor/rayleigh_plateau",
    #     "PCentered,NoP,PLIC,Laplacian"
    # ),
    # (
    #     "amr_rayleigh_plateau/Baseline_Runs_Level4_PCellForm_SameFluid/amr_PCST_NoP_PPIC_Laplacian/monitor/rayleigh_plateau",
    #     "PCentered,NoP,PPIC,Laplacian"
    # ),
    # (
    #     "ThreeD-rayleigh-plateau/NONAMR_RUNS/NONAMR_Cell14D_NoP/monitor/rayleigh_plateau",
    #     "Non-AMR,NoP,PLIC"
    # ),
    # (
    #     "ThreeD-rayleigh-plateau/NONAMR_RUNS/NONAMR_Cell14D_NoP_Laplacian/monitor/rayleigh_plateau",
    #     "Non-AMR,NoP,PLIC,Laplacian"
    # ),
    # (
    #     "ThreeD-rayleigh-plateau/NONAMR_RUNS/NONAMR_Cell14D_NoP_PPIC/monitor/rayleigh_plateau",
    #     "Non-AMR,NoP,PPIC"
    # ),
    # (
    #     "ThreeD-rayleigh-plateau/NONAMR_RUNS/NONAMR_Cell14D_NoP_PPIC_Laplacian/monitor/rayleigh_plateau",
    #     "Non-AMR,NoP,PPIC,Laplacian"
    # ),
    # (
    #     "ThreeD-rayleigh-plateau/monitor/rayleigh_plateau",
    #     "LIVE,NoP,PU,New"
    # ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/Debugging_PCST_P_Matching/monitor/rayleigh_plateau",
        "NoP,Matching"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/FaceParaboloid_Cell14D_PCST_NoP_NoMatching/monitor/rayleigh_plateau",
        "NoP,NoMatching"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/FaceParaboloid_Cell14D_PCST_P_Matching/monitor/rayleigh_plateau",
        "P,Matching"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/FaceParaboloid_Cell14D_PCST_P_NoMatching/monitor/rayleigh_plateau",
        "P,NoMatching"
    ),
    (
        "ThreeD-rayleigh-plateau/monitor/rayleigh_plateau",
        "LIVE,RP"
    ),

    # (
    #     "ThreeD-rayleigh-plateau/NONAMR_RUNS/NONAMR_Paraboloid_JibbenPU_Cell14D_NoP_PPIC/monitor/rayleigh_plateau",
    #     "NoP,ParaboloidForm,JibbenPU"
    # ),
    # (
    #     "ThreeD-rayleigh-plateau/NONAMR_RUNS/NONAMR_Paraboloid_PUppic_Cell14D_NoP_PPIC/monitor/rayleigh_plateau",
    #     "NoP,ParaboloidForm,PUppic"
    # ),



]

# Store simulation data here
sim_data = []
def safe_float(x):
    try:
        return float(x)
    except ValueError:
        return np.nan
    
for filename, label in sim_cases:

    data = np.loadtxt(filename, skiprows=2, converters=safe_float)

    step_data = data[:, 0].astype(int)
    tstar_data = data[:, 2]
    rmax_data = data[:, 3]
    rmin_data = data[:, 4]

    # Remove invalid points
    valid = (rmax_data > 0.0) & (rmin_data > 0.0)

    step_data = step_data[valid]
    tstar_data = tstar_data[valid]
    rmax_data = rmax_data[valid]
    rmin_data = rmin_data[valid]

    # Normalize
    rmax_data_norm = (rmax_data - r0) / (r0 * eps)
    rmin_data_norm = (r0 - rmin_data) / (r0 * eps)

    sim_data.append({
        "label": label,
        "step": step_data,
        "tstar": tstar_data,
        "rmax": rmax_data_norm,
        "rmin": rmin_data_norm,
    })


# ============================================================
# Max radius
# ============================================================

plt.figure(figsize=(7, 5))

# Analytical solutions
plt.plot(
    tstar,
    rmax_rayleigh_norm,
    label="Rayleigh"
)

plt.plot(
    tstar,
    rmax_weber_norm,
    "--",
    label="Weber"
)

# Simulation results
for case in sim_data:
    plt.plot(
        case["tstar"],
        case["rmax"],
        "o-",
        markersize=4,
        label=case["label"]
    )

plt.xlabel(r"Non-dimensional time, $t/t_{\mathrm{cap}}$")
plt.ylabel(r"$(r_{\max}-r_0)/(r_0\epsilon)$")

plt.yscale("log")
plt.grid(True, which="both")
plt.legend()

# plt.ylim(1, 10)
plt.xlim(0, 5)
plt.title('Maximum radius growth comparison')
plt.show()


# ============================================================
# Min radius
# ============================================================

plt.figure(figsize=(7, 5))

# Analytical solutions
plt.plot(
    tstar,
    rmin_rayleigh_norm,
    label="Rayleigh"
)

plt.plot(
    tstar,
    rmin_weber_norm,
    "--",
    label="Weber"
)

# Simulation results
for case in sim_data:
    plt.plot(
        case["tstar"],
        case["rmin"],
        "o",
        markersize=4,
        label=case["label"]
    )

plt.xlabel(r"Non-dimensional time, $t/t_{\mathrm{cap}}$")
plt.ylabel(r"$(r_0-r_{\min})/(r_0\epsilon)$")

plt.yscale("log")
plt.grid(True, which="both")
plt.legend()

# plt.ylim(1, 10)
plt.xlim(0, 5)
plt.title('Minimum radius growth comparison')
plt.show()