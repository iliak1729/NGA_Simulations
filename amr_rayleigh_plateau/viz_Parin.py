import numpy as np
import matplotlib.pyplot as plt
from scipy.special import iv

# ============================================================
# Parameters
# ============================================================

# ligament parameters
r0 = 0.1        # initial radius of the jet
#eps = 0.1       # perturbation amplitude
eps = 0.02

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
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/SliceBased_CSF_Jibben/monitor/rayleigh_plateau",
        "CSF,Jibben"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/SliceBased_CSF_LVIRA/monitor/rayleigh_plateau",
        "CSF,LVIRA"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/SliceBased_CSF_PLICNET/monitor/rayleigh_plateau",
        "CSF,PLICNET"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/SliceBased_PCST_P_Matching_Jibben/monitor/rayleigh_plateau",
        "PCST,P,M,Jibben"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/SliceBased_PCST_P_Matching_LVIRA/monitor/rayleigh_plateau",
        "PCST,P,M,LVIRA"
    ),
    (
        "ThreeD-rayleigh-plateau/FACE_PARABOLOID/SliceBased_PCST_P_Matching_PLICNET/monitor/rayleigh_plateau",
        "PCST,P,M,PLICNET"
    ),



    # (
    #     "ThreeD-rayleigh-plateau/monitor/rayleigh_plateau",
    #     "LIVE,RP"
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
# plt.xlim(0, 5)
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
# plt.xlim(0, 5)
plt.title('Minimum radius growth comparison')
plt.show()

# ============================================================
# VOF Mass Conservation Diagnostics
# ============================================================

# ============================================================
# VOF Mass Conservation Diagnostics - All Locations
# ============================================================

# Grid parameters
Nx, Ny, Nz = 64, 64, 64
Lx, Ly, Lz = 1.0, 1.0, 1.0

cell_volume = (Lx / Nx) * (Ly / Ny) * (Lz / Nz)

# Select which simulation to examine
case_index = 3 # 0=CSF Jibben, 3=PCST Jibben, etc.

filename, label = sim_cases[case_index]

data = np.loadtxt(
    filename,
    skiprows=2,
    converters=safe_float
)

# Time and timestep
steps = data[:, 0]
tstar_data = data[:, 2]

# VFint at Locations 1-5 (columns 9-13)
VF_locations = data[:, 9:14]

# Remove invalid/uninitialized data
valid = (
    np.all(np.isfinite(VF_locations), axis=1)
    & np.all(VF_locations > 0.0, axis=1)
)

steps = steps[valid]
tstar_data = tstar_data[valid]
VF_locations = VF_locations[valid]

# ============================================================
# Figure 1: Volume change between consecutive locations
# ============================================================

plt.figure(figsize=(9, 6))

location_labels = [
    "Location 1 → 2 (VOF Advection)",
    "Location 2 → 3 (dMOMdt)",
    "Location 3 → 4 (solveImplicit+Updates)",
    "Location 4 → 5 (add_surface_Tension_jump)"
]

for i in range(4):

    delta_V = VF_locations[:, i+1] - VF_locations[:, i]

    # Normalize by the volume of one grid cell
    delta_cells = delta_V / cell_volume

    plt.plot(
        tstar_data,
        delta_cells,
        "o-",
        markersize=3,
        label=location_labels[i]
    )

plt.axhline(0.0, color="black", linestyle="--", linewidth=1)

plt.xlabel(r"Non-dimensional time, $t/t_{\mathrm{cap}}$")
plt.ylabel(r"Volume change, $\Delta V/\Delta V_{\mathrm{cell}}$")

plt.title(f"Liquid Volume Change Between Locations\n{label}")

plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()


# ============================================================
# Figure 2: Cumulative volume error at all five locations
# ============================================================

plt.figure(figsize=(9, 6))

# Initial volume: first valid Location 1 value
V0 = VF_locations[0, 0]

for i in range(5):

    relative_error = (VF_locations[:, i] - V0) / V0

    plt.plot(
        tstar_data,
        relative_error,
        "o-",
        markersize=3,
        label=f"Location {i+1}"
    )

plt.axhline(0.0, color="black", linestyle="--", linewidth=1)

plt.xlabel(r"Non-dimensional time, $t/t_{\mathrm{cap}}$")
plt.ylabel(r"Relative volume error, $(V_L-V_{L,0})/V_{L,0}$")

plt.title(f"Cumulative Liquid Volume Error\n{label}")

plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# ============================================================
# VOF Mass Conservation Diagnostics
# ============================================================

# Grid parameters
Nx = 64
Ny = 64
Nz = 64

Lx = 1.0
Ly = 1.0
Lz = 1.0

cell_volume = (Lx / Nx) * (Ly / Ny) * (Lz / Nz)

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(9, 8),
    sharex=True
)

for filename, label in sim_cases:

    data = np.loadtxt(
        filename,
        skiprows=2,
        converters=safe_float
    )

    # Columns:
    # 0: timestep
    # 1: time
    # 2: t/t_cap
    # ...
    # 9:  VFint Location 1
    # 10: VFint Location 2
    # 11: VFint Location 3
    # 12: VFint Location 4
    # 13: VFint Location 5

    tstar_data = data[:, 2]

    VF_loc1 = data[:, 9]
    VF_loc2 = data[:, 10]

    # Ignore uninitialized monitoring values
    valid = (
        np.isfinite(VF_loc1)
        & np.isfinite(VF_loc2)
        & (VF_loc1 > 0.0)
        & (VF_loc2 > 0.0)
    )

    tstar_valid = tstar_data[valid]
    VF_loc1 = VF_loc1[valid]
    VF_loc2 = VF_loc2[valid]

    if len(VF_loc1) == 0:
        continue

    # --------------------------------------------------------
    # Per-timestep volume change
    # --------------------------------------------------------

    delta_V = VF_loc2 - VF_loc1

    # Express as equivalent number of grid-cell volumes
    delta_cells = delta_V / cell_volume

    ax1.plot(
        tstar_valid,
        delta_cells,
        "o-",
        markersize=3,
        label=label
    )

    # --------------------------------------------------------
    # Cumulative relative volume error
    # --------------------------------------------------------

    V0 = VF_loc1[0]

    relative_error = (VF_loc2 - V0) / V0

    ax2.plot(
        tstar_valid,
        relative_error,
        "o-",
        markersize=3,
        label=label
    )


# ------------------------------------------------------------
# Top: Volume gained per timestep
# ------------------------------------------------------------

ax1.axhline(0.0, color="black", linestyle="--", linewidth=1)

ax1.set_ylabel(
    r"Volume change per step, $\Delta V/\Delta V_{\mathrm{cell}}$"
)

ax1.set_title("VOF Mass Conservation Diagnostics")

ax1.grid(True, which="both", alpha=0.3)
ax1.legend()


# ------------------------------------------------------------
# Bottom: Cumulative relative volume error
# ------------------------------------------------------------

ax2.axhline(0.0, color="black", linestyle="--", linewidth=1)

ax2.set_xlabel(
    r"Non-dimensional time, $t/t_{\mathrm{cap}}$"
)

ax2.set_ylabel(
    r"Relative volume error, $(V_L-V_{L,0})/V_{L,0}$"
)

ax2.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.show()