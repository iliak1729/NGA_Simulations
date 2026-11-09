import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
# ============================================================
# Colors
# ============================================================

color_black = '#000000'
color_orange = '#E69F00'
color_skyblue = '#56B4E9'
color_teal = '#009E73'
color_yellow = '#F0E442'
color_blue = '#0072B2'
color_red = '#D55E00'
color_pink = '#CC79A7'


# ============================================================
# Physics
# ============================================================

sigma = 1.0
rho = 1.0
D = 1
La = 1200

mu = np.sqrt(sigma * rho * D / La)

Usigma = sigma / mu
Usigma = 1
Tmu = rho * D**2 / mu
print("mu = ", mu)
print("Usigma = ", Usigma)
print("Tmu = ", Tmu)
print("Tmax = ", Tmu)

# ============================================================
# Smoothing Options
# ============================================================

APPLY_SMOOTHING = True
SMOOTHING_WINDOW = 50

def smooth_data(y):
    if not APPLY_SMOOTHING or SMOOTHING_WINDOW <= 1:
        return y

    window = min(SMOOTHING_WINDOW, len(y))

    # Centered moving average
    kernel = np.ones(window) / window

    # Pad edges to avoid artificial drops near boundaries
    left = (window - 1) // 2
    right = window - 1 - left

    y_padded = np.pad(y, (left, right), mode="edge")

    return np.convolve(y_padded, kernel, mode="valid")
# ============================================================
# Simulation Data
# ============================================================

# List of simulation files, labels, and colors
# List of simulation files, labels, colors, and line types
sim_cases = [
    # ============================================================
    # D8 Cases - Original (Solid)
    # ============================================================
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200/StaticDrop_CSF_D8_SingleFluid/monitor/simulation",
    #     "CSF,D8",
    #     "xkcd:charcoal",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200/StaticDrop_PCST_D8_NoP_NoMatching_SingleFluid/monitor/simulation",
    #     "NoP,NoM,D8",
    #     "xkcd:azure",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200/StaticDrop_PCST_D8_NoP_Matching_SingleFluid/monitor/simulation",
    #     "NoP,M,D8",
    #     "xkcd:burnt orange",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200/StaticDrop_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "P,NoM,D8",
    #     "xkcd:red",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200/StaticDrop_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "P,M,D8",
    #     "xkcd:sienna",
    #     "-"
    # ),

    # ============================================================
    # D16 Cases - WallsBC (Dashed)
    # ============================================================
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_CSF_D16_SingleFluid/monitor/simulation",
    #     "CSF,D16,WallsBC",
    #     "xkcd:charcoal",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D16_NoP_NoMatching_SingleFluid/monitor/simulation",
    #     "NoP,NoM,D16,WallsBC",
    #     "xkcd:azure",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D16_NoP_Matching_SingleFluid/monitor/simulation",
    #     "NoP,M,D16,WallsBC",
    #     "xkcd:burnt orange",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D16_P_NoMatching_SingleFluid/monitor/simulation",
    #     "P,NoM,D16,WallsBC",
    #     "xkcd:teal",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D16_P_Matching_SingleFluid/monitor/simulation",
    #     "P,M,D16,WallsBC",
    #     "xkcd:purple",
    #     "--"
    # ),
    # ============================================================
    # D8 Cases - WallsBC (Dash-dot)
    # ============================================================
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_CSF_D8_SingleFluid/monitor/simulation",
    #     "CSF,D8,WallsBC",
    #     "xkcd:charcoal",
    #     "-."
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D8_NoP_NoMatching_SingleFluid/monitor/simulation",
    #     "NoP,NoM,D8,WallsBC",
    #     "xkcd:azure",
    #     "-."
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D8_NoP_Matching_SingleFluid/monitor/simulation",
    #     "NoP,M,D8,WallsBC",
    #     "xkcd:burnt orange",
    #     "-."
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "P,NoM,D8,WallsBC",
    #     "xkcd:teal",
    #     "-."
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC/StaticDrop_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "P,M,D8,WallsBC",
    #     "xkcd:purple",
    #     "-."
    # ),
    # ============================================================
    # D8 Cases - Refined Timestep (Dotted)
    # ============================================================
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_RefinedTimestep/StaticDrop_CSF_D8_SingleFluid/monitor/simulation",
    #     "CSF,D8,RefinedTime",
    #     "xkcd:charcoal",
    #     ":"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_RefinedTimestep/StaticDrop_PCST_D8_NoP_NoMatching_SingleFluid/monitor/simulation",
    #     "NoP,NoM,D8,RefinedTime",
    #     "xkcd:azure",
    #     ":"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_RefinedTimestep/StaticDrop_PCST_D8_NoP_Matching_SingleFluid/monitor/simulation",
    #     "NoP,M,D8,RefinedTime",
    #     "xkcd:burnt orange",
    #     ":"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC_RefinedTime/StaticDrop_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "P,NoM,D8,RefinedTime",
    #     "xkcd:red",
    #     ":"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_WallsBC_RefinedTime/StaticDrop_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "P,M,D8,RefinedTime",
    #     "xkcd:sienna",
    #     ":"
    # ),
    #     # ============================================================
    # # D8 Cases - Normal Weighting Comparison (Solid)
    # # ============================================================
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_DistanceSqWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "Distance Squared Weight",
    #     "xkcd:blue",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_DistanceWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "Distance Weight",
    #     "xkcd:orange",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_NoWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "No Weight",
    #     "xkcd:charcoal",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_SDWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "Surface Density Weight",
    #     "xkcd:teal",
    #     "-"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_VFWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
    #     "Volume Fraction Weight",
    #     "xkcd:purple",
    #     "-"
    # ),
    #     # ============================================================
    # # D8 Cases - Normal Weighting Comparison - Matching (Dashed)
    # # ============================================================
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_DistanceSqWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "Distance Squared Weight, Matching",
    #     "xkcd:blue",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_DistanceWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "Distance Weight, Matching",
    #     "xkcd:orange",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_NoWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "No Weight, Matching",
    #     "xkcd:charcoal",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_SDWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "Surface Density Weight, Matching",
    #     "xkcd:teal",
    #     "--"
    # ),
    # (
    #     "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/StaticDrop_VFWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
    #     "Volume Fraction Weight, Matching",
    #     "xkcd:purple",
    #     "--"
    # ),
        # ============================================================
    # D8 Cases - Translation Drop - Normal Weighting Comparison
    # Solid = No Matching | Dashed = Matching
    # ============================================================
    (
        "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/Translation_DistanceWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
        "Distance Weight, No Matching",
        "xkcd:orange",
        "-"
    ),
    (
        "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/Translation_DistanceSqWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
        "Distance Squared Weight, No Matching",
        "xkcd:blue",
        "-"
    ),
    (
        "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/Translation_NoWeight_PCST_D8_P_NoMatching_SingleFluid/monitor/simulation",
        "No Weight, No Matching",
        "xkcd:charcoal",
        "-"
    ),
    (
        "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/Translation_DistanceWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
        "Distance Weight, Matching",
        "xkcd:orange",
        "--"
    ),
    (
        "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/Translation_DistanceSqWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
        "Distance Squared Weight, Matching",
        "xkcd:blue",
        "--"
    ),
    (
        "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/Translation_NoWeight_PCST_D8_P_Matching_SingleFluid/monitor/simulation",
        "No Weight, Matching",
        "xkcd:charcoal",
        "--"
    ),
    # ============================================================
    # D8 Case - Translation Drop - CSF
    # ============================================================
    (
        "ThreeD-drop-translation/Basak_Jibben_1200_DistanceWeightNormals/Translation_DistanceSqWeight_CSF_D8_SingleFluid/monitor/simulation",
        "CSF, D8",
        "xkcd:red",
        "-"
    ),
]


# ============================================================
# Load Simulation Data
# ============================================================

sim_data = []

def safe_float(x):
    try:
        return float(x)
    except ValueError:
        return np.nan


for filename, label, color, linestyle in sim_cases:

    data = np.loadtxt(
        filename,
        skiprows=2,
        converters=safe_float,
        ndmin=2
    )

    # Columns
    time = data[:, 1]
    Vrms = data[:, 3]

    # Remove invalid data
    valid = (
        np.isfinite(time)
        & np.isfinite(Vrms)
        & (Vrms > 0.0)
    )

    time = time[valid]
    Vrms = Vrms[valid]

    # Normalize
    time_norm = time / Tmu
    Vrms_norm = Vrms / Usigma


    # Apply optional smoothing
    Vrms_norm = smooth_data(Vrms_norm)

    sim_data.append({
        "label": label,
        "time": time_norm,
        "Vrms": Vrms_norm,
        "color": color,
        "linestyle": linestyle
    })
# ============================================================
# Basak Reference Data
# ============================================================

basak_rms = np.loadtxt(
    "ThreeD-drop-translation/Basak_Jibben_1200/BasakRMSVelocity_D8.csv",
    delimiter=",",
    skiprows=1
)

basak_max = np.loadtxt(
    "ThreeD-drop-translation/Basak_Jibben_1200/BasakMaxVelocity_D8.csv",
    delimiter=",",
    skiprows=1
)


# ============================================================
# RMS Velocity
# ============================================================

plt.figure(figsize=(12, 9))

for case in sim_data:
    plt.semilogy(
        case["time"],
        case["Vrms"],
        label=case["label"],
        linewidth=5,
        color=case["color"],
        linestyle=case["linestyle"]
    )
# Basak reference RMS velocity
# plt.semilogy(
#     basak_rms[:, 0],
#     basak_rms[:, 1],
#     label="Basak Reference",
#     linewidth=4,
#     color="xkcd:dark grey",
#     linestyle=":"
# )
plt.title("Laplace Eq")

plt.xlabel(r"$t/T_\mu$")
plt.ylabel(r"$V_{\mathrm{rms}}/U_\sigma$")

plt.minorticks_on()

plt.grid(
    which='major',
    linestyle='-',
    linewidth=0.8,
    alpha=0.7
)

plt.grid(
    which='minor',
    linestyle=':',
    linewidth=0.5,
    alpha=0.4
)
plt.xlim(0,1)
# Get existing legend entries
handles, labels = plt.gca().get_legend_handles_labels()

# Add boundary condition legend entries
# handles += [
#     Line2D([0], [0], color='black', linewidth=3, linestyle='-',
#            label='Slip BC'),
#     Line2D([0], [0], color='black', linewidth=3, linestyle='--',
#            label='Periodic BC')
# ]

plt.legend(
    handles=handles,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    handlelength=5
)

plt.tight_layout()
plt.show()

# ============================================================
# Maximum Velocity
# ============================================================

plt.figure(figsize=(12, 9))

for filename, label, color, linestyle in sim_cases:

    data = np.loadtxt(
        filename,
        skiprows=2,
        converters=safe_float,
        ndmin=2
    )

    # Columns
    time = data[:, 1]
    Vmax = data[:, 4]

    # Remove invalid data
    valid = (
        np.isfinite(time)
        & np.isfinite(Vmax)
        & (Vmax > 0.0)
    )

    time = time[valid]
    Vmax = Vmax[valid]

    # Normalize
    time_norm = time / Tmu
    Vmax_norm = Vmax / Usigma

    # Apply optional smoothing
    Vmax_norm = smooth_data(Vmax_norm)

    plt.semilogy(
        time_norm,
        Vmax_norm,
        label=label,
        linewidth=5,
        color=color,
        linestyle=linestyle
    )
# Basak reference maximum velocity
# plt.semilogy(
#     basak_max[:, 0],
#     basak_max[:, 1],
#     label="Basak Reference",
#     linewidth=4,
#     color="xkcd:dark grey",
#     linestyle=":"
# )
plt.title("Laplace Eq")

plt.xlabel(r"$t/T_\mu$")
plt.ylabel(r"$V_{\mathrm{max}}/U_\sigma$")

plt.minorticks_on()

plt.grid(
    which='major',
    linestyle='-',
    linewidth=0.8,
    alpha=0.7
)

plt.grid(
    which='minor',
    linestyle=':',
    linewidth=0.5,
    alpha=0.4
)

plt.xlim(0, 1)

# Get existing legend entries
handles, labels = plt.gca().get_legend_handles_labels()

# Add boundary condition legend entries
# handles += [
#     Line2D([0], [0], color='black', linewidth=3, linestyle='-',
#            label='Slip BC'),
#     Line2D([0], [0], color='black', linewidth=3, linestyle='--',
#            label='Periodic BC')
# ]

plt.legend(
    handles=handles,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    handlelength=5
)

plt.tight_layout()
plt.show()