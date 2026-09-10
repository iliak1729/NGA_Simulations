import numpy as np
import matplotlib.pyplot as plt

# Colors
color_black = '#000000'
color_orange = '#E69F00'
color_skyblue = '#56B4E9'
color_teal = '#009E73'
color_yellow = '#F0E442'
color_blue = '#0072B2'
color_red = '#D55E00'
color_pink = '#CC79A7'

# Base directory
base = "amr_ppic_laplace/OLD_LapSweep_PLIC"
base2 = "amr_ppic_laplace/OLD_LapSweep_PPIC"
base3 = "amr_ppic_laplace/OLD_IntegrationTests"
# Load data

dataLIVE = np.loadtxt(
    f"amr_ppic_laplace/monitor/statistics",
    skiprows=2
)

data_1200_CSF = np.loadtxt(
    f"{base}/amr_LA1200CSF/monitor/statistics",
    skiprows=2
)

data_1200_NoP = np.loadtxt(
    f"{base}/amr_LA1200NoP/monitor/statistics",
    skiprows=2
)

data_1200_P = np.loadtxt(
    f"{base}/amr_LA1200P/monitor/statistics",
    skiprows=2
)

data_120_CSF = np.loadtxt(
    f"{base}/amr_LA120CSF/monitor/statistics",
    skiprows=2
)

data_120_NoP = np.loadtxt(
    f"{base}/amr_LA120NoP/monitor/statistics",
    skiprows=2
)

data_120_P = np.loadtxt(
    f"{base}/amr_LA120P/monitor/statistics",
    skiprows=2
)

# PPIC
data_120_CSF_PPIC = np.loadtxt(
    f"{base2}/amr_LA120CSF/monitor/statistics",
    skiprows=2
)

data_120_NoP_PPIC = np.loadtxt(
    f"{base2}/amr_LA120NoP/monitor/statistics",
    skiprows=2
)

data_120_NoP_PPIC_5Point = np.loadtxt(
    f"{base3}/amr_LA120NoP/monitor/statistics",
    skiprows=2
)


data_120_P_PPIC_5Point = np.loadtxt(
    f"{base3}/amr_LA120P/monitor/statistics",
    skiprows=2
)

data_120_P_PPIC = np.loadtxt(
    f"{base2}/amr_LA120P/monitor/statistics",
    skiprows=2
)

data_1200_CSF_PPIC = np.loadtxt(
    f"{base2}/amr_LA1200CSF/monitor/statistics",
    skiprows=2
)

data_1200_NoP_PPIC = np.loadtxt(
    f"{base2}/amr_LA1200NoP/monitor/statistics",
    skiprows=2
)

data_1200_P_PPIC = np.loadtxt(
    f"{base2}/amr_LA1200P/monitor/statistics",
    skiprows=2
)


plt.figure()

# La = 1200
plt.semilogy(
    data_1200_NoP[:, 2],
    data_1200_NoP[:, 6],
    label="PUST NoP, La1200",
    linewidth=4,
    linestyle="--",
    color=color_teal
)

plt.semilogy(
    data_1200_P[:, 2],
    data_1200_P[:, 6],
    label="PUST P, La1200",
    linewidth=4,
    linestyle="--",
    color=color_blue
)

plt.semilogy(
    data_1200_CSF[:, 2],
    data_1200_CSF[:, 6],
    label="CSF, La1200",
    linewidth=4,
    linestyle="--",
    color=color_yellow
)


# La = 120
plt.semilogy(
    data_120_NoP[:, 2],
    data_120_NoP[:, 6],
    label="PUST NoP, La120",
    linewidth=4,
    linestyle="-",
    color=color_teal
)

plt.semilogy(
    data_120_P[:, 2],
    data_120_P[:, 6],
    label="PUST P, La120",
    linewidth=4,
    linestyle="-",
    color=color_blue
)

plt.semilogy(
    data_120_CSF[:, 2],
    data_120_CSF[:, 6],
    label="CSF, La120",
    linewidth=4,
    linestyle="-",
    color=color_yellow
)
# PPIC
plt.semilogy(
    data_120_CSF_PPIC[:, 2],
    data_120_CSF_PPIC[:, 6],
    label="CSF_PPIC, La120",
    linewidth=4,
    linestyle="-",
    color=color_orange
)

plt.semilogy(
    data_1200_CSF_PPIC[:, 2],
    data_1200_CSF_PPIC[:, 6],
    label="CSF_PPIC, La1200",
    linewidth=4,
    linestyle="--",
    color=color_orange
)

plt.semilogy(
    data_120_P_PPIC[:, 2],
    data_120_P_PPIC[:, 6],
    label="PPIC_withP,La120",
    linewidth=2,
    linestyle="-",
    color=color_black
)

plt.semilogy(
    data_1200_P_PPIC[:, 2],
    data_1200_P_PPIC[:, 6],
    label="PPIC_withP,La1200",
    linewidth=2,
    linestyle="--",
    color=color_black
)

plt.semilogy(
    data_120_NoP_PPIC[:, 2],
    data_120_NoP_PPIC[:, 6],
    label="PPIC_NoP,La120",
    linewidth=2,
    linestyle="-",
    color=color_pink
)

plt.semilogy(
    data_1200_NoP_PPIC[:, 2],
    data_1200_NoP_PPIC[:, 6],
    label="PPIC_NoP,La1200",
    linewidth=2,
    linestyle="-",
    color=color_pink
)

# plt.ylim([3e-5, 1e-3])
plt.xlim([0, 1.0])
plt.title("Translating Drop")
plt.xlabel("t/Tv")
plt.ylabel("Vrms")

plt.legend()
plt.tight_layout()

plt.figure()

plt.semilogy(
    data_120_NoP_PPIC_5Point[:, 2],
    data_120_NoP_PPIC_5Point[:, 6],
    label="PPIC_NoP,La120, 5 Point",
    linewidth=2,
    linestyle="-",
    color=color_black
)

plt.semilogy(
    data_120_NoP_PPIC[:, 2],
    data_120_NoP_PPIC[:, 6],
    label="PPIC_NoP,La120, 50 Point",
    linewidth=2,
    linestyle="--",
    color=color_pink
)

plt.semilogy(
    data_120_P_PPIC_5Point[:, 2],
    data_120_P_PPIC_5Point[:, 6],
    label="PPIC_P,La120, 5 Point",
    linewidth=2,
    linestyle="-",
    color=color_red
)

plt.semilogy(
    data_120_P_PPIC[:, 2],
    data_120_P_PPIC[:, 6],
    label="PPIC_P,La120, 50 Point",
    linewidth=2,
    linestyle="--",
    color=color_blue
)
plt.ylim([1e-4,2.5e-4])
plt.xlim([0, 0.05])
plt.title("Translating Drop")
plt.xlabel("t/Tv")
plt.ylabel("Vrms")

plt.legend()
plt.tight_layout()


plt.show()