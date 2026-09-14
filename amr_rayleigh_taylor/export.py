import numpy as np

# =========================
# User inputs
# =========================
folder = "LargeDomain_results"
baseDirectory = f"amr_marangoni/{folder}/amr_Level5P_PPIC/monitor/"
input_file = f"{baseDirectory}statistics"
output_file = f"amr_marangoni/{folder}/amr_Level5NoP/NoP/output_P_PPIC.txt"
# Normalization
rho1 = .2
mu1 = 0.1
k1 = 0.001
sigma0 = 0.1
sigmaT = -0.1
a = 0.5
dT = 2/15

kr = 1
mur = 1
cp = 0.1

U0 = -sigmaT*a*dT/mu1
tNorm = a/U0

Vygb = -2*sigmaT*dT*a/(6*mu1+9*mu1)

muL = mu1 
muG = mu1*kr

kL = k1
kG = k1*kr

gradT = dT*a
R = a 
dSigmadT = sigmaT 

VygbB = -2 * sigmaT * gradT * R/(2*muL + 3*muG) * ((kG+2*kL)/(2*kL+kG))
# Keep every Nth row
subsample_rate = 10

# Number of header rows in the original file to skip
skiprows = 2

# =========================
# Load data
# =========================

data = np.loadtxt(input_file, skiprows=skiprows)

# Extract:
# column 2 -> index 1
# column 8 -> index 7
x = data[:, 1]/tNorm
y = data[:, 7]/Vygb


# =========================
# Subsample
# =========================

x_sub = x[::subsample_rate]
y_sub = y[::subsample_rate]


# =========================
# Write output
# =========================

with open(output_file, "w") as f:
    f.write("x, y\n")

    for x_val, y_val in zip(x_sub, y_sub):
        f.write(f"{x_val:.8g}, {y_val:.8g}\n")

print(f"Wrote {len(x_sub)} points to {output_file}")