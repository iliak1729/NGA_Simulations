import numpy as np

# =========================
# User inputs
# =========================
folder = "OLD_LapSweep_PLIC"
run = "amr_LA1200CSF"
baseDirectory = f"amr_ppic_laplace/{folder}/{run}/monitor/"
outputDirectory = f"amr_ppic_laplace/subsampled_results/"
input_file = f"{baseDirectory}statistics"
output_file = f"{outputDirectory}subsampled_{run}_PLIC.txt"

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
x = data[:, 2]
y = data[:, 6]


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