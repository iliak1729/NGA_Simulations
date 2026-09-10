import numpy as np

# =========================
# User inputs
# =========================
baseDirectory = "amr_marangoni/OLD_LargeDomain/amr_Level5P/monitor/"
input_file = f"{baseDirectory}statistics"
output_file = f"{baseDirectory}subsampled.txt"

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
x = data[:, 1]
y = data[:, 7]


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