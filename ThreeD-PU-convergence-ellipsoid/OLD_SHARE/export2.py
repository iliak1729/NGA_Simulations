from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


# Colors - Okabe & Ito
color_black = '#000000'
color_orange = '#E69F00'
color_skyblue = '#56B4E9'
color_teal = '#009E73'
color_yellow = '#F0E442'
color_blue = '#0072B2'
color_red = '#D55E00'
color_pink = '#CC79A7'

def safe_float(x):
    if isinstance(x, bytes):
        x = x.decode()

    x = x.strip().replace("D", "E")

    try:
        return float(x)
    except ValueError:
        return np.nan


# Get Data
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_2D_EllipsoidStress_SpreadSweepACOS"
d1 = np.loadtxt(datafile, skiprows=2,converters=(safe_float))
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_2D_EllipsoidStress_SpreadSweepACOS_Jibben"
d2 = np.loadtxt(datafile, skiprows=2,converters=(safe_float))
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_3D_EllipsoidStress_SpreadSweepACOS"
d3 = np.loadtxt(datafile, skiprows=2,converters=(safe_float))
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_3D_EllipsoidStress_SpreadSweepACOS_Jibben"
d4 = np.loadtxt(datafile, skiprows=2,converters=(safe_float))

# Extra Values for N=  16,32,64,128,256,512 for 2D and N= 16,32,64,128 for 3D
data2D16 = d1[d1[:, 0] == 16]
data2D32 = d1[d1[:, 0] == 32]
data2D64 = d1[d1[:, 0] == 64]
data2D128 = d1[d1[:, 0] == 128]
data2D256 = d1[d1[:, 0] == 256]
data2D512 = d1[d1[:, 0] == 512]

data2D16J = d2[d2[:, 0] == 16]
data2D32J = d2[d2[:, 0] == 32]
data2D64J = d2[d2[:, 0] == 64]
data2D128J = d2[d2[:, 0] == 128]
data2D256J = d2[d2[:, 0] == 256]
data2D512J = d2[d2[:, 0] == 512]

data3D16 = d3[d3[:, 0] == 16]
data3D32 = d3[d3[:, 0] == 32]
data3D64 = d3[d3[:, 0] == 64]
data3D128 = d3[d3[:, 0] == 128]

data3D16J = d4[d4[:, 0] == 16]
data3D32J = d4[d4[:, 0] == 32]
data3D64J = d4[d4[:, 0] == 64]
data3D128J = d4[d4[:, 0] == 128]

print("Data Loaded")
print("Check Data:")

print(f"2D Data Points (N=16): {len(data2D16)}")
print(f"2D Data Points (N=32): {len(data2D32)}")
print(f"2D Data Points (N=64): {len(data2D64)}")
print(f"2D Data Points (N=128): {len(data2D128)}")
print(f"2D Data Points (N=256): {len(data2D256)}")
print(f"2D Data Points (N=512): {len(data2D512)}")

print(f"2D Jibben Data Points (N=16): {len(data2D16J)}")
print(f"2D Jibben Data Points (N=32): {len(data2D32J)}")
print(f"2D Jibben Data Points (N=64): {len(data2D64J)}")
print(f"2D Jibben Data Points (N=128): {len(data2D128J)}")
print(f"2D Jibben Data Points (N=256): {len(data2D256J)}")
print(f"2D Jibben Data Points (N=512): {len(data2D512J)}")

print(f"3D Data Points (N=16): {len(data3D16)}")
print(f"3D Data Points (N=32): {len(data3D32)}")
print(f"3D Data Points (N=64): {len(data3D64)}")
print(f"3D Data Points (N=128): {len(data3D128)}")

print(f"3D Jibben Data Points (N=16): {len(data3D16J)}")
print(f"3D Jibben Data Points (N=32): {len(data3D32J)}")
print(f"3D Jibben Data Points (N=64): {len(data3D64J)}")
print(f"3D Jibben Data Points (N=128): {len(data3D128J)}")

# Print Shapes of one output
print(f"Shape of 2D Data (N=16): {data2D16.shape}")


# Columns to plot
columns = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]

names = [
    'Radius Error',
    'Tangent Error',
    'Curvature Error',
    'Sigma_xx Error',
    'Sigma_xy Error',
    'Sigma_xz Error',
    'Sigma_xy Error',
    'Sigma_yy Error',
    'Sigma_yz Error',
    'Sigma_xz Error',
    'Sigma_yz Error',
    'Sigma_zz Error',
    'Sigma Total Error',
    'X Force Error',
    'Y Force Error',
    'Z Force Error'
]

output_dir = Path("ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Output")
output_dir.mkdir(parents=True, exist_ok=True)

tikz_output_dir = output_dir / "tikz_data_b"
tikz_output_dir.mkdir(parents=True, exist_ok=True)


def write_xy(filename, x, y):
    """
    Write x-y data in a PGFPlots-friendly CSV format.
    """
    with open(filename, "w") as f:
        f.write("x, y\n")

        for x_val, y_val in zip(x, y):
            if np.isfinite(x_val) and np.isfinite(y_val):
                f.write(f"{x_val:.8g}, {y_val:.8g}\n")


count = 0

for col, name in zip(columns, names):

    count += 1

    # Safe filename prefix
    prefix = f"{count}_{name.replace(' ', '_')}"

    # =========================================================
    # 2D
    # =========================================================

    curves_2D = [
        (16,  data2D16),
        (32,  data2D32),
        (64,  data2D64),
        (128, data2D128),
        (256, data2D256),
        (512, data2D512)
    ]

    plt.figure()

    for N, data in curves_2D:

        x = data[:, 2]
        y = data[:, col]

        plt.semilogy(
            x,
            y,
            'o-',
            label=f'N={N}'
        )

        write_xy(
            tikz_output_dir / f"{prefix}_2D_N{N}.txt",
            x,
            y
        )

    plt.xlabel('Spread')
    plt.ylabel(name)
    plt.title(f'{name} vs Spread for 2D Ellipse')
    plt.legend()

    plt.savefig(
        output_dir / f"{prefix}_2D.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


    # =========================================================
    # 2D Jibben
    # =========================================================

    curves_2DJ = [
        (16,  data2D16J),
        (32,  data2D32J),
        (64,  data2D64J),
        (128, data2D128J),
        (256, data2D256J),
        (512, data2D512J)
    ]

    plt.figure()

    for N, data in curves_2DJ:

        x = data[:, 2]
        y = data[:, col]

        plt.semilogy(
            x,
            y,
            'o-',
            label=f'N={N}'
        )

        write_xy(
            tikz_output_dir / f"{prefix}_2D_Jibben_N{N}.txt",
            x,
            y
        )

    plt.xlabel('Spread')
    plt.ylabel(name)
    plt.title(f'{name} vs Spread for 2D Ellipse (Jibben)')
    plt.legend()

    plt.savefig(
        output_dir / f"{prefix}_2D_Jibben.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


    # =========================================================
    # 3D
    # =========================================================

    curves_3D = [
        (16,  data3D16),
        (32,  data3D32),
        (64,  data3D64),
        (128, data3D128)
    ]

    plt.figure()

    for N, data in curves_3D:

        x = data[:, 2]
        y = data[:, col]

        plt.semilogy(
            x,
            y,
            'o-',
            label=f'N={N}'
        )

        write_xy(
            tikz_output_dir / f"{prefix}_3D_N{N}.txt",
            x,
            y
        )

    plt.xlabel('Spread')
    plt.ylabel(name)
    plt.title(f'{name} vs Spread for 3D Ellipse')
    plt.legend()

    plt.savefig(
        output_dir / f"{prefix}_3D.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


    # =========================================================
    # 3D Jibben
    # =========================================================

    curves_3DJ = [
        (16,  data3D16J),
        (32,  data3D32J),
        (64,  data3D64J),
        (128, data3D128J)
    ]

    plt.figure()

    for N, data in curves_3DJ:

        x = data[:, 2]
        y = data[:, col]

        plt.semilogy(
            x,
            y,
            'o-',
            label=f'N={N}'
        )

        write_xy(
            tikz_output_dir / f"{prefix}_3D_Jibben_N{N}.txt",
            x,
            y
        )

    plt.xlabel('Spread')
    plt.ylabel(name)
    plt.title(f'{name} vs Spread for 3D Ellipse (Jibben)')
    plt.legend()

    plt.savefig(
        output_dir / f"{prefix}_3D_Jibben.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved plots and TikZ data for {name}")