from turtle import color

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import cm
import matplotlib.collections as mcoll
import matplotlib.path as mpath
import matplotlib.animation as animation
import math
from scipy import special
from pathlib import Path

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
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_2D_EllipsoidStress_SpreadSweep"
d1 = np.loadtxt(datafile, skiprows=2,converters=(safe_float))
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_2D_EllipsoidStress_SpreadSweep_Jibben"
d2 = np.loadtxt(datafile, skiprows=2,converters=(safe_float))
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_3D_EllipsoidStress_SpreadSweep"
d3 = np.loadtxt(datafile, skiprows=2,converters=(safe_float))
datafile = "ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Final/Error_Values_3D_EllipsoidStress_SpreadSweep_Jibben"
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

# Now we take the radius error trends in column 7 and plot with respect to spread in column 3
columns = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
names = ['Radius Error','Tangent Error','Curvature Error',
         'Sigma_xx Error', 'Sigma_xy Error','Sigma_xz Error',
         'Sigma_xy Error', 'Sigma_yy Error','Sigma_yz Error',
         'Sigma_xz Error', 'Sigma_yz Error','Sigma_zz Error',
         'Sigma Total Error','X Force Error','Y Force Error','Z Force Error']
output_dir = Path("ThreeD-PU-convergence-ellipsoid/OLD_SHARE/Output")
output_dir.mkdir(parents=True, exist_ok=True)
count = 0

# for col, name in zip(columns, names):
#     count = count + 1
#     # 2D 
#     plt.figure()
#     plt.semilogy(data2D16[:, 2], data2D16[:, col], 'o-', label='N=16', color=color_orange)
#     plt.semilogy(data2D32[:, 2], data2D32[:, col], 'o-', label='N=32', color=color_skyblue)
#     plt.semilogy(data2D64[:, 2], data2D64[:, col], 'o-', label='N=64', color=color_teal)
#     plt.semilogy(data2D128[:, 2], data2D128[:, col], 'o-', label='N=128', color=color_yellow)
#     plt.semilogy(data2D256[:, 2], data2D256[:, col], 'o-', label='N=256', color=color_blue)
#     plt.semilogy(data2D512[:, 2], data2D512[:, col], 'o-', label='N=512', color=color_red)

#     plt.xlabel('Spread')
#     plt.ylabel(name)
#     plt.title(f'{name} vs Spread for 2D Ellipse')
#     plt.legend()
#     plt.savefig(output_dir / f"{count}_{name.replace(' ', '_')}_2D.png", dpi=300, bbox_inches="tight")
#     # 2D Jibben
#     plt.figure()
#     plt.semilogy(data2D16J[:, 2], data2D16J[:, col], 'o-', label='N=16', color=color_orange)
#     plt.semilogy(data2D32J[:, 2], data2D32J[:, col], 'o-', label='N=32', color=color_skyblue)
#     plt.semilogy(data2D64J[:, 2], data2D64J[:, col], 'o-', label='N=64', color=color_teal)
#     plt.semilogy(data2D128J[:, 2], data2D128J[:, col], 'o-', label='N=128', color=color_yellow)
#     plt.semilogy(data2D256J[:, 2], data2D256J[:, col], 'o-', label='N=256', color=color_blue)
#     plt.semilogy(data2D512J[:, 2], data2D512J[:, col], 'o-', label='N=512', color=color_red)

#     plt.xlabel('Spread')
#     plt.ylabel(name)
#     plt.title(f'{name} vs Spread for 2D Ellipse (Jibben)')
#     plt.legend()
#     plt.savefig(output_dir / f"{count}_{name.replace(' ', '_')}_2D_Jibben.png", dpi=300, bbox_inches="tight")
#     # 3D 
#     plt.figure()
#     plt.semilogy(data3D16[:, 2], data3D16[:, col], 'o-', label='N=16', color=color_orange)
#     plt.semilogy(data3D32[:, 2], data3D32[:, col], 'o-', label='N=32', color=color_skyblue)
#     plt.semilogy(data3D64[:, 2], data3D64[:, col], 'o-', label='N=64', color=color_teal)
#     plt.semilogy(data3D128[:, 2], data3D128[:, col], 'o-', label='N=128', color=color_yellow)

#     plt.xlabel('Spread')
#     plt.ylabel(name)
#     plt.title(f'{name} vs Spread for 3D Ellipse')
#     plt.legend()
#     plt.savefig(output_dir / f"{count}_{name.replace(' ', '_')}_3D.png", dpi=300, bbox_inches="tight")
#     # 3D Jibben
#     plt.figure()
#     plt.semilogy(data3D16J[:, 2], data3D16J[:, col], 'o-', label='N=16', color=color_orange)
#     plt.semilogy(data3D32J[:, 2], data3D32J[:, col], 'o-', label='N=32', color=color_skyblue)
#     plt.semilogy(data3D64J[:, 2], data3D64J[:, col], 'o-', label='N=64', color=color_teal)
#     plt.semilogy(data3D128J[:, 2], data3D128J[:, col], 'o-', label='N=128', color=color_yellow)

#     plt.xlabel('Spread')
#     plt.ylabel(name)
#     plt.title(f'{name} vs Spread for 3D Ellipse (Jibben)')
#     plt.legend()
#     plt.savefig(output_dir / f"{count}_{name.replace(' ', '_')}_3D_Jibben.png", dpi=300, bbox_inches="tight")
    
#     plt.close('all')  # Close all figures to free memory
#     print(f"Saved plots for {name}")



# Grid sizes
N2D = np.array([16, 32, 64, 128, 256, 512])
N3D = np.array([16, 32, 64, 128])

datasets2D = [data2D16, data2D32, data2D64, data2D128, data2D256, data2D512]
datasets2DJ = [data2D16J, data2D32J, data2D64J, data2D128J, data2D256J, data2D512J]
datasets3D = [data3D16, data3D32, data3D64, data3D128]
datasets3DJ = [data3D16J, data3D32J, data3D64J, data3D128J]

# Columns to plot
cols = [6, 7, 8,9]
names = ["Radius Error", "Tangent Error", "Curvature Error","Sigma_xx Error"]

for col, name in zip(cols, names):

    def extract_values(dataset_list):
        values = []
        for data in dataset_list:
            row = data[np.isclose(data[:,2], 2.0)]
            if len(row) == 0:
                values.append(np.nan)
            else:
                values.append(row[0, col])
        return np.array(values)

    y2D  = extract_values(datasets2D)
    y2DJ = extract_values(datasets2DJ)
    y3D  = extract_values(datasets3D)
    y3DJ = extract_values(datasets3DJ)

    plt.figure(figsize=(6,5))

    plt.loglog(N2D, y2D,  'o-', label='2D', color=color_blue)
    plt.loglog(N2D, y2DJ, 's-', label='2D Jibben', color=color_red)
    plt.loglog(N3D, y3D,  '^-', label='3D', color=color_teal)
    plt.loglog(N3D, y3DJ, 'd-', label='3D Jibben', color=color_orange)



    anchor_N = {
    1: N2D[0],     # First-order line passes through N=64
    2: N2D[0],    # Second-order line passes through N=128
    3: N2D[0],     # Third-order line passes through N=256
    4: N2D[0],     
    5: N2D[0]     # Fourth-order line passes through N=512
    }

    anchor_error = {
        1: y2D[0],
        2: y2D[0],
        3: y2D[0],
        4: y2DJ[0],
        5: y2DJ[0]
    }

    Nref = np.array([16, 32, 64, 128, 256, 512])

    styles = {1: '--', 2: ':', 3: '-.',4: '--',5: '-'}
    orders = [1, 2, 3,4,5]
    if col >= 8:
        orders = [1, 2]
        anchor_error = {1: y2DJ[0], 2: y2DJ[0]}
    for p in orders:
        ref = anchor_error[p] * (anchor_N[p] / Nref)**p
        plt.loglog(
            Nref,
            ref,
            styles[p],
            color='k',
            linewidth=1.5,
            label=f'{p}$^{{{"st" if p==1 else "nd" if p==2 else "rd" if p==3 else "th" }}}$ Order'
    )


    plt.xlabel("Grid Resolution ($N$)")
    plt.ylabel(name)
    plt.title(f"{name} at Spread = 2.0")
    plt.grid(True, which='both', ls=':')
    plt.legend()
    plt.show()
    # plt.savefig(output_dir / f"{name.replace(' ','_')}_Spread2_Convergence.png",
    #             dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"Saved {name}")