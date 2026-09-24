import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Radial coordinate
# r = distance / support radius
# ============================================================

r = np.linspace(0.0, 1.0, 10000)


# ============================================================
# 1. Wendland C2
# ============================================================

def wendland_c2(r):
    return (1.0 - r)**4 * (4.0*r + 1.0)


# ============================================================
# 2. Wu-type compact kernel
#
# Compact polynomial kernel used here to represent a
# qualitatively different CSRBF profile.
# ============================================================

def wu(r):
    return (1.0 - r)**5 * (5.0*r**4 + 25.0*r**3
                            + 48.0*r**2 + 40.0*r + 8.0) / 8.0


# ============================================================
# 3. Buhmann-type compact kernel
#
# Example compact non-polynomial profile.
# ============================================================

def buhmann(r):
    # Avoid numerical issues exactly at r = 0
    rr = np.maximum(r, 1.0e-14)

    phi = (
        2.0*rr**4*np.log(rr)
        - (7.0/2.0)*rr**4
        + (16.0/3.0)*rr**3
        - 2.0*rr**2
        + 1.0/6.0
    )

    # Normalize to phi(0) = 1
    return 6.0*phi


# ============================================================
# 4. Gneiting-type compact kernel
#
# A simple member with a distinctly different compact profile.
# ============================================================

def gneiting(r):
    return (1.0 - r)**4 * (1.0 + 4.0*r)


# ============================================================
# 5. Cubic smoothstep
# ============================================================

def cubic_smoothstep(r):
    return 1.0 - 3.0*r**2 + 2.0*r**3


# ============================================================
# 6. Quintic smootherstep
# ============================================================

def quintic_smootherstep(r):
    return 1.0 - 10.0*r**3 + 15.0*r**4 - 6.0*r**5


# ============================================================
# 7. Quartic / biweight
# ============================================================

def biweight(r):
    return (1.0 - r**2)**2


# ============================================================
# 8. C-infinity compact bump
# ============================================================

def compact_bump(r):

    phi = np.zeros_like(r)

    mask = r < 1.0

    phi[mask] = np.exp(
        -r[mask]**2 / (1.0 - r[mask]**2)
    )

    return phi


# ============================================================
# 9. Modified Shepard
#
# Singular at r = 0. epsilon regularizes the function for
# plotting purposes.
# ============================================================

def modified_shepard(r, p=2, eps=0.05):

    rr = np.maximum(r, eps)

    phi = ((1.0 - rr) / rr)**p

    # Normalize by value at r = eps so maximum plotted value = 1
    phi /= ((1.0 - eps) / eps)**p

    return phi


# ============================================================
# 10. Cosine / Hann
# ============================================================

def cosine_weight(r):
    return 0.5 * (1.0 + np.cos(np.pi*r))


# ============================================================
# Evaluate weighting functions
# ============================================================

weights = {
    r"Wendland $C^2$":       wendland_c2(r),
    "Wu":                    wu(r),
    "Buhmann":               buhmann(r),
    "Gneiting":              gneiting(r),
    "Cubic Smoothstep":      cubic_smoothstep(r),
    "Quintic Smootherstep":  quintic_smootherstep(r),
    "Quartic / Biweight":    biweight(r),
    r"$C^\infty$ Bump":      compact_bump(r),
    "Modified Shepard":      modified_shepard(r),
    "Cosine / Hann":         cosine_weight(r),
}


# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(9, 6))

for label, phi in weights.items():
    plt.plot(
        r,
        phi,
        linewidth=2,
        label=label
    )

plt.xlabel(r"$r/R$")
plt.ylabel(r"$\psi(r)$")

plt.xlim(0.0, 1.0)
plt.ylim(0.0, 1.05)

plt.grid(True, alpha=0.25)

plt.legend(
    loc="best",
    fontsize=9,
    frameon=True
)

plt.tight_layout()
plt.show()