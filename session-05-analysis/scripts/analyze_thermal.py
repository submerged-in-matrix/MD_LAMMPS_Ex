import numpy as np
import matplotlib
matplotlib.use('Agg') # No display required
import matplotlib.pyplot as plt

# data from thermodynamic behavior from session 3 for Cu
T = np.array([100, 300, 500])
lx = np.array([18.170, 18.191, 18.193])  # Box length
a = lx / 5.0

# Compute thermal expansion co-efficient (alpha = (1/a0) * da/dT)
coeffs = np.polyfit(T, a, 1) #a = a0 + slope*T
a0 = coeffs[1]
slope = coeffs[0]
alpha = slope / a0

print(f"Lattice constant: {a}")
print(f"Linear fit: a = {a0:.4f} + {slope:.6f} * T")
print(f"Thermal expansion coefficient: {alpha*1e6:.1f} x 10^-6 K^-1")
print(f"Experimental value for Cu:     16.5 x 10^-6 K^-1")

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(T, a, 'ro', markersize=10, label='MD simulation (EAM)')
T_fit = np.linspace(0, 600, 100)
a_fit = np.polyval(coeffs, T_fit)
ax.plot(T_fit, a_fit, 'b--', label=f'Linear fit ($\\alpha$ = {alpha*1e6:.1f} x 10$^{{-6}}$ K$^{{-1}}$)')
ax.set_xlabel('Temperature (K)', fontsize=14)
ax.set_ylabel('Lattice parameter (Å)', fontsize=14)
ax.set_title('Thermal Expansion of FCC Copper (EAM)', fontsize=16)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('../outputs/thermal_expansion.png', dpi=300)
print("Plot saved to ../outputs/thermal_expansion.png")
