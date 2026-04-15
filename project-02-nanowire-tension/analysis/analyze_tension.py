import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Read data
data = []
with open('../outputs/stress_strain.dat') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                data.append([float(parts[1]), float(parts[2])])
            except ValueError:
                continue

data = np.array(data)
strain = data[:, 0]
stress = data[:, 1]

# Smooth with window=30 (R2 was already 0.994)
w = 50
stress_s = np.convolve(stress, np.ones(w)/w, mode='valid')
strain_s = strain[(w-1)//2 : (w-1)//2 + len(stress_s)]

# --- Scan different strain regions for E ---
print(f"{'Region':>16} {'E (GPa)':>10} {'R-squared':>10}")
print("-" * 40)

best_E = 0
best_region = ""
best_r2 = 0
best_coeffs = None

regions = [
    (0.005, 0.02,  "0.5-2%"),
    (0.01,  0.03,  "1-3%"),
    (0.02,  0.04,  "2-4%"),
    (0.03,  0.05,  "3-5%"),
    (0.04,  0.06,  "4-6%"),
    (0.05,  0.07,  "5-7%"),
    (0.06,  0.08,  "6-8%"),
    (0.03,  0.07,  "3-7%"),
    (0.04,  0.08,  "4-8%"),
    (0.05,  0.08,  "5-8%"),
]

for lo, hi, name in regions:
    mask = (strain_s > lo) & (strain_s < hi)
    if mask.sum() < 5:
        continue
    c = np.polyfit(strain_s[mask], stress_s[mask], 1)
    E = c[0]
    pred = np.polyval(c, strain_s[mask])
    ss_res = np.sum((stress_s[mask] - pred) ** 2)
    ss_tot = np.sum((stress_s[mask] - stress_s[mask].mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    print(f"{name:>16} {E:10.1f} {r2:10.4f}")

    # Best = highest E with R2 > 0.99 (truly linear region)
    if r2 > 0.99 and E > best_E:
        best_E = E
        best_region = name
        best_r2 = r2
        best_coeffs = c

# Fallback if nothing had R2 > 0.99
if best_coeffs is None:
    best_idx = -1
    for lo, hi, name in regions:
        mask = (strain_s > lo) & (strain_s < hi)
        if mask.sum() < 5:
            continue
        c = np.polyfit(strain_s[mask], stress_s[mask], 1)
        E = c[0]
        pred = np.polyval(c, strain_s[mask])
        ss_res = np.sum((stress_s[mask] - pred) ** 2)
        ss_tot = np.sum((stress_s[mask] - stress_s[mask].mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        if r2 > best_r2:
            best_r2 = r2
            best_E = E
            best_region = name
            best_coeffs = c

print(f"\nBest elastic region: {best_region}")
print(f"Young's modulus: {best_E:.1f} GPa (R2={best_r2:.4f})")

# Yield point
yield_idx = np.argmax(stress_s)
yield_strain = strain_s[yield_idx]
yield_stress = stress_s[yield_idx]

print(f"\n{'='*55}")
print(f"TENSILE TEST: COPPER NANOWIRE (EAM)")
print(f"{'='*55}")
print(f"Young's modulus:    {best_E:.1f} GPa  (fit region: {best_region})")
print(f"Cu [100] expt:      ~67 GPa")
print(f"Yield stress:       {yield_stress:.1f} GPa at {yield_strain*100:.1f}%")
print(f"Cu bulk yield:      ~0.07 GPa")
print(f"{'='*55}")

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

ax = axes[0]
ax.plot(strain * 100, stress, 'b-', linewidth=0.6, alpha=0.3, label='Raw')
ax.plot(strain_s * 100, stress_s, 'b-', linewidth=2, label='Smoothed')

# Parse best region bounds for plotting
for lo, hi, name in regions:
    if name == best_region:
        fit_x = np.linspace(lo, hi, 100)
        fit_y = np.polyval(best_coeffs, fit_x)
        # Extend the fit line further to show the slope clearly
        ext_x = np.linspace(0, hi + 0.02, 100)
        ext_y = np.polyval(best_coeffs, ext_x)
        ax.plot(ext_x * 100, ext_y, 'r--', linewidth=2.5,
                label=f'E = {best_E:.0f} GPa ({best_region})')
        break

ax.plot(yield_strain * 100, yield_stress, 'rv', markersize=14,
        label=f'Yield: {yield_stress:.1f} GPa at {yield_strain*100:.1f}%')

ax.set_xlabel('Engineering Strain (%)', fontsize=14)
ax.set_ylabel('Stress (GPa)', fontsize=14)
ax.set_title('(a) Full Stress-Strain Curve', fontsize=14)
ax.legend(fontsize=10, loc='lower right')
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# Right: all slopes overlaid
ax = axes[1]
colors = plt.cm.viridis(np.linspace(0, 1, len(regions)))
for i, (lo, hi, name) in enumerate(regions):
    mask = (strain_s > lo) & (strain_s < hi)
    if mask.sum() < 5:
        continue
    c = np.polyfit(strain_s[mask], stress_s[mask], 1)
    ax.plot([lo*100, hi*100], [np.polyval(c, lo), np.polyval(c, hi)],
            'o-', color=colors[i], linewidth=2, markersize=6,
            label=f'{name}: E={c[0]:.0f} GPa')

ax.set_xlabel('Engineering Strain (%)', fontsize=14)
ax.set_ylabel('Stress (GPa)', fontsize=14)
ax.set_title('(b) E from Different Strain Regions', fontsize=14)
ax.legend(fontsize=9, loc='upper left')
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/stress_strain_curve.png', dpi=300, bbox_inches='tight')
print(f"\nFigure saved")