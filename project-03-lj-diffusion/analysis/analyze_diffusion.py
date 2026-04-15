import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Read MSD data
data = []
with open('../outputs/msd_T1.0_rho0.8.dat') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                data.append([float(parts[0]), float(parts[1])])
            except ValueError:
                continue

data = np.array(data)
timestep_num = data[:, 0]
msd = data[:, 1]

# Convert timestep number to time in LJ units
dt = 0.005
time = timestep_num * dt

# --- Identify ballistic and diffusive regimes ---
# Ballistic: t < 1 tau (MSD ~ t^2)
# Diffusive: t > 50 tau (MSD ~ t, slope = 6D)

# Fit diffusive regime (t > 50 tau, skip early transient)
diff_mask = time > 50.0
coeffs = np.polyfit(time[diff_mask], msd[diff_mask], 1)
slope = coeffs[0]
D = slope / 6.0  # Einstein relation: MSD = 6*D*t in 3D

# Published reference value
D_ref = 0.06  # Meier et al., J. Chem. Phys. (2004)

print("=" * 55)
print("SELF-DIFFUSION: LJ FLUID")
print("=" * 55)
print(f"State point: T* = 1.0, rho* = 0.8")
print(f"System: 2048 atoms, production: 1000 tau")
print(f"Ensemble: NVE (after NVT equilibration)")
print("-" * 55)
print(f"MSD slope (d(MSD)/dt):  {slope:.4f} sigma^2/tau")
print(f"D* = slope/6:           {D:.4f}")
print(f"Published D* (ref):     {D_ref}")
print(f"Deviation:              {abs(D - D_ref)/D_ref*100:.1f}%")
print("=" * 55)

# --- Plot ---
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# Left: MSD vs time (linear scale)
ax = axes[0]
ax.plot(time, msd, 'b-', linewidth=1.5, label='MD data')
fit_t = np.linspace(50, time[-1], 100)
fit_msd = np.polyval(coeffs, fit_t)
ax.plot(fit_t, fit_msd, 'r--', linewidth=2.5,
        label=f'Linear fit (t > 50): D* = {D:.4f}')
ax.axvline(x=50, color='gray', linestyle=':', alpha=0.5, label='Fit start')
ax.set_xlabel('Time ($\\tau$)', fontsize=14)
ax.set_ylabel('MSD ($\\sigma^2$)', fontsize=14)
ax.set_title('(a) Mean-Squared Displacement', fontsize=14)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# Right: log-log plot (shows ballistic -> diffusive crossover)
ax = axes[1]
# Skip first point if MSD is zero/tiny
pos_mask = (time > 0.5) & (msd > 0.1)
ax.loglog(time[pos_mask], msd[pos_mask], 'b-', linewidth=1.5, label='MD data')

# Reference slopes
t_ref = np.logspace(-1, 1, 50)
ax.loglog(t_ref, 3 * t_ref**2, 'g--', linewidth=1.5, alpha=0.6,
          label='Ballistic: MSD ~ $t^2$')
t_ref2 = np.logspace(1, 3, 50)
ax.loglog(t_ref2, 6 * D * t_ref2, 'r--', linewidth=1.5, alpha=0.6,
          label=f'Diffusive: MSD = 6D*t')

ax.set_xlabel('Time ($\\tau$)', fontsize=14)
ax.set_ylabel('MSD ($\\sigma^2$)', fontsize=14)
ax.set_title('(b) Log-Log: Ballistic to Diffusive Crossover', fontsize=14)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
os.makedirs('../figures', exist_ok=True)
plt.savefig('../figures/diffusion_lj.png', dpi=300, bbox_inches='tight')
print(f"\nFigure saved")
