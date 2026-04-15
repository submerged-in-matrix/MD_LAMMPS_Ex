import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- Read fix ave/time output files ---
def read_ave_file(filename):
    """Read LAMMPS fix ave/time output, skip comment lines."""
    data = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    vals = [float(x) for x in parts]
                    data.append(vals)
                except ValueError:
                    continue
    return np.array(data)

# --- Collect data from all temperatures ---
temperatures = [100, 200, 300, 400, 500, 600, 700, 800, 900]
ncells = 6  # 6x6x6 unit cells

results = {'T': [], 'a_mean': [], 'a_err': [], 'density_mean': [], 'density_err': [], 'pe_mean': [], 'pe_err': []}

for T in temperatures:
    fname = f'../outputs/ave_{T}K.dat'
    data = read_ave_file(fname)

    # Columns: timestep, lx, vol, pe_per_atom, density, press
    lx = data[:, 1]
    dens = data[:, 4]
    pe = data[:, 3]

    a = lx / ncells  # lattice parameter

    # Block averaging (5 blocks)
    nblocks = 5
    block_size = len(a) // nblocks
    a_blocks = [a[i*block_size:(i+1)*block_size].mean() for i in range(nblocks)]
    d_blocks = [dens[i*block_size:(i+1)*block_size].mean() for i in range(nblocks)]
    pe_blocks = [pe[i*block_size:(i+1)*block_size].mean() for i in range(nblocks)]

    results['T'].append(T)
    results['a_mean'].append(np.mean(a_blocks))
    results['a_err'].append(np.std(a_blocks) / np.sqrt(nblocks))
    results['density_mean'].append(np.mean(d_blocks))
    results['density_err'].append(np.std(d_blocks) / np.sqrt(nblocks))
    results['pe_mean'].append(np.mean(pe_blocks))
    results['pe_err'].append(np.std(pe_blocks) / np.sqrt(nblocks))

T = np.array(results['T'])
a = np.array(results['a_mean'])
a_err = np.array(results['a_err'])
dens = np.array(results['density_mean'])
pe = np.array(results['pe_mean'])

# --- Compute thermal expansion coefficient ---
coeffs = np.polyfit(T, a, 1)
a0 = coeffs[1]
slope = coeffs[0]
alpha = slope / a0

print("=" * 60)
print("THERMAL EXPANSION OF FCC COPPER (EAM - Cu_u3)")
print("=" * 60)
print(f"{'T (K)':>8} {'a (A)':>10} {'err':>10} {'density':>10} {'PE/atom':>10}")
print("-" * 60)
for i in range(len(T)):
    print(f"{T[i]:8.0f} {a[i]:10.4f} {a_err[i]:10.6f} {dens[i]:10.4f} {pe[i]:10.4f}")
print("-" * 60)
print(f"Linear fit: a(T) = {a0:.4f} + {slope:.6f} * T")
print(f"Thermal expansion coeff: alpha = {alpha*1e6:.1f} x 10^-6 K^-1")
print(f"Experimental value:      alpha = 16.5 x 10^-6 K^-1")
print(f"Deviation from experiment: {abs(alpha*1e6 - 16.5)/16.5*100:.1f}%")
print("=" * 60)

# --- Figure 1: Lattice parameter vs Temperature ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

ax = axes[0]
ax.errorbar(T, a, yerr=a_err, fmt='ro', markersize=8, capsize=4, label='MD (EAM Cu_u3)')
T_fit = np.linspace(0, 1000, 200)
a_fit = np.polyval(coeffs, T_fit)
ax.plot(T_fit, a_fit, 'b--', linewidth=1.5,
        label=f'Linear fit ($\\alpha$ = {alpha*1e6:.1f} $\\times$ 10$^{{-6}}$ K$^{{-1}}$)')
ax.set_xlabel('Temperature (K)', fontsize=14)
ax.set_ylabel('Lattice parameter ($\\AA$)', fontsize=14)
ax.set_title('(a) Thermal Expansion', fontsize=14)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# --- Figure 2: Density vs Temperature ---
ax = axes[1]
ax.plot(T, dens, 'gs-', markersize=8, label='MD (EAM Cu_u3)')
ax.set_xlabel('Temperature (K)', fontsize=14)
ax.set_ylabel('Density (g/cm$^3$)', fontsize=14)
ax.set_title('(b) Density', fontsize=14)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

# --- Figure 3: PE/atom vs Temperature ---
ax = axes[2]
ax.plot(T, pe, 'b^-', markersize=8, label='MD (EAM Cu_u3)')
ax.set_xlabel('Temperature (K)', fontsize=14)
ax.set_ylabel('PE per atom (eV)', fontsize=14)
ax.set_title('(c) Cohesive Energy', fontsize=14)
ax.legend(fontsize=11)
ax.tick_params(labelsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('../figures/thermal_expansion_copper.png', dpi=300, bbox_inches='tight')
print("\nFigure saved: ../figures/thermal_expansion_copper.png")
