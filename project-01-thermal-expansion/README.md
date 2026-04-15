# Thermal Expansion of FCC Copper via Molecular Dynamics

## Summary

Computed the linear thermal expansion coefficient of FCC copper
from 100–900 K using classical MD with the EAM potential (Cu_u3).
Result: α = 18.1 × 10⁻⁶ K⁻¹ (experimental: 16.5 × 10⁻⁶ K⁻¹,
deviation: 9.7%).

## Method

- **Software:** LAMMPS (stable release, serial)
- **Potential:** EAM (Cu_u3.eam, Foiles/Baskes/Daw parameterization)
- **System:** 864 atoms (6×6×6 FCC unit cells), periodic boundaries
- **Protocol per temperature:**
  1. Energy minimization (CG, etol=1e-6)
  2. NPT equilibration at target T and P=0 bar (20,000 steps, dt=1 fs)
  3. NPT production (30,000 steps), lx averaged via fix ave/time
- **Temperatures:** 100, 200, 300, 400, 500, 600, 700, 800, 900 K
- **Each temperature run is independent** (starts from perfect lattice)

## Key Result

| T (K) | a (Å)  | Density (g/cm³) | PE/atom (eV) |
|--------|--------|-----------------|--------------|
| 100    | 3.6223 | 8.881           | -3.527       |
| 200    | 3.6266 | 8.849           | -3.513       |
| 300    | 3.6349 | 8.789           | -3.501       |
| 400    | 3.6369 | 8.774           | -3.488       |
| 500    | 3.6424 | 8.735           | -3.473       |
| 600    | 3.6596 | 8.613           | -3.458       |
| 700    | 3.6612 | 8.601           | -3.443       |
| 800    | 3.6714 | 8.529           | -3.428       |
| 900    | 3.6751 | 8.504           | -3.404       |

Linear fit: a(T) = 3.6134 + 6.5 × 10⁻⁵ T

**α = 18.1 × 10⁻⁶ K⁻¹** (experiment: 16.5 × 10⁻⁶ K⁻¹)

## Files

- `scripts/in.cu_thermal` — Parameterized LAMMPS input (pass -var T)
- `scripts/run_all.sh` — Runs all 9 temperatures
- `analysis/analyze_expansion.py` — Extracts α with error bars
- `figures/thermal_expansion_copper.png` — Publication figure
- `outputs/` — Raw LAMMPS logs and fix ave/time data

## How to Reproduce

```bash
cd scripts
cp /usr/share/lammps/potentials/Cu_u3.eam .
./run_all.sh
cd ../analysis
python3 analyze_expansion.py
```

## Discussion

The 9.7% overestimate of α is consistent with known EAM limitations:
the Cu_u3 potential was fitted primarily to 0 K elastic and structural
properties. Anharmonic effects at high temperature are captured only
approximately. The 600 K data point shows a slight jump, possibly
indicating the onset of enhanced anharmonicity or insufficient
equilibration at that specific temperature. Longer production runs
(100,000+ steps) and a quadratic fit would improve the high-T
accuracy.

## References

- Thompson et al. (2022) Comput. Phys. Commun. 271, 108171 (LAMMPS)
- Foiles, Baskes, Daw (1986) Phys. Rev. B 33, 7983 (EAM potential)
- Touloukian et al. (1975) Thermophysical Properties of Matter (exp. α)
