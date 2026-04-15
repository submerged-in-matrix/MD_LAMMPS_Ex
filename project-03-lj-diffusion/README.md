# Self-Diffusion Coefficient of a Lennard-Jones Fluid

## Summary

Computed the self-diffusion coefficient of a LJ fluid at T* = 1.0,
rho* = 0.8 using the Einstein relation (MSD slope). Result:
D* = 0.064 (published reference: D* = 0.06, deviation: 7%).

## Method

- **Software:** LAMMPS (serial)
- **Potential:** LJ/cut, epsilon=1.0, sigma=1.0, cutoff=2.5
- **System:** 2048 atoms, periodic boundaries
- **Protocol:**
  1. Energy minimization
  2. NVT equilibration at T*=1.0 (50,000 steps, melts FCC into liquid)
  3. NVE production (200,000 steps = 1000 tau)
- **Critical decision:** NVE for production, not NVT. Thermostat
  corrupts velocity autocorrelation and MSD, giving wrong D.

## Key Result

| Property | MD Result | Published Reference |
|----------|-----------|-------------------|
| D* (reduced units) | 0.064 | 0.06 (Meier 2004) |
| Deviation | 7% | — |

## Physics

MSD grows linearly with time in the diffusive regime:
MSD(t) = 6*D*t (3D). The slope of the linear fit for t > 50 tau
gives D. At short times (t < 1 tau), motion is ballistic (MSD ~ t^2)
because atoms haven't collided yet. The crossover from ballistic
to diffusive is the emergence of macroscopic transport from
microscopic collisions.

## Files

- `scripts/in.lj_diffusion` — Parameterized input (-var TEMP, -var DENS)
- `analysis/analyze_diffusion.py` — MSD analysis + Einstein relation
- `figures/diffusion_lj.png` — MSD + log-log crossover figure
- `outputs/msd_T1.0_rho0.8.dat` — Raw MSD data

## How to Reproduce

```bash
cd scripts
lmp -in in.lj_diffusion -var TEMP 1.0 -var DENS 0.8 > ../outputs/log.txt 2>&1
cd ../analysis
python3 analyze_diffusion.py
```

## References

- Meier et al. (2004) J. Chem. Phys. 121, 3671 (reference D* values)
- Allen & Tildesley (2017) Computer Simulation of Liquids, Ch. 2
- Thompson et al. (2022) Comput. Phys. Commun. 271, 108171 (LAMMPS)
