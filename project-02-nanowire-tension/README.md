# Tensile Deformation of Copper Nanowire via MD

## Summary

Simulated uniaxial tensile deformation of a [001] copper nanowire
at 300 K using classical MD with EAM potential. Extracted Young's
modulus (51.5 GPa) and yield stress (4.6 GPa). Results consistent
with published nanowire MD studies and size-dependent strengthening.

## Method

- **Software:** LAMMPS (serial)
- **Potential:** EAM (Cu_u3.eam)
- **System:** 3804 atoms, cylindrical wire, radius ~18 A, axis [001]
- **Boundaries:** s s p (free surfaces in x,y; periodic along wire z)
- **Protocol:**
  1. Energy minimization (surface relaxation)
  2. NVT equilibration at 300 K (20 ps)
  3. Constant strain rate deformation (erate = 0.001/ps = 10^9 /s)
  4. Total strain: 20%, 200,000 steps

## Key Results

| Property | MD Result | Bulk Cu (expt) | Notes |
|----------|-----------|----------------|-------|
| Young's modulus | 51.5 GPa | ~67 GPa [100] | 23% lower due to surface effects |
| Yield stress | 4.6 GPa | ~0.07 GPa | 65x higher: no pre-existing dislocations |
| Yield strain | 9.1% | ~0.1% | Pristine crystal resists until surface nucleation |

## Discussion

The reduced modulus compared to bulk is expected for nanowires
where a significant fraction of atoms are under-coordinated surface
atoms with lower stiffness. The dramatically elevated yield stress
reflects the absence of pre-existing dislocations — the crystal
must nucleate them from the free surface, requiring stress near the
theoretical shear strength. The serrated plastic flow region shows
repeated dislocation burst events typical of nanoscale deformation.

The strain rate (10^9 /s) is ~12 orders of magnitude faster than
typical experiments (~10^-3 /s). This artificially elevates both
yield stress and flow stress. Lower strain rates would require
prohibitively long simulations but would give lower yield stresses.

## Files

- `scripts/in.cu_nanowire` — LAMMPS input with deformation protocol
- `analysis/analyze_tension.py` — Systematic E extraction + plotting
- `figures/stress_strain_curve.png` — Publication figure
- `outputs/stress_strain.dat` — Raw stress-strain data
- `outputs/dump.tension.*` — Trajectory for OVITO visualization

## How to Reproduce

```bash
cd scripts
cp /usr/share/lammps/potentials/Cu_u3.eam .
lmp -in in.cu_nanowire > ../outputs/log_tension.txt 2>&1
cd ../analysis
python3 analyze_tension.py
```

## References

- Thompson et al. (2022) Comput. Phys. Commun. 271, 108171 (LAMMPS)
- Foiles, Baskes, Daw (1986) Phys. Rev. B 33, 7983 (EAM)
- Pal & Reddy (2024) Molecular Dynamics for Materials Modeling (Routledge)
