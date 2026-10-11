# DMSO GAFF2 topology validation

## Purpose

Validate the conversion of a single DMSO molecule from AMBER/GAFF2
topology to GROMACS topology by comparing single-point potential-energy
terms using the same molecular coordinates.

This is a topology-conversion check, not a validation of bulk DMSO
thermodynamic properties or transport properties.

## Molecular model

- Molecule: dimethyl sulfoxide (DMSO), C2H6OS
- Atom count: 10
- Bond count: 9
- Force field: AMBER GAFF2
- Partial charges: Antechamber AM1-BCC
- Additional parameters: generated using parmchk2
- AMBER topology and coordinates: dmso.prmtop, dmso.inpcrd
- GROMACS topology and coordinates: dmso_gmx.top, dmso_gmx.gro

The MOL2 partial charges were adjusted to make their rounded sum exactly
zero at the stored precision. This adjustment was applied to the first
carbon atom. It is a numerical charge-neutralization correction and
should not be interpreted as a re-derived or independently validated
partial-charge model.

## Validation procedure

1. Generate the GAFF2 MOL2 file using Antechamber and AM1-BCC charges.
2. Generate the additional parameter file using parmchk2.
3. Build the AMBER topology with tleap.
4. Convert the topology to GROMACS format using ParmEd.
5. Run single-point energy calculations with AMBER sander and GROMACS.
6. Compare bonded and nonbonded energy terms after unit conversion.

The precise-coordinate comparison uses coordinates converted from the
AMBER inpcrd file to nm, written with five decimal places, and translated
by +5 nm in each dimension in a 10 nm box. This avoids the larger bonded
energy discrepancy caused by rounding coordinates to 0.001 nm in a
standard GRO file.

## Energy comparison

AMBER energies were reported in kcal/mol and converted to kJ/mol using
1 kcal/mol = 4.184 kJ/mol.

| Energy term | AMBER (kJ/mol) | GROMACS (kJ/mol) |
|---|---:|---:|
| Bond | 2.191161 | 2.19137 |
| Angle | 1.518374 | 1.51844 |
| Proper dihedral | 0.401246 | 0.401219 |
| 1-4 Lennard-Jones | -0.372376 | -0.372497 |
| 1-4 Coulomb | -118.972040 | -118.976 |
| Short-range Lennard-Jones | -0.340578 | -0.340426 |
| Short-range Coulomb | 25.869672 | 25.8707 |

The listed terms agree closely for this single-molecule comparison.
Small differences may arise from output precision and implementation
details. The comparison does not establish equivalence for all possible
configurations or simulation conditions.

## Important limitations

- The model has not been validated against experimental bulk DMSO
  density, dielectric properties, diffusion, or other liquid properties.
- The single-molecule calculation does not test periodic bulk-liquid
  behavior.
- The energy comparison is specific to the recorded coordinates and
  settings.
- The numerical charge correction should be reviewed before using this
  model for production simulations.
- The GROMACS single-point setup uses a 4.0 nm cutoff; the AMBER
  single-point reference setup may not use identical nonbonded cutoff
  conventions. Interpret nonbonded comparisons in light of the actual
  input settings.

## Included files

- `dmso_seed.mol2`, `dmso_test.sdf`: starting structure files
- `dmso_gaff2.mol2`: GAFF2 atom types and partial charges
- `dmso.frcmod`: additional AMBER parameters
- `dmso.prmtop`, `dmso.inpcrd`, `dmso.pdb`: AMBER topology and coordinates
- `dmso_gmx.top`, `dmso_gmx.gro`: GROMACS topology and coordinates
- `dmso_sp.out`: AMBER single-point output
- `dmso_gmx_sp_precise.log`: GROMACS precise-coordinate single-point log
- `dmso_gmx_sp_precise_mdrun.log`: GROMACS run log
- `dmso_energy_precise.xvg`: GROMACS energy data
- `leap.in`, `leap.log`, `antechamber.log`, `parmchk2.log`: parameterization logs
- `*.mdp`, `*grompp.log`, `grompp_check_pass.log`: GROMACS input and preprocessing logs
