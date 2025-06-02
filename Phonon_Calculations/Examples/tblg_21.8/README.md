# Example: 21.8° Twisted Bilayer Graphene Phonon Calculations

This folder contains input and output files for phonon calculations of twisted bilayer graphene (tBLG) with a twist angle of 21.8°. The calculations are performed to generate phonon band structures along high-symmetry paths and to sample the entire Brillouin zone for density of states (DOS) computations.

## Contents

- **in.phonon**: Input file for phonon calculations using the REBO potential.
- **in.phonon_tersoff**: Input file for phonon calculations using the Tersoff potential.
- **lammps.dat**: LAMMPS data file for the tBLG structure (REBO potential).
- **lammps.dat_tersoff**: LAMMPS data file for the tBLG structure (Tersoff potential).
- **C.drip, CH.rebo, BNC.tersoff**: Potential parameter files for LAMMPS.

## Workflow

1. **Generate FORCE_CONSTANTS**: Run the force constant generation script or workflow to produce the `FORCE_CONSTANTS` file in this folder. This is the first and essential step for all subsequent calculations.
2. **Generate q-point files**: Create the required q-point files for band structure and DOS calculations (e.g., `q_points_band.dat`, `q_points.dat`, etc.).
3. **Prepare Input Files**: Use the provided LAMMPS data and potential files to set up the system. Refer to `inp.dat` for creating the input files for `parphom.x`.
4. **Run Phonon Calculations**: Execute the phonon calculation using the appropriate input file (`in.phonon` or `in.phonon_tersoff`).
5. **Band Structure**: Calculate the phonon dispersion along a high-symmetry path in the Brillouin zone to obtain the phonon band structure.
6. **Brillouin Zone Sampling**: Use a dense grid to sample the entire Brillouin zone for DOS calculations.
7. **For further analysis, you can use the utility routines.**

## Purpose

- **Band Structure**: Compute the phonon dispersion along high-symmetry directions.
- **DOS Calculation**: Compute the phonon density of states by sampling the full Brillouin zone or an irreducible wedge of the Brillouin zone.

## Notes
- Ensure all required potential files are present before running calculations.
- The workflow is compatible with LAMMPS and the phonon calculation tools used in this project.

---
For more details, refer to the main project documentation or contact the authors.
