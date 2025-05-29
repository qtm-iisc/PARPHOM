! Package: PARPHOM
! Authors: Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
! License: GPL-3.0
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! Program: compute_phonon_bands (Main driver for phonon dispersion and band structure calculations)
!
!> \file   main.F90
!> \brief  Main program to orchestrate phonon dispersion and electronic band structure calculations.
!> \details
!> This program initializes MPI and BLACS, reads inputs (structure, q-points, force constants),
!> distributes data across processes, diagonalizes dynamical matrices, and writes phonon and electronic band outputs.

program compute_phonon_bands
  implicit none

  !> Initialize MPI environment: set up communicator, ranks, and sizes
  call initialize_mpi()

  !> Read user-defined input file and populate global parameters
  call read_input()

  !> Set up the BLACS grid for parallel linear algebra routines
  call blacs_grid_initialization()

  !> Load atomic structure from LAMMPS data files into distributed arrays
  call read_lammps_data()

  !> Read list of q-points for phonon dispersion calculations
  call read_q_file()

  !> Allocate arrays for dynamical matrices and wavefunctions across BLACS grid
  call allocate_distributed_arrays()

  !> Read precomputed force constants for dynamical matrix construction
  call read_force_constants()

  !> Diagonalize dynamical matrices and output phonon and electronic band structures
  call diagonalize_and_write()

  !> Finalize MPI and release all resources
  call close_mpi()

end program compute_phonon_bands
