#!/bin/bash

module swap PrgEnv-cray PrgEnv-intel
module load intel
module load cray-hdf5

ftn -o PhonFreq dynamical_matrix.f90 reading_lammps_data.f90 convert_to_ij.f90 create_dynamical_matrix.f90 read_input.f90 read_q_points.f90
