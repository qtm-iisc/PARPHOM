#!/bin/bash


# Cray
# ----

#module swap PrgEnv-cray PrgEnv-intel
#module load intel
#module load cray-hdf5

#ftn -o PhonFreq dynamical_matrix.f90 reading_lammps_data.f90 convert_to_ij.f90 create_dynamical_matrix.f90 read_input.f90 read_q_points.f90


# QTM
# ---

LIB='-L${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a -L${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -L${MKLROOT}/lib/intel64/libmkl_avx512.so -L${MKLROOT}/lib/intel64/libmkl_def.so -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl'

COMP='-I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include'



h5pfc -o PhonFreq dynamical_matrix.f90 reading_lammps_data.f90 convert_to_ij.f90 create_dynamical_matrix.f90 read_input.f90 read_q_points.f90 $LIB $COMP 

