#!/bin/bash


# Cray
# ----

module swap PrgEnv-cray PrgEnv-intel/6.0.9
module load intel
#module load cray-hdf5
module load cray-hdf5-parallel/1.12.0.0

ftn -o PhonFreq_cray dynamical_matrix.f90 reading_lammps_data.f90 convert_to_ij.f90 create_dynamical_matrix.f90 read_input.f90 read_q_points.f90


# QTM
# ---

#LIB='-L${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a -L${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -L${MKLROOT}/lib/intel64 -L${MKLROOT}/lib/intel64/libmkl_avx512.so -L${MKLROOT}/lib/intel64/libmkl_def.so -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl'

#COMP='-I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include'

#h5pfc -o PhonFreq_QTM dynamical_matrix.f90 reading_lammps_data.f90 convert_to_ij.f90 create_dynamical_matrix.f90 read_input.f90 read_q_points.f90 $LIB $COMP 

