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
!------------------------------------------------------------------------------
!> \file mod_variables.F90
!> \brief Defines global variables, parameters, and derived types for PARPHOM.
!>
!> This module contains all global variables, parameters, and custom types used
!> throughout the PARPHOM package for phonon calculations. It centralizes the
!> definitions of BLACS/ScaLAPACK grid information, system properties, file
!> metadata, and MPI communication structures.
!>
!> - Provides derived types for BLACS grid, LAMMPS input, Brillouin zone points,
!>   ScaLAPACK variables, system structure, force constants, and distributed arrays.
!> - Contains global variables for distributed matrix storage, eigenvalues,
!>   eigenvectors, group velocities, and MPI communication.
!> - All variables and types are accessible via the `global_variables` module.
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This module must be used (`use global_variables`) in any subroutine or function
!>   that requires access to global parameters or data structures.
!>   All allocations and deallocations of distributed arrays are handled elsewhere.
!>
!> \warning
!>   Ensure that all variables are properly initialized before use.
!>   Incorrect initialization may lead to runtime errors or incorrect results.
!>
!> \remarks
!>   This file is part of the PARPHOM package for phonon calculations.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!------------------------------------------------------------------------------
module global_variables
    
    use hdf5

    integer , parameter :: char_len=2500
    character(len=char_len) :: err_msg, debug_str

    type blacs_info
#ifdef __QPOOL
        integer , allocatable, dimension(:) :: context
#else
        integer  :: context
#endif
        integer  :: nprow, npcol, myprow, mypcol, rank, size_
    end type blacs_info

    type lammps
        character(len=char_len) :: location, name_, atom_style
        integer  :: at_types
    end type lammps

    type bz_points
        character(len=char_len) :: location, name_
        integer  :: npt, start, finish
        double precision, allocatable, dimension(:,:) :: points
    end type bz_points

    type scalapack_variables
        integer  :: mb, nb, il, iu, comp_num_eval, comp_num_evec
        double precision :: vl, vu, abstol, orfac
        character(1) :: range_, comp_evec
    end type scalapack_variables

    type system
        integer  :: natom
        double precision, allocatable, dimension(:) :: mass
        double precision, allocatable, dimension(:,:) :: real_pos, crys, sup_cell_info
        integer , allocatable, dimension(:) :: at_types_i, lay_types
        double precision, dimension(3,3) :: lat , rec_lat
    end type system

    type comparr
        double complex, allocatable, dimension(:) :: mat
        integer, dimension(9) :: desca
        integer :: size_, lld, locq
    end type comparr


    type fc
        double precision, allocatable, dimension(:) :: mat
        integer, dimension(9) :: desca
        integer :: size_, lld, locq
        character(len=char_len) :: location, name_, dset_name
    end type fc

    

    type(blacs_info) :: grid
    type(lammps) :: lammps_file
    type(bz_points) :: q_file
    type(system) :: moire
    type(fc) :: force_const
    type(comparr), target :: dyn_mat, evec, vel
    logical :: evec_comp, comp_vel, print_progress
    character(len=1) :: vel_method
    double precision, allocatable, dimension(:) :: eval
    type(scalapack_variables) :: pzheevx_vars
    integer , parameter :: BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,     &
                          CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,    &
                          RSRC_ = 7, CSRC_ = 8, LLD_ = 9
    


    type mpi_group
        integer  :: comm
        integer  :: color
        integer  :: key
        integer  :: size_
        integer  :: rank
    end type mpi_group

    type(mpi_group) :: mpi_global
    integer  :: mpierr

#ifdef __QPOOL
    type(mpi_group) :: mpi_local
    integer  :: num_pools
#endif

    character(len=char_len) :: output_file_name, output_file_location

    integer , dimension(8) :: date_time
    character(len=char_len), parameter :: date_format = '(A,X,2(I0,A),I4,3(A,I2.2),A,SP,I0,A,SS,I0)'

end module
