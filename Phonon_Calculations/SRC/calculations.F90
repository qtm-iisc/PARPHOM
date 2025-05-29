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
!> \file calculations.F90
!> \brief Main driver for diagonalizing dynamical matrices and writing phonon calculation results.
!>
!> This file contains routines for iterating over q-points, constructing and diagonalizing the dynamical matrix,
!> and writing the resulting phonon properties to output files. It synchronizes all MPI processes and ensures
!> that results are written in a distributed parallel environment.
!>
!> - Loops over all q-points in the Brillouin zone
!> - Calls routines to construct the dynamical matrix, diagonalize it, and write results
!> - Uses MPI barriers to synchronize processes
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This file is part of the PARPHOM package for phonon calculations.
!>   All routines assume that global variables and MPI have been initialized.
!>
!> \warning
!>   Ensure that all input files and global variables are set before calling these routines.
!>   Incorrect setup may lead to runtime errors or incorrect results.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain

subroutine diagonalize_and_write()
    use global_variables
    use mpi
    implicit none  
    
    integer :: q_loc    

    do q_loc = q_file%start, q_file%finish
        call create_dynamical_matrix(q_loc,0)
        call diagonalize_dynamical_matrix()
        call write_output(q_loc)    
    end do  

    call mpi_barrier(mpi_global%comm, mpierr)
    
    return

end subroutine
