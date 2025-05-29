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
!> \file mpi_end.F90
!> \brief Finalizes the MPI environment and outputs completion messages for PARPHOM.
!>
!> This subroutine synchronizes all processes, prints a completion message with a timestamp, and finalizes the MPI environment.
!> If MPI fails to close, an error message is printed and the program exits.
!>
!> - Calls MPI barrier to synchronize all ranks
!> - Prints completion message with date and time
!> - Finalizes MPI and handles errors
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This file is part of the PARPHOM package for phonon calculations.
!>
!> \warning
!>   Ensure that all calculations and file outputs are complete before calling this routine.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!

subroutine close_mpi()

    use mpi
    use global_variables

    implicit none

    call mpi_barrier(mpi_global%comm, mpierr)
    call date_time_message('Calculations completed on ')     
    call mpi_finalize(mpierr)
    
    if (mpierr .ne. 0) then
        debug_str = 'MPI failed to close'    
        call debug_output(mpierr)
        call exit
    end if

    return

end subroutine
