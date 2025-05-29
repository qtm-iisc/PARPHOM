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
!> \file mpi_start.F90
!> \brief Initializes MPI environment, sets communicator, rank, and size, and prints start message.
!>
!> This subroutine initializes the MPI environment for PARPHOM, sets up the global communicator, determines the number of processes and the rank of each process, checks for initialization errors, and prints a startup logo and date-time message.
!>
!> - Calls `mpi_init` to initialize MPI
!> - Sets `mpi_global%comm` to `MPI_COMM_WORLD`
!> - Retrieves the total number of processes and the rank of the current process
!> - Checks for initialization errors and outputs debug message if any
!> - Prints a startup logo and a date-time stamped message indicating the number of MPI processes
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This file is part of the PARPHOM package for phonon calculations.
!>
!> \warning
!>   Ensure that MPI is not already initialized before calling this routine.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> 
!> \subroutine initialize_mpi
!> \brief  Initializes the MPI environment and sets up global variables.
!>
!> This routine:
!>   - Calls `mpi_init` to initialize MPI.
!>   - Sets `mpi_global%comm` to `MPI_COMM_WORLD`.
!>   - Retrieves the total number of processes and the rank of the current process.
!>   - Checks for initialization errors and outputs debug message if any.
!>   - Prints a startup logo and a date-time stamped message indicating the number of MPI processes.
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This file is part of the PARPHOM package for phonon calculations.
!>
!> \warning
!>   Ensure that MPI is not already initialized before calling this routine.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
subroutine initialize_mpi()
  use mpi
  use global_variables
  implicit none
  character(char_len) :: start_msg

  !> Initialize the MPI environment
  call mpi_init(mpierr)
  mpi_global%comm = MPI_COMM_WORLD

  !> Get the total number of MPI processes
  call mpi_comm_size(mpi_global%comm, mpi_global%size_, mpierr)

  !> Get the rank of this MPI process
  call mpi_comm_rank(mpi_global%comm, mpi_global%rank, mpierr)

  !> Check for MPI initialization errors
  if (mpierr .ne. 0) then
    debug_str = 'MPI failed to initialize'
    call debug_output(mpierr)
    call exit
  end if

  !> Print the program start logo
  call print_start_logo()
  write(start_msg,'(A,I0,A)') "Program started with ", mpi_global%size_," MPI processes on "
  !> Print date and time with start message
  call date_time_message(trim(adjustl(start_msg)))

  return
end subroutine initialize_mpi
