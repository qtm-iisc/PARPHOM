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
!> \file q_read.F90
!> \brief Reads and distributes q-point lists for phonon calculations in PARPHOM.
!>
!> This file contains routines for reading the list of q-points from input files and distributing them
!> across MPI processes for parallel phonon calculations. It handles both single-pool and multi-pool
!> (QPOOL) parallelization strategies and provides debug output for q-point distribution.
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This file is part of the PARPHOM package for phonon calculations.
!>
!> \warning
!>   Ensure that the q-point file exists and is formatted correctly before running.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain

!> \brief Reads the q-point list from file and distributes q-points among MPI processes.
!> \details
!> This subroutine reads the list of q-points from the specified file, determines the range of q-points
!> assigned to each process (or pool), and allocates the q-point array accordingly. It provides debug output
!> for the distribution and contents of q-points read by each process or pool.
subroutine read_q_file()
    use global_variables
    implicit none
    integer :: i, error
    character(len=char_len) :: file_name_ , temp
#ifdef __DEBUG
    integer :: j
#endif
#ifdef __QPOOL

    call distribution_length(q_file%npt, mpi_local%color-1, num_pools, &
                             q_file%start, q_file%finish)
#else
    call distribution_length(q_file%npt, 0, 1, q_file%start, q_file%finish)
#endif


#ifdef __DEBUG
    write(debug_str, '(A)') '\r\nRank      Start Q read        End Q read'
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (mpi_global%rank==i) then
            write(*,'(3I8)') i, q_file%start, q_file%finish
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif


    write(file_name_,'(2A)') trim(adjustl(q_file%location)),trim(adjustl(q_file%name_))

    open(unit=25, file=trim(adjustl(file_name_)),action='read',iostat=error)
    
    if (error.ne.0) then
        write(err_msg, '(2A)') 'Error reading ', trim(adjustl(file_name_))
        call error_message()
        call exit
    end if

    allocate(q_file%points(q_file%finish-q_file%start+1,3))

    if (q_file%start == 1) then
        do i=q_file%start,q_file%finish
            read(25,*) q_file%points(i,1), q_file%points(i,2), q_file%points(i,3)
        end do
    else
        do i=1,q_file%start-1
            read(25,*) temp
        end do
        do i=1,q_file%finish-q_file%start+1
            read(25,*) q_file%points(i,1), q_file%points(i,2), q_file%points(i,3)
        end do
    end if
    close(unit=25)

#ifdef __DEBUG
#ifdef __QPOOL
    debug_str = '\r\n\r\nLocal Q points read in each group'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=1,num_pools
        if (mpi_local%color == i) then
            if (mpi_local%rank == 0) then
                do j=1,q_file%finish-q_file%start+1
                    write(*,'(I0,3F16.6)') i, q_file%points(j,1), &
                                              q_file%points(j,2), &
                                              q_file%points(j,3)
                end do
            end if
        end if
    end do
#else
    debug_str = '\r\n\r\nQ points read from file'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    if (mpi_global%rank==0) then
        do j = 1, q_file%finish-q_file%start+1
            write(*,'(I0,3F16.6)')  j, q_file%points(j,1), &
                                       q_file%points(j,2), &
                                       q_file%points(j,3)
        end do
    end if
    call mpi_barrier(mpi_global%comm, mpierr)
    
#endif
#endif


    return

end subroutine
