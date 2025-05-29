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
!> \file print_messages.F90
!> \brief Utility routines for printing messages, errors, and timestamps in PARPHOM.
!>
!> This file contains subroutines for printing the program logo, debug messages, error messages,
!> and date/time information in a parallel MPI environment. All output is managed so that only the root
!> process prints to the screen, ensuring clean output in distributed runs.
!>
!> \author Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
!> \ingroup phonon_allocation
!>
!> \note
!>   This file is part of the PARPHOM package for phonon calculations.
!>
!> \copyright GPL-3.0 Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain

!> \brief Prints the PARPHOM program start logo and author credits (root process only).
subroutine print_start_logo()
    use global_variables
    implicit none
    if (mpi_global%rank == 0) then
        write(6,*) " "
        write(6,*) " "
        write(6,'(A)') "██████╗  █████╗ ██████╗ ██████╗ ██╗  ██╗ ██████╗ ███╗   ███╗"
        write(6,'(A)') "██╔══██╗██╔══██╗██╔══██╗██╔══██╗██║  ██║██╔═══██╗████╗ ████║"
        write(6,'(A)') "██████╔╝███████║██████╔╝██████╔╝███████║██║   ██║██╔████╔██║"
        write(6,'(A)') "██╔═══╝ ██╔══██║██╔══██╗██╔═══╝ ██╔══██║██║   ██║██║╚██╔╝██║"
        write(6,'(A)') "██║     ██║  ██║██║  ██║██║     ██║  ██║╚██████╔╝██║ ╚═╝ ██║"
        write(6,'(A)') "╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝"
        write(6,*) " "
        write(6,*) " "
        write(6,'(A)') "               Phonon Spectrum Calculator v1.0"
        write(6,*) " "
        write(6,*) " "
        write(6,'(A)') "        S. Mandal, I.Maity, H R Krishnamurthy, M. Jain"
        write(6,*) " "
        write(6,'(A)') "        Please cite the following paper for any publication using this code:"
        write(6,'(A)') "        S. Mandal, I. Maity, H. R. Krishnamurthy, M. Jain, arXiv:2410.21075"
        write(6,*) " "
        write(6,*) " "
        write(6,'(A)') "        This program is free software: you can redistribute it and/or modify"
        write(6,'(A)') "        it under the terms of the GNU General Public License as published by"
        write(6,'(A)') "        the Free Software Foundation, either version 3 of the License, or"
        write(6,'(A)') "        (at your option) any later version."
        write(6,*) " "
        write(6,'(A)') "        This program is distributed in the hope that it will be useful,"
        write(6,'(A)') "        but WITHOUT ANY WARRANTY; without even the implied warranty of"
        write(6,'(A)') "        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
        write(6,'(A)') "        GNU General Public License for more details."
        write(6,*) " "
        write(6,'(A)') "        You should have received a copy of the GNU General Public License"
        write(6,'(A)') "        along with this program.  If not, see <https://www.gnu.org/licenses/>."
        write(6,*) " "
        write(6,*) " "
    end if  
end subroutine

!> \brief Prints a debug or error message, optionally with an error code (root process only).
!> \param code Integer error code. If zero, prints only the debug string; otherwise, prints error code as well.
subroutine debug_output(code)
    use global_variables
    implicit none
    integer, intent(in) :: code 
    if (mpi_global%rank == 0) then
        if (code.ne.0) then
            write(6,'( 2A, I5 )') trim(adjustl(debug_str)),' Error Code = ', code
        else
            write(6,'( A )') trim(adjustl(debug_str))
        end if
    end if
    return    
end subroutine

!> \brief Prints the current error message string (root process only).
subroutine error_message()
    use global_variables
    implicit none
    if (mpi_global%rank == 0) then
            write(6,'( A )') trim(adjustl(err_msg))
    end if
    return    
end subroutine

!> \brief Prints a message with the current date and time (root process only).
!> \param input_str The message to print before the date/time.
subroutine date_time_message(input_str)
    use global_variables
    implicit none

    character(len=*), intent(in) :: input_str

    call date_and_time(VALUES=date_time)

    write(debug_str, trim(date_format)) trim(adjustl(input_str)), date_time(3), &
                                        "/", date_time(2), "/", date_time(1),   &
                                        " at ", date_time(5), ':',date_time(6), &
                                        ':',date_time(7), " UTC: ", &
                                        date_time(4)/60, ':', mod(date_time(4),60)
    call debug_output(0)
    return
end subroutine

#ifdef __QPOOL
!> \brief Prints a message with the current date and time for local MPI pool (pool root only).
!> \param input_str The message to print before the date/time.
subroutine date_time_message_local(input_str)
    use global_variables
    implicit none
    
    character(len=*), intent(in) :: input_str
    character(len=2000) :: local_output

    call mpi_barrier(mpi_local%comm, mpierr)

    call date_and_time(VALUES=date_time)

    write(local_output,trim(date_format)) trim(adjustl(input_str)), date_time(3), &
                                        "/", date_time(2), "/", date_time(1),   &
                                        " at ", date_time(5), ':',date_time(6), &
                                        ':',date_time(7), " UTC: ", &
                                        date_time(4)/60, ':', mod(date_time(4),60)

    if (mpi_local%rank==0) then 
        write(*,*) trim(local_output)
    end if

    call mpi_barrier(mpi_local%comm, mpierr)

    return 
end subroutine
#endif

