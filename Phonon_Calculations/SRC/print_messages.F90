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
    end if
end subroutine


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




subroutine error_message()
    use global_variables
    implicit none
    if (mpi_global%rank == 0) then
            write(6,'( A )') trim(adjustl(err_msg))
    end if
    return    
end subroutine


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

