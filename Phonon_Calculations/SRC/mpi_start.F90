subroutine initialize_mpi()

    use mpi
    use global_variables

    implicit none
    character(char_len) :: start_msg

    call mpi_init(mpierr)
    mpi_global%comm = MPI_COMM_WORLD
    call mpi_comm_size(mpi_global%comm, mpi_global%size_, mpierr)
    call mpi_comm_rank(mpi_global%comm, mpi_global%rank, mpierr)


    if (mpierr .ne. 0) then
        debug_str = 'MPI failed to initialize'    
        call debug_output(mpierr)
        call exit
    end if
    
    write(start_msg,'(A,I0,A)') "Program started with ", mpi_global%size_," MPI processes on "
    call date_time_message(trim(adjustl(start_msg)))

    return

end subroutine
