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
