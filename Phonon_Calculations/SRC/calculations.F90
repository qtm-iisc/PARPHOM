subroutine diagonalize_and_write()
    use global_variables
    use mpi
    implicit none  
    
    integer :: q_loc    

    do q_loc = q_file%start, q_file%finish
        call create_dynamical_matrix(q_loc)
        call diagonalize_dynamical_matrix()
        call write_output(q_loc)    
    end do  

    call mpi_barrier(mpi_global%comm, mpierr)
    
    return

end subroutine
