program Phonon_frequencies

    use mpi
    use global_variables

    implicit none

    character(500) :: input_file
    integer :: i,j

    call mpi_init(mpierr)
    mpi_global%comm = MPI_COMM_WORLD
    call mpi_comm_size(mpi_global%comm, mpi_global%size_, mpierr)
    call mpi_comm_rank(mpi_global%comm, mpi_global%rank, mpierr)

    ! Read input file
    ! ---------------

    call get_command_argument(1,input_file)
    if (len_trim(input_file)==0) stop "ERROR! No input file entered"

    call read_input_file(input_file)
   
    if (no_of_groups.gt.mpi_global%size_) stop "ERROR! Invalid number of diagonalizations"

    ! Read LAMMPS data
    ! ----------------

    allocate(moire%mass(lammps_file%at_types))
    allocate(moire%real_pos(moire%natom,3))
    allocate(moire%crys(moire%natom,3))
    allocate(moire%at_types_i(moire%natom))

    call read_lammps_data()

    if (mpi_global%rank.eq.0) then
        write(*,*) "LAMMPS data read"
    end if

    ! Split the processors into the number of required groups. 
    ! Each group does a seperate diagonalization.
    ! --------------------------------------------------------
  
    allocate(blacs_grid%context(no_of_groups))
    call get_blacs_icontexts_single_q()
    call blacs_gridinfo( blacs_grid%context, blacs_grid%nprow, blacs_grid%npcol, &
                         blacs_grid%myprow,  blacs_grid%mypcol)

    if (mpi_global%rank==0) then
        write(*,*) "BLACS grid initiated"
    end if

    if (debug.eq.True) then
    write(*,'(8I6,X,4I8)') mpi_global%rank, mpi_global%color, mpi_global%rank,     &
                           mpi_global%size_, blacs_grid%nprow, blacs_grid%npcol,  &
                           blacs_grid%myprow, blacs_grid%mypcol, blacs_grid%context
    end if
    ! Read q points in each group 
    ! ---------------------------
    call distribution_length(q_file%nqpt, mpi_global%color-1, no_of_groups, st_q, en_q)
    write(*,*) st_q, en_q
    allocate(q_file%points((en_q-st_q+1),3))
    call read_q_file()
    call mpi_barrier(mpi_global%comm, mpierr)
    if (mpi_global%rank.eq.0) then
        write(*,*) "Q file read"
    end if
    write(*,*) q_file%nqpt
    call mpi_barrier(mpi_global%comm, mpierr)
    if (mpi_global%rank.eq.0) then
        write(*,'(3F10.2)') ((q_file%points(j,i),i=1,3),j=1,q_file%nqpt) 
    end if
    
    
   
    ! Read Force Constants
    ! --------------------
    
    ! Find size of local array to be allocated in each group

    call get_array_size(3*moire%natom, 3*moire%natom, fc%size_, fc%lld, fc%ncol)
    allocate(fc%a(fc%size_),STAT=allocate_status)
    if (allocate_status.ne.0) stop "Insufficient memory for force constants."

    ! Allocate array descriptors for force constants

    call descinit(fc%desca,3*moire%natom,3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, 0, 0, blacs_grid%context, fc%lld, info)  
    if (mpi_global%rank==0) then
        write(*,*) "Force array descriptor initiated"
    end if
    call read_force_constants()

    call mpi_barrier(mpi_global%comm,mpierr)

   
    ! Create dynamical matrix and diagonalize for every q points in a group
    ! ---------------------------------------------------------------------
    
    call get_array_size(3*moire%natom,3*moire%natom,dynmat%size_,dynmat%lld, dynmat%ncol)   
    allocate(dynmat%a(dynmat%size_),STAT=allocate_status)
    if (allocate_status.ne.0) stop "Insufficient memory for force constants."
    call descinit(dynmat%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, 0, 0, blacs_grid%context, dynmat%lld, info)
    
    call get_array_size(3*moire%natom,3*moire%natom,evec%size_,evec%lld,evec%ncol)
    allocate(evec%a(evec%size_),STAT=allocate_status)
    if (allocate_status.ne.0) stop "Insufficient memory for force constants."
    call descinit(evec%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb,   &
                  pzheevx_vars%nb, 0, 0, blacs_grid%context, evec%lld, info)
   
    allocate(W(3*moire%natom))
    do q_index = 1,(en_q-st_q)
        call create_dynamical_matrix()
        call mpi_barrier(mpi_global%comm, mpierr)
        if (mpi_global%rank.eq.0) then
            write(*,*) "Dynamical Matrix for q = ", q_index, "created "
        end if
        call diagonalize_matrix()
        call write_output()
        dynmat%a = cmplx(0.0,0.0)
        evec%a = cmplx(0.0,0.0)
        call mpi_barrier(mpi_global%comm, mpierr) 
    end do

    call blacs_exit(1) 
    call mpi_barrier(mpi_global%comm, mpierr) 
!    call mpi_comm_free(mpi_local%comm, mpierr)

    ! Deallocate arrays
    ! -----------------
    
    deallocate(W)
    deallocate(evec%a)
    deallocate(dynmat%a)
    deallocate(moire%mass)
    deallocate(moire%real_pos)
    deallocate(moire%crys)
    deallocate(moire%at_types_i)
    deallocate(fc%a)
    deallocate(q_file%points)
   
    call mpi_finalize(mpierr)

end program
