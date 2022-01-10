subroutine allocate_distributed_arrays()

    use global_variables
    use mpi
    implicit none

    integer :: rsrc, csrc, info
    integer, external :: numroc
#ifdef __DEBUG
    integer :: i
#endif

    rsrc = 0
    csrc = 0

    force_const%locq = numroc(3*moire%natom,pzheevx_vars%mb,grid%mypcol,csrc,grid%npcol)
    force_const%locq = max(1,force_const%locq)
    force_const%lld  = numroc(3*moire%natom,pzheevx_vars%nb,grid%myprow,rsrc,grid%nprow)
    force_const%lld  = max(1,force_const%lld)

    force_const%size_ = force_const%lld*force_const%locq
    allocate(force_const%mat(force_const%size_))

    write(debug_str,'(A)') "Force Constant Matrix Allocated"

    call debug_output(0)

    call descinit(force_const%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, force_const%lld, info)
#ifdef __DEBUG
    write(debug_str,'(A)') "\r\nForce Constant Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", force_const%size_, &
                                   " Local Rows: ", force_const%lld, &
                                   " Local Cols: ", force_const%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nDynamical Matrix descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", force_const%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", force_const%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", force_const%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", force_const%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", force_const%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", force_const%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", force_const%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", force_const%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", force_const%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do


#endif

    dyn_mat%locq = numroc(3*moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    dyn_mat%locq = max(dyn_mat%locq,1)
    dyn_mat%lld = numroc(3*moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    dyn_mat%lld = max(dyn_mat%lld,1)

    dyn_mat%size_ = dyn_mat%locq*dyn_mat%lld

    allocate(dyn_mat%mat(dyn_mat%size_))
    write(debug_str,'(A)') "Dynamical Matrix Allocated"
    call debug_output(0)

    call descinit(dyn_mat%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, dyn_mat%lld, info)

#ifdef __DEBUG
    write(debug_str,'(A)') "\r\nDynamical Matrix Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", dyn_mat%size_, &
                                   " Local Rows: ", dyn_mat%lld, &
                                   " Local Cols: ", dyn_mat%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nDynamical Matrix descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", dyn_mat%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", dyn_mat%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", dyn_mat%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", dyn_mat%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", dyn_mat%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", dyn_mat%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", dyn_mat%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", dyn_mat%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", dyn_mat%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif
    
    evec%locq = numroc(3*moire%natom, pzheevx_vars%mb, grid%mypcol, csrc, grid%npcol)
    evec%locq = max(evec%locq,1)
    evec%lld = numroc(3*moire%natom, pzheevx_vars%nb, grid%myprow, rsrc, grid%nprow)
    evec%lld = max(evec%lld,1)

    evec%size_ = evec%locq*evec%lld
    call descinit(evec%desca, 3*moire%natom, 3*moire%natom, pzheevx_vars%mb, &
                  pzheevx_vars%nb, rsrc, csrc, grid%context, evec%lld, info)
    
    if (pzheevx_vars%comp_evec=='V') then
        allocate(evec%mat(evec%size_))
        write(debug_str,'(A)') "Eigenvectors Allocated"
        call debug_output(0)

#ifdef __DEBUG
    write(debug_str,'(A)') "\r\nEigenvector Parameters: "
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A,I8,3(A,I0))') "Rank: ", mpi_global%rank, &
                                   "    |  Size: ", evec%size_, &
                                   " Local Rows: ", evec%lld, &
                                   " Local Cols: ", evec%locq
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str, '(A)') "\r\nEigenvector descriptor:"
    call debug_output(0)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(*,'(A)') " ---------"
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  DTYPE_: ", evec%desca(DTYPE_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CTXT_: ", evec%desca(CTXT_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  M_: ", evec%desca(M_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  N_: ", evec%desca(N_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  MB_: ", evec%desca(MB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  NB_: ", evec%desca(NB_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  RSRC_: ", evec%desca(RSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  CSRC_: ", evec%desca(CSRC_)
            write(*,'(A,I8,A,I0)') "Rank: ", mpi_global%rank, &
                                   "    |  LLD_: ", evec%desca(LLD_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif
    else
        allocate(evec%mat(1))
    end if

    
    allocate(eval(3*moire%natom))
    write(debug_str,'(A)') "Eigenvalue Allocated"
    call debug_output(0)


end subroutine
