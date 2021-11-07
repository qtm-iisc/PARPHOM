subroutine diagonalize_matrix()

    use mpi
    use global_variables

    implicit none
    integer :: ia, ja
    double precision, allocatable, dimension(:) :: gap, rwork
    double complex, allocatable, dimension(:) :: work
    integer, allocatable, dimension(:) :: iwork, ifail, iclustr
    integer :: lwork, liwork, lrwork, lgap, lifail, liclustr
    integer, parameter :: abstol = -1
    integer, parameter :: orfac = -1
    double precision :: st_time, en_time

    st_time = MPI_WTIME()

    ia = 1
    ja = 1
    lgap = blacs_grid%nprow*blacs_grid%npcol
    lifail = 3*moire%natom
    liclustr = 2*lgap
    
    allocate(gap(lgap))
    allocate(ifail(lifail))
    allocate(iclustr(liclustr))
    lwork = -1
    lrwork = -1
    liwork = -1
    allocate(work(1))
    allocate(rwork(1))
    allocate(iwork(1))
    
    if (pzheevx_vars%comp_evec) then
        call pzheevx('V', pzheevx_vars%range_, 'U', 3*moire%natom, dynmat%a,   &
                     ia, ja, dynmat%desca, pzheevx_vars%vl_, pzheevx_vars%vu_, &
                     pzheevx_vars%il_, pzheevx_vars%iu_, abstol,               &
                     pzheevx_vars%num_eval_comp,  pzheevx_vars%num_evec_comp,  &
                     W, orfac, evec%a, 1, 1, evec%desca, work, lwork, rwork,   &
                     lrwork, iwork, liwork, ifail, lifail, iclustr, gap, info)
    else
        call pzheevx('N', pzheevx_vars%range_, 'U', 3*moire%natom, dynmat%a,   &
                     ia, ja, dynmat%desca, pzheevx_vars%vl_, pzheevx_vars%vu_, &
                     pzheevx_vars%il_, pzheevx_vars%iu_, abstol,               &
                     pzheevx_vars%num_eval_comp,  pzheevx_vars%num_evec_comp,  &
                     W, orfac, evec%a, 1, 1, evec%desca, work, lwork, rwork,   &
                     lrwork, iwork, liwork, ifail, lifail, iclustr, gap, info)
    end if

    lwork = int(work(1)+1)
    lrwork = int(rwork(1)+1)
    liwork = int(iwork(1)+1)

    deallocate(work)
    deallocate(iwork)
    deallocate(rwork)

    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(rwork(lrwork))

    if (pzheevx_vars%comp_evec) then
        call pzheevx('V', pzheevx_vars%range_, 'U', 3*moire%natom, dynmat%a,   &
                     ia, ja, dynmat%desca, pzheevx_vars%vl_, pzheevx_vars%vu_, &
                     pzheevx_vars%il_, pzheevx_vars%iu_, abstol,               &
                     pzheevx_vars%num_eval_comp,  pzheevx_vars%num_evec_comp,  &
                     W, orfac, evec%a, 1, 1, evec%desca, work, lwork, rwork,   &
                     lrwork, iwork, liwork, ifail, lifail, iclustr, gap, info)
    else
        call pzheevx('N', pzheevx_vars%range_, 'U', 3*moire%natom, dynmat%a,   &
                     ia, ja, dynmat%desca, pzheevx_vars%vl_, pzheevx_vars%vu_, &
                     pzheevx_vars%il_, pzheevx_vars%iu_, abstol,               &
                     pzheevx_vars%num_eval_comp,  pzheevx_vars%num_evec_comp,  &
                     W, orfac, evec%a, 1, 1, evec%desca, work, lwork, rwork,   &
                     lrwork, iwork, liwork, ifail, lifail, iclustr, gap, info)
    end if
   
    if (info.ne.0) then
        write(*,*) "PZHEEVX completed with info = ", info
    end if

    if (info.gt.0) then
        if (mod(info,2).ne.0) then
            write(*,*) "One or more eigenvalues failed to converge"
        end if
        if (mod(info/2, 2).ne.0) then
            write(*,*) "Eigenvectors corresponding to the following indices could not be orthogonalized"
            write(*,*) iclustr
        end if
        if (mod(info/4, 2).ne.0) then
            write(*,*) "All eigenvectors in the specified interval could not be computed due to insufficient memory"
        end if
        if (mod(info/8,2).ne.0) then
            write(*,*) "One or more eigenvalues were not computed"
        end if
    end if

    deallocate(work)
    deallocate(iwork)
    deallocate(rwork)
    deallocate(gap)
    deallocate(ifail)
    deallocate(iclustr)

    en_time = MPI_WTIME()

    if (mpi_global%rank==0) then
        write(*,'(A,I4,A,F10.4,A)') "Diagoanlizations completed in group ",mpi_global%color, &
                   " in ", en_time-st_time, " sec."
    end if

    return

end subroutine
