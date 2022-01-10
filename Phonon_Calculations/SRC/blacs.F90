subroutine blacs_grid_initialization()

    use global_variables
    use mpi

    implicit none
#ifdef __DEBUG
    integer :: i
#endif
#ifdef __KPOOL
    integer, allocatable, dimension(:,:) :: map, dim_
    integer, allocatable, dimension(:) :: icontxts
    integer :: global_icontxt, st_group, j, k
    character(500) :: format_
#endif

#ifdef __KPOOL

    ! With KPools
    ! -----------
    !  
    ! Split the processors in groups with each group handling a k-point
    ! The number of pools is taken as an input from the user : num_pools
    ! 
    ! mpi_local%color = group number (1,2, ...)
    ! mpi_local%key = rank of a process within the group
    ! 
   
    if (num_pools.gt.mpi_global%size_) then
        debug_str = 'Invalid k-pools. k-pools set at present is the value of Error code.'
        call debug_output(num_pools)
        call close_mpi()
        call exit
    end if

    call get_color(mpi_global%rank, mpi_global%size_, num_pools,  &
                   mpi_local%color, mpi_local%key, mpi_local%size_)
    call mpi_comm_split(mpi_global%comm, mpi_local%color, mpi_local%key, &
                        mpi_local%comm, mpierr)
    call mpi_comm_rank(mpi_local%comm, mpi_local%rank, mpierr)

#ifdef __DEBUG
    call mpi_barrier(mpi_global%comm, mpierr)
    debug_str = ' '
    do i=1,2
        call debug_output(0)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)
    debug_str = 'MPI split completed.'
    call debug_output(mpierr)
    debug_str = 'Local Group Information'
    call debug_output(0)
    debug_str = 'GlobRank  GlobSize   LocRank   LocSize     Color       Key'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do k=0,mpi_global%size_-1
        if (mpi_global%rank==k) then
            write(*,'( I8,5I10 )') mpi_global%rank, mpi_global%size_, mpi_local%rank, & 
                                   mpi_local%size_, mpi_local%color, mpi_local%key
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif
    
    allocate(dim_(num_pools,2))
    allocate(icontxts(num_pools))

    dim_= 0
    call get_optimum_grid(mpi_local%size_, grid%npcol, grid%nprow)

    if (mpi_local%rank==0) then
        dim_(mpi_local%color, 1) = grid%nprow
        dim_(mpi_local%color, 2) = grid%npcol
    end if

    call mpi_allreduce(MPI_IN_PLACE, dim_, num_pools*2, MPI_INT, MPI_SUM, &
                       mpi_global%comm, mpierr)

    call blacs_get(-1,0,global_icontxt)
    icontxts = global_icontxt

    st_group = 0
#ifdef __DEBUG    
    call mpi_barrier(mpi_global%comm, mpierr)
    debug_str = ' '
    do i=1,2
        call debug_output(0)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)
    debug_str = 'GlobRank    Icontxt     Group'
    call debug_output(0)
#endif
    do i=1,num_pools
        allocate(map(dim_(i,1),dim_(i,2)))
        do j=1,dim_(i,1)
            do k=1,dim_(i,2)
                map(j,k) = st_group
                st_group = st_group+1
            end do
        end do
        call blacs_gridmap(icontxts(i),map,dim_(i,1), dim_(i,1), dim_(i,2))
#ifdef __DEBUG
        do k=0,mpi_global%size_-1
            if (mpi_global%rank==k) then
                write(format_,'("(I8,2I10,5X,A,",I6,"(I0,X))")') dim_(i,1)*dim_(i,2)
                write(*,trim(adjustl(format_))) mpi_global%rank, icontxts(i), &
                                                i, " Map = ", map
            end if
            call mpi_barrier(mpi_global%comm, mpierr)
        end do
        call mpi_barrier(mpi_global%comm, mpierr)
#endif
        deallocate(map)
    end do
    
    call mpi_barrier(mpi_global%comm, mpierr)
    grid%context = icontxts(mpi_local%color)
    grid%nprow = dim_(mpi_local%color,1)
    grid%npcol = dim_(mpi_local%color,2)
    call blacs_gridinfo(grid%context, grid%nprow, grid%npcol, grid%myprow, grid%mypcol)
#ifdef __DEBUG
    call mpi_barrier(mpi_global%comm, mpierr)
    debug_str = 'Localized BLACS grid initiated'
    call debug_output(0)
    debug_str = 'GlobRank    Color   LocRank  Myprow   Nprow  Mypcol  Npcol         Context'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do k=0,mpi_global%size_-1
        if (mpi_global%rank==k) then
            write(*,'(7I8,I16)') mpi_global%rank, mpi_local%color, mpi_local%rank, &
                                 grid%myprow, grid%nprow, grid%mypcol, grid%npcol, &
                                 grid%context
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif
    deallocate(dim_)
    deallocate(icontxts)


#else

    ! No Kpools activated
    !
    ! The processors are divided into a blacs grid of nprow*npcol
    ! Each processor is assigned an unique row index (myprow) and 
    ! an unique column index (mypcol)
    ! --------------------------------

    call blacs_pinfo(grid%rank, grid%size_)
    call get_optimum_grid(mpi_global%size_, grid%npcol, grid%nprow)
    call blacs_get(-1,0,grid%context)
    call blacs_gridinit(grid%context, 'C', grid%nprow, grid%npcol)
    call blacs_gridinfo(grid%context, grid%nprow, grid%npcol, grid%myprow, grid%mypcol) 


#ifdef __DEBUG
    debug_str = ' '
    do i=1,2
        call debug_output(0)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)
    debug_str = 'BLACS grid initiated successfully'
    call debug_output(0)
    debug_str = 'BLACS grid info:'
    call debug_output(0)
    debug_str = 'Rank       Size     myprow     nprow    mypcol    npcol'
    call debug_output(0)
    call mpi_barrier(mpi_global%comm, mpierr)
    do i=0,mpi_global%size_-1
        if (mpi_global%rank==i) then
            write(*,'(I4,5I10 )') grid%rank, grid%size_, grid%myprow, & 
                                  grid%nprow, grid%mypcol, grid%npcol
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif

#endif

    call mpi_barrier(mpi_global%comm, mpierr)
    return

end subroutine








subroutine get_optimum_grid(size_, npcol, nprow)

    implicit none
    integer, intent(in) :: size_
    integer, intent(out) :: npcol, nprow

    do nprow = int(sqrt(real(size_)))+1, 1, -1
        npcol = size_ / nprow
        if (nprow*npcol == size_ ) exit
    end do

    return

end subroutine






#ifdef __KPOOL
subroutine get_color(proc_id, total_procs, no_of_groups, val, key_val, num_each)

    ! returns the color to split the communicators depending
    ! upon the total number of processes and the number of
    ! groups to be created
    !
    ! color = group no. that the proc belongs to
    ! ------------------------------------------------------

    implicit none

    integer, intent(in)  :: proc_id, total_procs, no_of_groups
    integer, intent(out) :: val, key_val, num_each
    integer :: rem_each

    num_each = total_procs/no_of_groups
    rem_each = mod(total_procs,no_of_groups)

    if (proc_id.le.(num_each*(no_of_groups-rem_each)-1)) then
        val = proc_id/num_each + 1
        key_val = 5+val+proc_id
    else
        num_each = num_each + 1
        val = (proc_id+no_of_groups-rem_each)/(num_each) + 1
        key_val = 877+val+proc_id
    end if

    return

end subroutine
#endif
