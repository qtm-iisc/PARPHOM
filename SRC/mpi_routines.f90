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






!subroutine get_optimum_grid_dist()
!
!    ! gets the optimum distribution of processes in each group 
!    !
!    ! Tries to distribute grids as close to a square as possible
!    ! ----------------------------------------------------------
!
!    use global_variables
!
!    implicit none
!    integer :: nprow, npcol
!    integer :: col, sz, flag
!
!    do nprow = int(sqrt(real(mpi_local%size_)))+1,1,-1
!        npcol = mpi_local%size_/nprow
!        if (nprow*npcol.eq.mpi_local%size_) goto 11
!    end do
!11  continue
!    if (mpi_local%rank==0) then
!        dim_(mpi_local%color,1) = nprow
!        dim_(mpi_local%color,2) = npcol
!    end if
!    return
!end subroutine


subroutine get_blacs_icontexts()

    use mpi
    use global_variables

    ! Subroutine splitting global communicators into groups and 
    ! obtaining blacs context for each group.
    ! ---------------------------------------------------------

    implicit none
    integer, allocatable, dimension(:,:) :: map, dim_
    integer, allocatable, dimension(:)   :: icontxts
    integer :: st_group,i,j,k, global_contxt
    character(500) :: format_

    call get_color(mpi_global%rank, mpi_global%size_, no_of_groups, &
                   mpi_local%color, mpi_local%key, mpi_local%size_)
    
    call mpi_comm_split(mpi_global%comm, mpi_local%color, mpi_local%key, &
                        mpi_local%comm, mpierr)
    call mpi_comm_rank(mpi_local%comm, mpi_local%rank, mpierr)

    allocate(dim_(no_of_groups,2))
    allocate(icontxts(no_of_groups))
    dim_ = 0
    !call get_optimum_grid_dist()

    do i = int(sqrt(real(mpi_local%size_)))+1,1,-1
        blacs_grid%nprow = i
        blacs_grid%npcol = mpi_local%size_/blacs_grid%nprow
        if (blacs_grid%nprow*blacs_grid%npcol.eq.mpi_local%size_) goto 11
    end do
11  continue
    if (mpi_local%rank==0) then
        dim_(mpi_local%color,1) = blacs_grid%nprow
        dim_(mpi_local%color,2) = blacs_grid%npcol
    end if

    call mpi_allreduce(MPI_IN_PLACE,dim_,no_of_groups*2,MPI_INT,MPI_SUM, & 
                       mpi_global%comm, mpierr)
    call blacs_get(-1,0,global_contxt)
    icontxts = global_contxt 
    call mpi_barrier(mpi_global%comm,mpierr)
    st_group = 0
    do i=1,no_of_groups
        allocate(map(dim_(i,1),dim_(i,2)))
        do j=1,dim_(i,1)
            do k=1,dim_(i,2)
                map(j,k) = st_group
                st_group = st_group+1
            end do  
        end do
        call blacs_gridmap(icontxts(i), map, dim_(i,1), dim_(i,1), dim_(i,2))

!        ! Uncomment the following lines for debugging
!        ! -------------------------------------------
!
!        do k=0,mpi_global%nprocs-1
!            if (mpi_global%rank==k) then
!                write(format_,'("(2I4,A,I4,5X,A,",I4,"I4)")') dim_(i,1)*dim_(i,2)
!                write(*,TRIM(ADJUSTL(format_))) mpi_global%rank, icontxts(i), &
!                                                "   Group = ", i, "Map = ", map
!            end if
!            call mpi_barrier(mpi_global%comm,mpierr)
!        end do
        
        deallocate(map)
    end do
    
    !blacs_grid%context = icontxts(mpi_local%color)
    blacs_grid%context = icontxts(mpi_local%color)
    blacs_grid%nprow = dim_(mpi_local%color,1)
    blacs_grid%npcol = dim_(mpi_local%color,2)
    deallocate(icontxts)
    deallocate(dim_)
    return
end subroutine
