subroutine read_force_constants()
    
    use mpi
    use hdf5
    use global_variables

    implicit none

    integer :: hdf5err
    integer(HID_T) :: file_id, dset_id, space, memspace
    integer(SIZE_T), parameter :: num_ele = 2
    integer(HSIZE_T), dimension(4, num_ele) :: coord
    integer(HSIZE_T), dimension(1) :: count_, data_dims
    double precision, dimension(num_ele) :: phi_ij
    integer :: ia_first, ja_first, jastart, iastart, iaend, jaend, ia, ja
    integer :: lroffset, lcoffset, rsrc, csrc, lrindx, lcindx, ipos
    double precision :: st_time, en_time

    rsrc = 0
    csrc = 0

    st_time = MPI_Wtime()

    count_(1) = num_ele
    data_dims(1) = num_ele

    call h5open_f(hdf5err)
    call h5fopen_f(TRIM(ADJUSTL(fc%location))//TRIM(ADJUSTL(fc%file_name)), &
                   H5F_ACC_RDONLY_F, file_id, hdf5err)
    call h5dopen_f(file_id,TRIM(ADJUSTL(fc%dataset)), dset_id, hdf5err)
    call h5dget_space_f(dset_id, space, hdf5err)

    call first_array_index(fc%desca, ia_first, ja_first)

    do jastart = ja_first, fc%desca(N_), fc%desca(NB_)*blacs_grid%npcol
        do iastart = ia_first, fc%desca(M_), fc%desca(MB_)*blacs_grid%nprow
            iaend = min(fc%desca(M_), iastart+fc%desca(MB_)-1)
            jaend = min(fc%desca(N_), iastart+fc%desca(NB_)-1)

            call infog2l(iastart,jastart,fc%desca,blacs_grid%nprow, blacs_grid%npcol, &
                         blacs_grid%myprow, blacs_grid%mypcol, lroffset, lcoffset,    &
                         rsrc, csrc)
            do ja = jastart, jaend
                do ia = iastart, iaend
                    lrindx = lroffset+(ia-iastart)
                    lcindx = lcoffset+(ja-jastart)
                    ipos = lrindx + (lcindx-1)*fc%desca(LLD_)

                    call twod_to_fourd(ia,ja,coord)
                    call h5sselect_elements_f(space, H5S_SELECT_SET_F,4,    &
                                              num_ele, coord, hdf5err)
                    call h5screate_simple_f(1,count_,memspace,hdf5err)
                    call h5dread_f(dset_id,H5T_IEEE_F64LE,phi_ij,data_dims, &
                                   hdf5err, memspace, file_space_id = space)
                    fc%a(ipos) = (phi_ij(1)+phi_ij(2))*0.5 
                    call h5sclose_f(memspace,hdf5err)
                end do
            end do
        end do
    end do

    call h5sclose_f(space,hdf5err)
    call h5dclose_f(dset_id,hdf5err)
    call h5fclose_f(file_id,hdf5err)
    call h5close_f(hdf5err)

    en_time = MPI_Wtime()

    if ((blacs_grid%myprow==0).and.(blacs_grid%mypcol==0)) then
        write(*,*) "Force constant read in group ", mpi_local%color," in ", en_time-st_time, " sec."
    end if

    return
end subroutine
