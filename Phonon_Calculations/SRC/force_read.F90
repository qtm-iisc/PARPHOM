subroutine read_force_constants()
    
    use hdf5 
    use mpi
    use global_variables

    implicit none

    integer, external :: numroc
    integer :: rsrc, csrc
    character(len=char_len) :: file_name_, done_line
    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend, ia, ja
    integer :: lrindx, lcindx, lroffset, lcoffset, ipos
    integer(4) :: hdf5_error

    integer(hid_t) :: file_id, dset_id
    integer(hid_t) :: dataspace, memspace
    integer(hsize_t), dimension(4,1) :: coord
    integer(hsize_t), dimension(1) :: coord_mem, data_dims
    integer(size_t), parameter :: ONE_ = 1D0

    double precision, parameter :: alpha = 0.500000000000000

#ifdef __DEBUG
    integer :: dataspace_npoints, memspace_npoints, i, j
    character(len=char_len) :: format_
#endif

    rsrc = 0
    csrc = 0

    if (grid%myprow.ge.force_const%desca(RSRC_)) then
        ia_first = (grid%myprow - force_const%desca(RSRC_))*force_const%desca(MB_)+1
    else
        ia_first = (grid%myprow + (grid%nprow - force_const%desca(RSRC_)))* &
                    force_const%desca(MB_) + 1
    endif

    if (grid%mypcol.ge.force_const%desca(CSRC_)) then
        ja_first = (grid%mypcol - force_const%desca(CSRC_))*force_const%desca(NB_)+1
    else
        ja_first = (grid%mypcol + (grid%npcol - force_const%desca(CSRC_)))* &
                    force_const%desca(NB_) + 1
    endif

    write(file_name_, '(2A)') trim(adjustl(force_const%location)), &
                              trim(adjustl(force_const%name_))

    call h5open_f(hdf5_error)
    call h5fopen_f(trim(adjustl(file_name_)),H5F_ACC_RDONLY_F, file_id, hdf5_error)
    if (hdf5_error.ne.0) then
        write(debug_str,'(2A)') "HDF5 could not open the Force constant file", &
                                trim(adjustl(file_name_))
        call debug_output(hdf5_error)
        call exit
    end if
    call h5dopen_f(file_id, trim(adjustl(force_const%dset_name)), dset_id, hdf5_error)
    if (hdf5_error.ne.0) then
        write(debug_str,'(4A)') "HDF5 could not open dataset",        &
                                trim(adjustl(force_const%dset_name)), & 
                                "in the Force constant file", &
                                trim(adjustl(file_name_))
        call debug_output(hdf5_error)
        call exit
    end if
    call h5dget_space_f(dset_id, dataspace, hdf5_error)

    data_dims(1) = force_const%size_
    
    call h5screate_simple_f(1, data_dims, memspace, hdf5_error)

    do jastart = ja_first, force_const%desca(N_), grid%npcol*force_const%desca(NB_)
        do iastart = ia_first, force_const%desca(M_), grid%nprow*force_const%desca(MB_)
            iaend = min(force_const%desca(M_), iastart+force_const%desca(MB_)-1)
            jaend = min(force_const%desca(N_), jastart+force_const%desca(NB_)-1)
            ia = iastart
            ja = jastart
            call infog2l(ia, ja, force_const%desca, grid%nprow, grid%npcol, &
                         grid%myprow, grid%mypcol, lroffset, lcoffset, rsrc, csrc)
            do ja = jastart, jaend
                do ia = iastart, iaend
                    lrindx = lroffset + ia - iastart
                    lcindx = lcoffset + ja - jastart
                    ipos = lrindx + (lcindx-1)*force_const%desca(LLD_)
                    call twod_to_fourd(ia,ja,coord)
                    call h5sselect_elements_f(dataspace, H5S_SELECT_APPEND_F, 4, &
                        ONE_, coord, hdf5_error)
                    coord_mem(1) = ipos
                    call h5sselect_elements_f(memspace, H5S_SELECT_APPEND_F, 1, &
                        ONE_, coord_mem, hdf5_error)
                end do
            end do
        end do
    end do

#ifdef __DEBUG
    call h5sget_select_elem_npoints_f(dataspace, dataspace_npoints, hdf5_error)
    call h5sget_select_elem_npoints_f(memspace, memspace_npoints, hdf5_error)
    do i=0,mpi_global%size_-1
        if (i==mpi_global%rank) then
            write(debug_str,"(3(A,I0),A)") "Processor (",grid%myprow,",",grid%mypcol, &
                                           ") selected ",dataspace_npoints, &
                                           " points in the force constant file."
            write(*,'(A)') trim(adjustl(debug_str))
            write(debug_str,"(3(A,I0),A)") "Processor (",grid%myprow,",",grid%mypcol, &
                                           ") selected ",memspace_npoints,&
                                           " points in the memory space to write data."   
            write(*,'(A)') trim(adjustl(debug_str))
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
#endif

    call h5dread_f(dset_id, H5T_IEEE_F64LE, force_const%mat, data_dims, hdf5_error, &
                   mem_space_id = memspace, file_space_id = dataspace)

    call h5sclose_f(memspace, hdf5_error)
    call h5sclose_f(dataspace, hdf5_error)
    call h5dclose_f(dset_id, hdf5_error)
    call h5fclose_f(file_id, hdf5_error)
    call h5close_f(hdf5_error)

#ifdef __DEBUG
    write(debug_str,'(A)') "The local force constant elements in each rank are: "
    call debug_output(0)
    do i=1,mpi_global%size_
        write(format_, '(A,I0,A)') '(A,I0,A,',force_const%size_,'F12.6)'
        if (mpi_global%rank==i-1) then
            write(*,format_) "Rank :", mpi_global%rank, &
                             "\r\nForce Constant : \r\n", (force_const%mat(j), &
                             j=1,force_const%size_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)
#endif


    if (print_progress) then
        write(done_line, "(A)") "Force constants read on "
        call date_time_message(trim(done_line))
    end if
    
    call pdgeadd('T', 3*moire%natom, 3*moire%natom, &
                  alpha, force_const%mat, 1, 1, force_const%desca, &
                  alpha, force_const%mat, 1, 1, force_const%desca)

#ifdef __DEBUG
    write(debug_str,'(A)') "Local force constants in each rank after symmetrization are: "
    call debug_output(0)
    do i=1,mpi_global%size_
        write(format_, '(A,I0,A)') '(A,I0,A,',force_const%size_,'F12.6)'
        if (mpi_global%rank==i-1) then
            write(*,format_) "Rank :", mpi_global%rank, &
                             "\r\nForce Constant : \r\n", (force_const%mat(j), &
                             j=1,force_const%size_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)

#endif

    if (print_progress) then
        write(done_line, "(A)") "Force constants symmetrized on " 
        call date_time_message(trim(done_line))
    end  if
    return

end subroutine




SUBROUTINE twod_to_fourd (i,j,coord)
  USE HDF5
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, j
  INTEGER(HSIZE_T), DIMENSION(4,1), INTENT(OUT) :: coord

  IF (MOD(i,3).ne.0) THEN
      coord(1,1) = MOD(i,3)
  ELSE
      coord(1,1) = 3
  END IF
  IF (MOD(j,3).ne.0) THEN
      coord(2,1) = MOD(j,3)
  ELSE
      coord(2,1) = 3
  END IF
  IF (MOD(i,3).ne.0) THEN
      coord(3,1) = i/3 + 1
  ELSE
      coord(3,1) = i/3
  END IF
  IF (MOD(j,3).ne.0) THEN
      coord(4,1) = j/3 + 1
  ELSE
      coord(4,1) = j/3
  END IF

  RETURN 

END SUBROUTINE
