! Package: PARPHOM
! Authors: Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain
! License: GPL-3.0
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!> \brief Write phonon calculation results (eigenvalues, eigenvectors, velocities) to an HDF5 file in parallel.
!>
!> This subroutine manages all output operations for a single q-point in a parallel phonon calculation.
!> It performs the following major tasks:
!>   - Checks if the HDF5 output file exists. If not, creates a new file. If it exists, opens it for parallel access.
!>   - Creates a group for the current q-point and writes the q-vector to the file.
!>   - Stores the computed phonon eigenvalues for the q-point, converting them to physical units (cm^-1).
!>   - If eigenvectors are computed, redistributes the distributed matrix blocks and writes the complex eigenvectors
!>     to the file using a compound HDF5 datatype for real and imaginary parts.
!>   - If group velocities are computed, constructs the velocity matrix for each Cartesian direction by evaluating
!>     derivatives of the dynamical matrix, and writes the diagonal elements (group velocities) to the file.
!>   - All I/O is performed in parallel using MPI and HDF5 collective operations for efficiency and scalability.
!>   - Handles all HDF5 resource management, including creation and closing of property lists, dataspaces, datasets,
!>     and groups, as well as MPI communicator splitting for selective writing.
!>   - Provides progress messages if enabled.
!>
!> \param[in] q_indx Integer index of the q-point for which results are written.
!>
!> \details
!>   - The subroutine is designed for distributed-memory parallel environments (MPI).
!>   - It assumes the presence of global variables and distributed arrays for eigenvalues, eigenvectors, and velocities.
!>   - The output file structure is hierarchical: each q-point is a group containing datasets for q-vector, eigenvalues,
!>     eigenvectors, and velocities.
!>   - Eigenvectors are stored as compound datasets with separate real and imaginary fields.
!>   - Velocities are written for each Cartesian direction and band, normalized by the phonon frequency.
!>   - All error handling is performed via HDF5 and MPI error codes.
!>
!> \authors Shinjan Mandal, Indrajit Maity, H R Krishnamurthy, Manish Jain


subroutine write_output(q_indx)

    use hdf5
    use global_variables
    use mpi

    implicit none
    
    integer, intent(in) :: q_indx
    character(len=char_len) :: file_name, group_name, done_line
    logical :: file_exists
    
    ! HDF5 Variables
    ! --------------
    integer(4) :: hdf5_error
    integer(hid_t) :: plist_id, glist_id, dlist_id
    integer(hid_t) :: file_id
    integer(hid_t) :: dset_id !, dset_id2
    integer(hid_t) :: group_id
    integer(hid_t) :: filespace, memspace, dataspace_id !, dataspace_id2
    integer(hsize_t), dimension(1) :: dim_eval, dim_q, dim_mem
    integer(hsize_t), dimension(2) :: dim_evec
    integer(hsize_t), dimension(2) :: dim_vel

    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend, ia, ja
    integer :: lrindx, lcindx, lroffset, lcoffset, rsrc, csrc
    integer(hsize_t), dimension(1) :: ipos
    integer(hsize_t), dimension(2,1) :: coord
    integer(hsize_t), dimension(3,1) :: coord_vel
    integer(size_t), parameter :: ONE_ = 1

    integer(size_t) :: type_size, type_size_r, type_size_i, offset
    integer(hid_t) :: dtr_id, dti_id, dtype_id

    integer :: proc_should_write, i
    integer :: write_comm
   
    integer, external :: numroc

    double precision, allocatable, dimension(:) :: temp
    double complex, parameter :: alpha = cmplx(1,0)
    double complex, parameter :: beta = cmplx(0,0)

    ! Step 1: Open HDF5 interface
    call h5open_f(hdf5_error)

    ! Step 2: Prepare output file name and check if it exists
    write(file_name, '(2A)') trim(adjustl(output_file_location)), &
                             trim(adjustl(output_file_name))
    inquire(file=file_name, exist=file_exists)
    write(group_name,'(I12.12)') q_indx 

    ! Step 3: Create or open the HDF5 file (serial, root process only)
    if (file_exists) then
        if (mpi_global%rank == 0) then
            call h5fopen_f(trim(adjustl(file_name)), H5F_ACC_RDWR_F, file_id, hdf5_error)
        end if
    else
        if (mpi_global%rank==0) then
            call h5fcreate_f(trim(adjustl(file_name)), H5F_ACC_TRUNC_F, file_id, hdf5_error)
        end if
    end if

    ! Step 4: Create group for this q-point and write q-vector (root only)
    if (mpi_global%rank==0) then
        call h5gcreate_f(file_id, group_name, group_id, hdf5_error)
        dim_q(1) = 3
        call h5screate_simple_f(1, dim_q, filespace, hdf5_error)
        call h5dcreate_f(group_id, 'q_vec', H5T_IEEE_F64LE, filespace, dset_id, hdf5_error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, q_file%points(q_indx,:),dim_q,hdf5_error)
        call h5dclose_f(dset_id,hdf5_error)
        call h5sclose_f(filespace, hdf5_error)
        call h5gclose_f(group_id,hdf5_error)
        call h5fclose_f(file_id, hdf5_error)
    end if

    ! Step 5: Open file and group in parallel (all processes)
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_error)
    call h5pset_fapl_mpio_f(plist_id, int(mpi_global%comm,kind=4) , MPI_INFO_NULL, hdf5_error)
    call h5fopen_f(trim(adjustl(file_name)), H5F_ACC_RDWR_F, file_id, hdf5_error, &
                   access_prp = plist_id)    
    call h5pcreate_f(H5P_GROUP_ACCESS_F, glist_id, hdf5_error)
    call h5pset_all_coll_metadata_ops_f(glist_id, .true. , hdf5_error)
    call h5gopen_f(file_id, group_name, group_id, hdf5_error, gapl_id = glist_id)

    ! Step 6: Store eigenvalues (all processes)
    dim_eval(1) = pzheevx_vars%comp_num_eval
    eval = eval*15.633302*33.35641   ! Convert to cm^-1
    call h5screate_simple_f(1, dim_eval, filespace, hdf5_error)
    call h5dcreate_f(group_id, 'eigenvalues', H5T_IEEE_F64LE, filespace, dset_id, hdf5_error)
    call h5pcreate_f(H5P_DATASET_XFER_F, dlist_id, hdf5_error)
    call h5pset_dxpl_mpio_f(dlist_id, H5FD_MPIO_COLLECTIVE_F, hdf5_error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, eval, dim_eval, hdf5_error, xfer_prp = dlist_id)
    call h5sclose_f(filespace,hdf5_error)
    call h5dclose_f(dset_id, hdf5_error)
    call h5pclose_f(dlist_id, hdf5_error)

    ! Step 7: Close group, file, and property lists after eigenvalue write
    call h5gclose_f(group_id,hdf5_error)
    call h5pclose_f(glist_id, hdf5_error)
    call h5pclose_f(plist_id, hdf5_error)
    call h5fclose_f(file_id, hdf5_error)
    
    ! Step 8: Print progress if enabled
    if (print_progress) then
        write(done_line,"(I0,3A)") pzheevx_vars%comp_num_eval, " eigenvalues written to ", &
                                  trim(adjustl(file_name)), " on "
        call date_time_message(trim(done_line))
    end if

    ! Step 9: Store eigenvectors if requested
    ! - Split communicators so only relevant processes write
    ! - Open file/group in parallel for writing
    ! - Create compound datatype for complex numbers
    ! - Redistribute and select elements for writing
    ! - Write real and imaginary parts separately
    ! - Close all HDF5 resources
    ! - Print progress if enabled

    ! Since the local arrays contain elements distributed in block-cyclic fashion
    ! we have to reverse this distribution and select the coordinates in the 
    ! global file that correspond to the elements in the distributed array 
    ! 
    ! If j eigenvectors are computed, only the first j columns of the distributed
    ! matrix store the eigenvectors.
    ! The processes containing the relevant columns writes them out to the file.


    !if (pzheevx_vars%comp_evec=='V') then
    
    if (evec_comp) then

        rsrc = 0
        csrc = 0

        dim_evec(1) = 3*moire%natom
        dim_evec(2) = pzheevx_vars%comp_num_evec

        if (grid%myprow.ge.evec%desca(RSRC_)) then
            ia_first = (grid%myprow - evec%desca(RSRC_))*evec%desca(MB_)+1
        else
            ia_first = (grid%myprow + (grid%nprow - evec%desca(RSRC_)))* &
                        evec%desca(MB_) +1
        endif
        if (grid%mypcol.ge.evec%desca(CSRC_)) then
            ja_first = (grid%mypcol - evec%desca(CSRC_))*evec%desca(NB_)+1
        else
            ja_first = (grid%mypcol + (grid%npcol - evec%desca(CSRC_)))* &
                      evec%desca(NB_) +1
        endif

        proc_should_write = 0

        if (ja_first.le.pzheevx_vars%comp_num_evec) then
            proc_should_write = 1
        end if
        
        call mpi_comm_split(mpi_global%comm, proc_should_write, mpi_global%rank, &
                            write_comm, mpierr)


        if (proc_should_write==1) then

            ! Open file in parallel
            
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_error)
            call h5pset_fapl_mpio_f(plist_id, int(write_comm,kind=4) , MPI_INFO_NULL, &
                                    hdf5_error)
            call h5fopen_f(trim(adjustl(file_name)), H5F_ACC_RDWR_F, file_id, hdf5_error, &
                   access_prp = plist_id)

            ! Open group in parallel
              
            call h5pcreate_f(H5P_GROUP_ACCESS_F, glist_id, hdf5_error)
            call h5pset_all_coll_metadata_ops_f(glist_id, .true. , hdf5_error)
            call h5gopen_f(file_id, group_name, group_id, hdf5_error, gapl_id = glist_id)

            ! Create compoud datatype for storing complex values

            call h5tget_size_f(H5T_IEEE_F64LE, type_size_r, hdf5_error)
            call h5tget_size_f(H5T_IEEE_F64LE, type_size_i, hdf5_error)
            type_size = type_size_r+type_size_i

            call h5tcreate_f(H5T_COMPOUND_F, type_size, dtype_id, hdf5_error)
            offset = 0
            call h5tinsert_f(dtype_id, 'r', offset, H5T_IEEE_F64LE, hdf5_error)
            offset = offset+type_size_r
            call h5tinsert_f(dtype_id, 'i', offset, H5T_IEEE_F64LE, hdf5_error)

            ! Create dataspace

            dim_mem(1) = evec%size_
            
            call h5screate_simple_f(2, dim_evec, dataspace_id, hdf5_error)
            
            ! Create dataset

            call h5dcreate_f(group_id,'evec', dtype_id,dataspace_id,dset_id,hdf5_error)
            call h5pcreate_f(H5P_DATASET_XFER_F, dlist_id, hdf5_error)
            call h5pset_dxpl_mpio_f(dlist_id, H5FD_MPIO_COLLECTIVE_F, hdf5_error)

            call h5screate_simple_f(1, dim_mem, memspace, hdf5_error)
            
            do jastart = ja_first, pzheevx_vars%comp_num_evec, grid%npcol*evec%desca(NB_)
                do iastart = ia_first, evec%desca(M_), grid%nprow*evec%desca(MB_)
                    iaend = min(evec%desca(M_), iastart+evec%desca(MB_)-1)
                    jaend = min(evec%desca(N_), jastart+evec%desca(NB_)-1)
                    ia = iastart
                    ja = jastart
                    call infog2l (ia, ja, evec%desca, grid%nprow, grid%npcol, &
                                  grid%myprow, grid%mypcol, lroffset, lcoffset, &
                                  rsrc, csrc)
                    do ja=jastart,jaend
                        do ia=iastart,iaend
                            lrindx = lroffset + (ia-iastart)
                            lcindx = lcoffset + (ja-jastart)
                            ipos(1) = lrindx + (lcindx-1)*evec%desca(LLD_)
                            
                            coord(1,1) = ia
                            coord(2,1) = ja
                            
                            if (ja.le.pzheevx_vars%comp_num_evec) then
                             call h5sselect_elements_f(memspace, H5S_SELECT_APPEND_F, 1, &
                                    ONE_, ipos, hdf5_error) 
                             call h5sselect_elements_f(dataspace_id,H5S_SELECT_APPEND_F,2,&
                                    ONE_, coord, hdf5_error)  
                            end if
                        end do
                    end do
                end do
            end do
            
            call h5tcreate_f(H5T_COMPOUND_F, type_size_r, dtr_id, hdf5_error)
            offset = 0
            call h5tinsert_f(dtr_id, 'r', offset, H5T_IEEE_F64LE, hdf5_error)
            call h5tcreate_f(H5T_COMPOUND_F, type_size_i, dti_id, hdf5_error)
            offset = 0
            call h5tinsert_f(dti_id, 'i', offset, H5T_IEEE_F64LE, hdf5_error) 
            
            allocate(temp(evec%size_))
            temp = real(evec%mat)
            call h5dwrite_f(dset_id, dtr_id, temp, &
                            dim_evec, hdf5_error, mem_space_id = &
                            memspace, file_space_id = dataspace_id, xfer_prp=dlist_id) 
            temp = aimag(evec%mat)
            call h5dwrite_f(dset_id, dti_id, temp, &
                            dim_evec, hdf5_error, mem_space_id = &
                            memspace, file_space_id = dataspace_id, xfer_prp=dlist_id)  
            deallocate(temp)
            call h5pclose_f(dlist_id, hdf5_error)
            call h5sclose_f(memspace,hdf5_error)
            call h5dclose_f(dset_id, hdf5_error)
            call h5sclose_f(dataspace_id, hdf5_error)
            call h5gclose_f(group_id,hdf5_error)
            call h5pclose_f(glist_id, hdf5_error)
            call h5pclose_f(plist_id, hdf5_error)
            call h5fclose_f(file_id, hdf5_error)
        end if
        call mpi_comm_free(write_comm, mpierr)
        if (print_progress) then
            write(done_line,"(I0,3A)") pzheevx_vars%comp_num_evec," eigenvectors written to ",&
                                          trim(adjustl(file_name)), " on "
            call date_time_message(trim(done_line))
        end if
    end if
    


    ! ------------------------------------------------------------------------------
    !
    !   Store velocities
    !   ----------------
    !   
    !   To store the velocities, a dataspace of dimension (natom,natom) is created
    !   Then the derivatives of the hamiltonians are constructed and the velocity matrix 
    !   is computed as
    !
    !      v_{\alpha,q}^{m,n} = 1/(2\omega_q) < \psi_{m,q} | dD/d\alpha | \psi_{n,q} > 
    !
    !   This generates a distributed velocity matrix of which (num_evec,num_evec)
    !   elements are populated.
    !
    !   The velocities of the different bands are the diagonal elements of this matrix.
    !
    ! --------------------------------------------------------------------------------

    if (comp_vel) then

        rsrc = 0
        csrc = 0


        if (grid%myprow.ge.vel%desca(RSRC_)) then
            ia_first = (grid%myprow - vel%desca(RSRC_))*vel%desca(MB_)+1
        else
            ia_first = (grid%myprow + (grid%nprow - evec%desca(RSRC_)))* &
                        evec%desca(MB_) +1
        endif
        if (grid%mypcol.ge.vel%desca(CSRC_)) then
            ja_first = (grid%mypcol - vel%desca(CSRC_))*vel%desca(NB_)+1
        else
            ja_first = (grid%mypcol + (grid%npcol - vel%desca(CSRC_)))* &
                      vel%desca(NB_) +1
        endif

        proc_should_write = 0

        do jastart = ja_first, pzheevx_vars%comp_num_evec, grid%npcol*vel%desca(NB_)
            do iastart = ia_first, pzheevx_vars%comp_num_evec, grid%nprow*vel%desca(MB_)
                iaend = min(vel%desca(M_), iastart+vel%desca(MB_)-1)
                jaend = min(vel%desca(N_), jastart+vel%desca(NB_)-1)
                do ja=jastart,jaend
                    do ia=iastart,iaend
                        if (ia==ja) then
                            proc_should_write=1
                        end if
                    end do
                end do
            end do
        end do

        ! Splitting the processes into groups which write and do not write


        call mpi_comm_split(mpi_global%comm, proc_should_write, mpi_global%rank, &
                            write_comm, mpierr)
        
        
                        
        if (proc_should_write==1) then

            ! Create dataset in the group that will write to the file

            ! Open file in parallel
            
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, hdf5_error)
            call h5pset_fapl_mpio_f(plist_id, int(write_comm,kind=4) , MPI_INFO_NULL, &
                                 hdf5_error)
            call h5fopen_f(trim(adjustl(file_name)), H5F_ACC_RDWR_F, file_id, hdf5_error, &
                 access_prp = plist_id)

            ! Open group in parallel

            call h5pcreate_f(H5P_GROUP_ACCESS_F, glist_id, hdf5_error)
            call h5pset_all_coll_metadata_ops_f(glist_id, .true. , hdf5_error)
            call h5gopen_f(file_id, group_name, group_id, hdf5_error, gapl_id = glist_id)

            ! Create dataset

            dim_vel(2) = pzheevx_vars%comp_num_evec
            dim_vel(1) = 3
            
            call h5screate_simple_f(2, dim_vel, dataspace_id, hdf5_error)
            call h5dcreate_f(group_id,'vel',H5T_IEEE_F64LE,dataspace_id,dset_id,hdf5_error)
            call h5pcreate_f(H5P_DATASET_XFER_F, dlist_id, hdf5_error)
            call h5pset_dxpl_mpio_f(dlist_id, H5FD_MPIO_COLLECTIVE_F, hdf5_error)
            dim_mem(1) = vel%size_
            call h5screate_simple_f(1, dim_mem, memspace, hdf5_error)
        end if


        
        do i=1,3
          
          ! Compute the velocity for each cartesian direction

          call create_dynamical_matrix(q_indx,i)
          call pzgemm('N','N',3*moire%natom, 3*moire%natom, 3*moire%natom, alpha, &
                       dyn_mat%mat, 1, 1, dyn_mat%desca,            &
                       evec%mat,1,1,evec%desca,beta, vel%mat,1,1,vel%desca)

          call pzgemm('C','N',3*moire%natom, 3*moire%natom, 3*moire%natom, alpha,  &
                       evec%mat, 1, 1, evec%desca,vel%mat,1,1,vel%desca,beta,&
                       vel%mat, 1, 1, vel%desca) 
            
          

          if (proc_should_write==1) then
              do jastart = ja_first, pzheevx_vars%comp_num_evec, grid%npcol*vel%desca(NB_)
                  do iastart = ia_first, pzheevx_vars%comp_num_evec, grid%nprow*vel%desca(MB_)
                      iaend = min(vel%desca(M_), iastart+vel%desca(MB_)-1)
                      jaend = min(vel%desca(N_), jastart+vel%desca(NB_)-1)
                      ia = iastart
                      ja = jastart
                      call infog2l (ia, ja, vel%desca, grid%nprow, grid%npcol, &
                                    grid%myprow, grid%mypcol, lroffset, lcoffset, &
                                    rsrc, csrc)
                      do ja=jastart,jaend
                          do ia=iastart,iaend
                              if ((ia.eq.ja).and.(ia.le.pzheevx_vars%comp_num_evec)) then
                                  lrindx = lroffset + (ia-iastart)
                                  lcindx = lcoffset + (ja-jastart)
                                  ipos(1) = lrindx + (lcindx-1)*vel%desca(LLD_)
                                  coord_vel(2,1) = ja
                                  coord_vel(1,1) = i
                                  vel%mat(ipos) = vel%mat(ipos)/ &
                                                max(2*abs(eval(ja)), 1d-3) &
                                                *33.35641 * (15.633302)**2 * 10**2 &
                                                *2*3.141592653589793 
                                  call h5sselect_elements_f(memspace,H5S_SELECT_APPEND_F,&
                                         1,ONE_, ipos, hdf5_error)
                                  call h5sselect_elements_f(dataspace_id, &
                                         H5S_SELECT_APPEND_F,2,ONE_,coord_vel,hdf5_error)
                              end if
                          end do
                      end do
                  end do
              end do

              allocate(temp(vel%size_))
              temp = real(vel%mat)
              call h5dwrite_f(dset_id, H5T_IEEE_F64LE, temp, &
                              dim_vel, hdf5_error, mem_space_id = &
                              memspace, file_space_id = dataspace_id, xfer_prp=dlist_id)
              deallocate(temp)
              call h5sselect_none_f(memspace, hdf5_error)
              call h5sselect_none_f(dataspace_id,hdf5_error)
          end if
        end do
    
        if (proc_should_write==1) then
            call h5sclose_f(memspace,hdf5_error)
            call h5pclose_f(dlist_id, hdf5_error)
            call h5dclose_f(dset_id, hdf5_error)
            call h5sclose_f(dataspace_id, hdf5_error)
            call h5gclose_f(group_id,hdf5_error)
            call h5pclose_f(glist_id, hdf5_error)
            call h5pclose_f(plist_id, hdf5_error)
            call h5fclose_f(file_id, hdf5_error)
        end if
        call mpi_comm_free(write_comm, mpierr)
        if (print_progress) then 
            write(done_line,"(3A)") "Group velocities written to ",trim(adjustl(file_name)),&
                                    " on "
            call date_time_message(trim(done_line))
        end if
    end if

    ! Step 11: Finalize and close HDF5 interface
    call mpi_barrier(mpi_global%comm, mpierr)
    call h5close_f(hdf5_error)

    return

end subroutine
