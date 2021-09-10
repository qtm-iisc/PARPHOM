subroutine write_output()

    use mpi
    use hdf5
    use global_variables

    implicit none
    integer :: i
    character(500) :: f_name

    if (mpi_local%rank==0) then
        write(f_name,43) mpi_local%color-1, q_index-1
        open(unit=20, file=trim(adjustl(f_name)),action='write')
        write(20,'(3(F8.5X))') q_file%points(q_index,:)
        do i=1,pzheevx_vars%num_eval_comp
            write(20,'(A,I6,A,F16.6)') 'W(',i,') =',sqrt(abs(W(i)))*15.633302* &
                                                   sign(1.0,W(i))*33.3564
        end do
        close(20)
    end if

43  FORMAT('P_vec_ele_',i4.4,'_',i3.3)

end subroutine
