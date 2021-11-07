SUBROUTINE read_input_file(inp_file)
 
  USE global_variables

  IMPLICIT NONE
  CHARACTER(500), INTENT(IN) :: inp_file
  CHARACTER(500) :: temp

  OPEN(unit=12, file=TRIM(ADJUSTL(inp_file)), action='read')

  READ(12,'(A)') temp
  READ(12,'(A)') temp
  READ(12,'(A)') temp
  READ(12,'(A)') temp
  READ(12,'(A)') temp
  READ(12,'(A)') lammps_file%location
  fc%location = lammps_file%location
  q_file%location = fc%location
  READ(12,'(A)') temp
  READ(12,'(A)') lammps_file%name
  READ(12,'(A)') temp
  READ(12,'(A)') lammps_file%suffix
  READ(12,'(A)') temp
  READ(12,*)     moire%natom
  READ(12,'(A)') temp
  READ(12,*)     lammps_file%at_types
  READ(12,'(A)') temp
  READ(12,'(A)') lammps_file%atom_style
  READ(12,'(A)') temp
  READ(12,'(A)') q_file%name
  READ(12,'(A)') temp
  READ(12,*)     q_file%nqpt
  READ(12,'(A)') temp
  READ(12,*)     no_of_groups
  READ(12,'(A)') temp
  READ(12,*)     pzheevx_vars%mb
  READ(12,'(A)') temp
  READ(12,*)     pzheevx_vars%nb
  READ(12,'(A)') temp
  READ(12,'(A)') fc%file_name
  READ(12,'(A)') temp
  READ(12,'(A)') fc%dataset
  READ(12,'(A)') temp
  READ(12,*)     pzheevx_vars%comp_evec
  READ(12,'(A)') temp
  READ(12,'(A)') pzheevx_vars%range_
  READ(12,'(A)') temp
  READ(12,*)     pzheevx_vars%vl_, pzheevx_vars%vu_
  READ(12,'(A)') temp
  READ(12,*)     pzheevx_vars%il_, pzheevx_vars%iu_

  CALL MPI_BARRIER(mpi_global%comm,mpierr)

  IF ((pzheevx_vars%range_.ne.'A').and. &
      (pzheevx_vars%range_.ne.'V').and. &
      (pzheevx_vars%range_.ne.'I')) THEN
      IF (mpi_global%rank.eq.0) THEN
          WRITE(*,*) "Invalid operation"
          WRITE(*,*) "Entered value for range is illegal."
          WRITE(*,*) "Input 'A' for all eigenvalues to be found"
          WRITE(*,*) "Input 'V' for all eigenvalues to be found in the interval specified"
          WRITE(*,*) "Input 'I' for all eigenvalues to be found within the indices specified."
      END IF
      STOP
  END IF


  CLOSE(unit=12)

  RETURN

END SUBROUTINE 

