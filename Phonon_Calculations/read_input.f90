SUBROUTINE read_input(inp_file, location, fil, suffix, natom, at_types, atoms_style, mb, nb, &
                      num_q_pts, q_file, FC_file, FC_dataset_name)
  
  IMPLICIT NONE
  
  CHARACTER(500), INTENT(OUT)  :: location, fil, suffix, q_file, FC_file, FC_dataset_name
  CHARACTER(1),   INTENT(OUT)  :: atoms_style
  INTEGER,        INTENT(OUT)  :: natom, at_types, mb, nb, num_q_pts

  CHARACTER(500), INTENT(IN) :: inp_file
  INTEGER :: i

  OPEN(unit=12, file=TRIM(ADJUSTL(inp_file)), action='read')

  READ(12,*) location
  READ(12,*) fil
  READ(12,*) suffix
  READ(12,*) natom
  READ(12,*) at_types
  READ(12,*) atoms_style
  READ(12,*) mb
  READ(12,*) nb
  READ(12,*) num_q_pts
  READ(12,*) q_file
  READ(12,*) FC_file
  READ(12,*) FC_dataset_name

  CLOSE(unit=12)

  RETURN
END SUBROUTINE 

!PROGRAM TEST
!  
!  IMPLICIT NONE
!  
!  CHARACTER(500) :: location, fil, suffix, q_file, FC_file, FC_dataset_name, inp_file
!  CHARACTER(1)   :: atoms_style
!  INTEGER        :: natom, at_types, mb, nb, num_q_pts
!
!  CALL GET_COMMAND_ARGUMENT(1,inp_file)
!  IF (LEN_TRIM(inp_file) == 0) STOP "Enter proper input file"
!
!  CALL read_input(inp_file, location, fil, suffix, natom, at_types, atoms_style, &
!                  mb, nb, num_q_pts, q_file, FC_file, FC_dataset_name) 
!  WRITE(*,*) "LOCATION= ", TRIM(ADJUSTL(location)) 
!  WRITE(*,*) TRIM(ADJUSTL(fil)) 
!  WRITE(*,*) TRIM(ADJUSTL(suffix)) 
!  WRITE(*,*) TRIM(ADJUSTL(FC_file)) 
!  WRITE(*,*) TRIM(ADJUSTL(FC_dataset_name)) 
!  WRITE(*,*) TRIM(ADJUSTL(q_file)) 
!  WRITE(*,*) atoms_style
!  WRITE(*,*) mb, nb, num_q_pts, at_types, natom  
!
!END PROGRAM


   

