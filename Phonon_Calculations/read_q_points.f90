SUBROUTINE read_q_file(q_file, location, nqpt, q_array)
  IMPLICIT NONE
  CHARACTER(500), INTENT(IN) :: q_file, location
  INTEGER,        INTENT(IN) :: nqpt
  DOUBLE PRECISION, DIMENSION(nqpt,3) :: q_array
  INTEGER :: i
  
  OPEN(unit=1, file=TRIM(ADJUSTL(location))//TRIM(ADJUSTL(q_file)), action="read")

  DO i=1,nqpt
    READ(1,*) q_array(i,1), q_array(i,2), q_array(i,3)
  END DO

  RETURN

END SUBROUTINE 

!PROGRAM TEST
!  IMPLICIT NONE
!  CHARACTER(500) :: q_file, location
!  INTEGER, PARAMETER :: nqpt = 101
!  DOUBLE PRECISION, ALLOCATABLE :: q_array(:,:)
!  INTEGER :: i  
!
! ALLOCATE(q_array(nqpt,3))
!
!  q_file = "q_points.dat"
!  location = "/mnt/lustre/phy3/physhinj/Parallel_Phonon/FINAL_VERSION_2/"
! 
!  CALL read_q_file(q_file, location, nqpt, q_array)
!
!  DO i=1,nqpt
!    WRITE(*,"(3F16.6)") q_array(i,1), q_array(i,2), q_array(i,3)
!  END DO
!  
!  DEALLOCATE(q_array)
!
!END PROGRAM 
  
