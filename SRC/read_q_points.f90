SUBROUTINE read_q_file()
  
  USE global_variables

  IMPLICIT NONE
  INTEGER :: i
  
  OPEN(unit=1, file=TRIM(ADJUSTL(q_file%location))//TRIM(ADJUSTL(q_file%name)), action="read")

  if (st_q==1) then
    DO i=st_q,en_q
        READ(1,*) q_file%points(i,1), q_file%points(i,2), q_file%points(i,3)
    END DO
  else
    DO i=1,st_q-1
        READ(1,*)
    END DO
    DO i=1,en_q-st_q+1
        READ(1,*) q_file%points(i,1), q_file%points(i,2), q_file%points(i,3)
    END DO
  end if

  CLOSE(1)

  RETURN

END SUBROUTINE 
