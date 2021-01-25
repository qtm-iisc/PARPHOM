SUBROUTINE twod_to_fourd (i,j,coord)
  USE HDF5
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i, j
  INTEGER(HSIZE_T), DIMENSION(4,2), INTENT(OUT) :: coord

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
   

  IF (MOD(j,3).ne.0) THEN
      coord(1,2) = MOD(j,3)
  ELSE
      coord(1,2) = 3
  END IF
  IF (MOD(i,3).ne.0) THEN
      coord(2,2) = MOD(i,3)
  ELSE
      coord(2,2) = 3
  END IF
  IF (MOD(j,3).ne.0) THEN
      coord(3,2) = j/3 + 1
  ELSE
      coord(3,2) = j/3
  END IF
  IF (MOD(i,3).ne.0) THEN
      coord(4,2) = i/3 + 1
  ELSE
      coord(4,2) = i/3
  END IF


  RETURN
END SUBROUTINE


