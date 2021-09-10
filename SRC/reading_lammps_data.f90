SUBROUTINE read_lammps_data()

  USE global_variables

  IMPLICIT NONE
  INTEGER :: i, temp
  DOUBLE PRECISION :: x1, x2, y1, y2, z1, z2, xy, xz, yz 
  
  OPEN(unit=1,file=TRIM(ADJUSTL(lammps_file%location))//TRIM(ADJUSTL(lammps_file%name))//TRIM(ADJUSTL(lammps_file%suffix)),action="read")
  
  ! Read lattice parameters

  DO i=1,5
    READ(1,*)
  END DO

  READ(1,*) x1,x2
  READ(1,*) y1,y2
  READ(1,*) z1,z2
  READ(1,*) xy, xz, yz
  
  moire%lat(1,1) = x2-x1
  moire%lat(1,2) = 0
  moire%lat(1,3) = 0
  moire%lat(2,1) = xy
  moire%lat(2,2) = y2-y1
  moire%lat(2,3) = 0
  moire%lat(3,1) = 0
  moire%lat(3,2) = 0
  moire%lat(3,3) = z2-z1
  
  ! Read masses

  DO i=1,3
    READ(1,*)
  END DO
  
  DO i=1,lammps_file%at_types
    READ(1,*) temp, moire%mass(i)
  END DO

  ! Read real coordinates

  DO i=1,3
    READ(1,*)
  END DO 

  IF (lammps_file%atom_style.eq."A") THEN
    DO i=1,moire%natom
      READ(1,*) temp, moire%at_types_i(i), moire%real_pos(i,1), &
                      moire%real_pos(i,2), moire%real_pos(i,3)
    END DO
  ELSE IF (lammps_file%atom_style.eq."M") THEN
    DO i=1,moire%natom
      READ(1,*) temp, temp, moire%at_types_i(i), moire%real_pos(i,1), &
                            moire%real_pos(i,2), moire%real_pos(i,3)
    END DO
  ELSE
    IF (mpi_global%rank==0) THEN
        WRITE(*,*) "Atom style not supported"
    END IF
    STOP
  END IF

  ! Get crys

  DO i=1,moire%natom
    CALL linsolve(moire%real_pos(i,:),moire%lat,moire%crys(i,:))
  END DO
  
  CALL matinvtransp3(moire%lat, moire%rec_lat)
  
  RETURN
END SUBROUTINE


SUBROUTINE matinvtransp3(A, B)
        !! Performs a direct calculation of the inverse of a 3×3 matrix.
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: A(3,3)   !! Matrix
  DOUBLE PRECISION, INTENT(OUT):: B(3,3)   !! Inverse matrix
  DOUBLE PRECISION             :: detinv
  DOUBLE PRECISION, PARAMETER  :: PI2 = 6.283185307179586

  ! Calculate the inverse determinant of the matrix
  
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix

  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))*PI2
  B(1,2) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))*PI2
  B(1,3) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))*PI2
  B(2,1) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))*PI2
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))*PI2
  B(2,3) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))*PI2
  B(3,1) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))*PI2
  B(3,2) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))*PI2
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))*PI2
       
  RETURN
END SUBROUTINE

SUBROUTINE linsolve(x, A, u)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: x
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: A
  DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: u
  DOUBLE PRECISION :: detinv
  DOUBLE PRECISION, DIMENSION(3,3) :: B

  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
  
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  
  u = MATMUL(TRANSPOSE(B),x)
  
  RETURN
END SUBROUTINE
