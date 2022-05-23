SUBROUTINE read_lammps_data()

  USE global_variables

  IMPLICIT NONE
  INTEGER :: i, temp, error !,j
  DOUBLE PRECISION :: x1, x2, y1, y2, z1, z2, xy, xz, yz, temp_d
  CHARACTER(CHAR_LEN) :: file_name_

  ALLOCATE(moire%mass(lammps_file%at_types))
  ALLOCATE(moire%at_types_i(moire%natom))
  !ALLOCATE(moire%lay_types(moire%natom))
  ALLOCATE(moire%real_pos(moire%natom,3))
  ALLOCATE(moire%crys(moire%natom,3))
  ALLOCATE(moire%sup_cell_info(moire%natom,3))
  
  WRITE(file_name_,'(2A)') TRIM(ADJUSTL(lammps_file%location)), &
                           TRIM(ADJUSTL(lammps_file%name_))
  OPEN(unit=1,file=TRIM(ADJUSTL(file_name_)),action="read",iostat=error)

  IF (error.ne.0) THEN
     WRITE(err_msg,'(2A)') 'Error reading LAMMPS file ', trim(adjustl(file_name_))
     CALL error_message()
     CALL EXIT
  END IF

#ifdef __DEBUG
  WRITE (debug_str,'(A)') '\r\nStarting to read LAMMPS Data'
  CALL debug_output(0)
#endif

  ! -----------------------
  ! Read lattice parameters
  ! -----------------------

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
  moire%lat(3,1) = xz
  moire%lat(3,2) = yz
  moire%lat(3,3) = z2-z1

#ifdef __DEBUG
  WRITE (debug_str,'(A)') '\r\nLattice vectors: '
  CALL debug_output(0)
  DO i=1,3
    WRITE (debug_str, '(3F12.4)') moire%lat(i,1), moire%lat(i,2), moire%lat(i,3)
    CALL debug_output(0)
  END DO
#endif
  
  ! -----------
  ! Read masses
  ! -----------

  DO i=1,3
    READ(1,*)
  END DO
 
  DO i=1,lammps_file%at_types
    READ(1,*) temp, moire%mass(i)
  END DO

#ifdef __DEBUG
  WRITE(debug_str, '(A)') '\r\nMasses'
  CALL debug_output(0)
  DO temp=1,lammps_file%at_types
    WRITE(debug_str, '(A,I0,A,F10.6)') "Atom type ", temp, ": ", moire%mass(temp)
    CALL debug_output(0)
  END DO
#endif

  ! ---------------------
  ! Read real coordinates
  ! ---------------------

  DO i=1,3
    READ(1,*)
  END DO 


  SELECT CASE (lammps_file%atom_style)
    CASE ("atomic","Atomic","ATOMIC")
      DO i=1,moire%natom
        READ(1,*) temp, moire%at_types_i(i), moire%real_pos(i,1), &
                        moire%real_pos(i,2), moire%real_pos(i,3), &
                        moire%sup_cell_info(i,1), moire%sup_cell_info(i,2), &
                        moire%sup_cell_info(i,3)
      END DO
#ifdef __DEBUG
      WRITE(debug_str, '(A)') '\r\n Atom Number      Atom type         X          Y          Z'
      CALL debug_output(0)
      DO i = 1,moire%natom
        WRITE(debug_str, '(2I8, 3F12.6, 3I6)') i, moire%at_types_i(i), moire%real_pos(i,1), &
                                          moire%real_pos(i,2), moire%real_pos(i,3), &
                        moire%sup_cell_info(i,1), moire%sup_cell_info(i,2), &
                        moire%sup_cell_info(i,3)
        CALL debug_output(0)
      END DO
#endif
    CASE ("molecular","Molecular","MOLECULAR")
      ALLOCATE(moire%lay_types(moire%natom))
      DO i=1,moire%natom
        READ(1,*) temp, moire%lay_types(i), moire%at_types_i(i), moire%real_pos(i,1), &
                        moire%real_pos(i,2), moire%real_pos(i,3), &
                        moire%sup_cell_info(i,1), moire%sup_cell_info(i,2), &
                        moire%sup_cell_info(i,3)
      END DO 
#ifdef __DEBUG
      WRITE(debug_str, '(A)') '\r\nLayer    Atom type         X          Y          Z'
      CALL debug_output(0)
      DO i = 1,moire%natom
        WRITE(debug_str, '(2I8, 3F12.6, 3I6)') moire%lay_types(i), moire%at_types_i(i), &
                              moire%real_pos(i,1), moire%real_pos(i,2), moire%real_pos(i,3), &
                        moire%sup_cell_info(i,1), moire%sup_cell_info(i,2), &
                        moire%sup_cell_info(i,3)
        CALL debug_output(0)
      END DO
#endif
    CASE ("full", "Full", "FULL")
      ALLOCATE(moire%lay_types(moire%natom))
      DO i=1,moire%natom
        READ(1,*) temp, moire%lay_types(i), moire%at_types_i(i), temp_d, &
                        moire%real_pos(i,1), moire%real_pos(i,2), moire%real_pos(i,3), &
                        moire%sup_cell_info(i,1), moire%sup_cell_info(i,2), &
                        moire%sup_cell_info(i,3)
      END DO
#ifdef __DEBUG
      WRITE(debug_str, '(A)') '\r\nLayer    Atom type         X          Y          Z'
      CALL debug_output(0)
      DO i = 1,moire%natom
        WRITE(debug_str, '(2I8, 3F12.6, 3I6)') moire%lay_types(i), moire%at_types_i(i), &
                              moire%real_pos(i,1), moire%real_pos(i,2), moire%real_pos(i,3), &
                        moire%sup_cell_info(i,1), moire%sup_cell_info(i,2), &
                        moire%sup_cell_info(i,3)
        CALL debug_output(0)
      END DO
#endif
    CASE DEFAULT
        write(err_msg,'(A)') "Atom style not supported. Supported styles are atomic and molecular"
        CALL error_message()
        write(err_msg,'(A)') "Cannot read LAMMPS file"
        CALL error_message()
        CALL EXIT
  END SELECT

  ! --------
  ! Get crys
  ! --------

  DO i=1,moire%natom
    CALL linsolve(moire%real_pos(i,:),moire%lat,moire%crys(i,:))
  END DO

  !DO i=1,moire%natom
  !  DO j=1,3
  !      moire%crys(i,j) = moire%crys(i,j) + moire%sup_cell_info(i,j)
  !      moire%real_pos(i,j) = moire%real_pos(i,j) + &
  !                            moire%sup_cell_info(i,1)*moire%lat(1,j) + &
  !                            moire%sup_cell_info(i,2)*moire%lat(2,j) + &
  !                            moire%sup_cell_info(i,3)*moire%lat(3,j)
  !  END DO
  !END DO

#ifdef __DEBUG
      WRITE(debug_str, '(A)') '\r\n Atom type         X          Y          Z'
      CALL debug_output(0)
      DO i = 1,moire%natom
        WRITE(debug_str, '(I8, 3F12.6)') moire%at_types_i(i), moire%crys(i,1), &
                                         moire%crys(i,2), moire%crys(i,3)
        CALL debug_output(0)
      END DO
#endif
  

  
  ! -----------------------
  ! Get crystal coordinates
  ! -----------------------

  CALL matinvtransp3(moire%lat, moire%rec_lat)

#ifdef __DEBUG
  WRITE (debug_str,'(A)') '\r\nReciprocal Lattice vectors: '
  CALL debug_output(0)
  DO i=1,3
    WRITE (debug_str, '(3F12.4)') moire%rec_lat(i,1), moire%rec_lat(i,2), moire%rec_lat(i,3)
    CALL debug_output(0)
  END DO
#endif

  CLOSE(unit=1)
   
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
