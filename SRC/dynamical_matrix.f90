subroutine create_dynamical_matrix()

    use global_variables

    implicit none
    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend, ia, ja
    integer :: lroffset, lcoffset, rsrc, csrc, lrindx, lcindx, ipos
    double precision, dimension(3) :: rc_i, rc_j, q_pt
    double precision, dimension(1) :: phi_ij
    double precision :: m_ij
    double complex :: d_ij

    q_pt = q_file%points(q_index,:)
    
    rsrc = 0
    csrc = 0

    call first_array_index(dynmat%desca, ia_first, ja_first)

    do jastart = ja_first, dynmat%desca(N_), dynmat%desca(NB_)*blacs_grid%npcol
        do iastart = ia_first, dynmat%desca(M_), dynmat%desca(MB_)*blacs_grid%nprow
            iaend = min(dynmat%desca(M_), iastart+dynmat%desca(MB_)-1)
            jaend = min(dynmat%desca(N_), iastart+dynmat%desca(NB_)-1)
            call infog2l(iastart, jastart, dynmat%desca,blacs_grid%nprow, &
                         blacs_grid%npcol, blacs_grid%myprow, blacs_grid%mypcol, &
                         lroffset, lcoffset, rsrc, csrc)
            do ja = jastart,jaend
                rc_j = moire%crys((ja-1)/3+1,:)
                do ia = iastart,iaend
                    rc_i = moire%crys((ia-1)/3+1,:)
                    lrindx = lroffset+(ia-iastart)
                    lcindx = lcoffset+(ja-jastart)
                    ipos = lrindx+(lcindx-1)*dynmat%desca(LLD_)
                    phi_ij(1) = fc%a(ipos)
                    if (ia.eq.ja) then
                        dynmat%a(ipos) = phi_ij(1)/moire%mass(moire%at_types_i((ia-1)/3+1))
                    else
                        call create_dynamical_ij(phi_ij, q_pt, rc_i, rc_j, moire%lat, d_ij) 
                        m_ij = sqrt(moire%mass(moire%at_types_i((ia-1)/3+1)) * &
                                    moire%mass(moire%at_types_i((ja-1)/3+1)))
                        dynmat%a(ipos) = d_ij/m_ij
                    end if
                end do
            end do
        end do
    end do

end subroutine



SUBROUTINE create_dynamical_ij(f_ij, q_pt, rci, rcj, lat, d_ij)
  
  ! Subroutine that creates the dynamical matrix element. 
  ! 1/sqrt(m1*m2) needs to be multiplied with the output
  
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(1), INTENT(IN) :: f_ij
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: rci, rcj, q_pt
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: lat
  DOUBLE COMPLEX, INTENT(OUT) :: d_ij
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793
  DOUBLE PRECISION, PARAMETER :: tol = 0.0000001
  DOUBLE COMPLEX, PARAMETER :: IM = (0,1)
  INTEGER :: s,t,l,i,j,ix, multiplicity
  DOUBLE PRECISION :: dist, minimum
  DOUBLE PRECISION, DIMENSION(3) :: rcj_sup, rri, rrcj_sup
  DOUBLE PRECISION, DIMENSION(9) :: neigh
  INTEGER, ALLOCATABLE, DIMENSION(:) :: b
  DOUBLE COMPLEX :: phase
  DOUBLE PRECISION :: prephase

  rri = MATMUL(TRANSPOSE(lat),rci)

  DO l=1,9
    i = (l-1)/3 - 1
    j = l - 2 - 3*(i+1)
    rcj_sup(1) = rcj(1) + i
    rcj_sup(2) = rcj(2) + j
    rcj_sup(3) = rcj(3)
    rcj_sup = MATMUL(TRANSPOSE(lat),rcj_sup)
    dist = sqrt((rri(1)-rcj_sup(1))**2 + &
                (rri(2)-rcj_sup(2))**2 + &
                (rri(3)-rcj_sup(3))**2)
    neigh(l) = dist
  END DO
  
  minimum = MINVAL(neigh)
  b = PACK([(ix, ix=1,SIZE(neigh))], neigh<minimum+tol)

  phase = (0,0)
  multiplicity = 0
  DO l=1,SIZE(b)
    i = (b(l)-1)/3 - 1
    j = b(l) - 2 - 3*(i+1)
    rcj_sup(1) = rcj(1) + i - rci(1)
    rcj_sup(2) = rcj(2) + j - rci(2)
    rcj_sup(3) = rcj(3) - rci(3)
    prephase = DOT_PRODUCT(q_pt, rcj_sup)
    phase  = phase + EXP(-2*PI*IM*prephase)
    multiplicity = multiplicity + 1
  END DO

  d_ij = f_ij(1)*phase / multiplicity

  RETURN  
END SUBROUTINE

