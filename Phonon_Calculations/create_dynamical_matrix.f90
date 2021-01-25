! Subroutine that creates the dynamical matrix element. 
! 1/sqrt(m1*m2) needs to be multiplied with the output


SUBROUTINE create_dynamical_ij(f_ij, q_pt, rci, rcj, lat, d_ij)
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

  DO l=-1,9
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

