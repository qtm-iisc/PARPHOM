subroutine create_dynamical_matrix(q_indx, derivative)

    use mpi
    use global_variables

    implicit none
    integer, intent(in) :: derivative, q_indx
    double precision, dimension(3) :: delta, delta_crys, q_shifted
    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend , rsrc, csrc
    integer :: lroffset, lcoffset, ia, ja, lrindx, lcindx, ipos
    double complex :: dij, dijplus, dijminus
    double precision, parameter :: DEL_INCREMENT = 1e-6

#ifdef __DEBUG
    integer :: i,j
    character(len=50) :: format_
#endif


    ! -------------------------------------------------------------------------------
    ! Subroutine that creates the dynamical matrix in a block cyclically distributed 
    ! manner across the processors. The variable derivative is passed to compute the
    ! derivative of the dynamical matrix.
    !
    !           derivative = 0 :: no derivative is computed
    !           derivative = 1 :: dD/dkx is computed
    !           derivative = 2 :: dD/dky is computed
    !           derivative = 3 :: dD/dkz is computed
    !   
    !     We provide two options to compute the derivatives of the dynamical matrix.
    !
    !     The first option is to use the central finite difference is used for 
    !     computing the derivatives.
    !     The delta is taken as (0.00001 1/Ang) in the reciprocal lattice,
    !     in each direction.
    !
    !     The second option is to use the analytical derivative of the dynamical matrix
    !     
    !           (dD/dq)_(i,j) = 1/sqrt(m_i*m_j)*\phi_{\alpha,\beta}^{i,j}*
    !                           exp(1j*dot(q,(r_j - r_i)) * (i*(r_j - r_i))
    !
    !       
    ! -------------------------------------------------------------------------------


    rsrc = 0
    csrc = 0
    
    if (grid%myprow.ge.dyn_mat%desca(RSRC_)) then
        ia_first = (grid%myprow - dyn_mat%desca(RSRC_))*dyn_mat%desca(MB_)+1
    else
        ia_first = (grid%myprow + (grid%nprow - dyn_mat%desca(RSRC_)))* &
                    dyn_mat%desca(MB_) + 1
    endif

    if (grid%mypcol.ge.dyn_mat%desca(CSRC_)) then
        ja_first = (grid%mypcol - dyn_mat%desca(CSRC_))*dyn_mat%desca(NB_)+1
    else
        ja_first = (grid%mypcol + (grid%npcol - dyn_mat%desca(CSRC_)))* &
                    dyn_mat%desca(NB_) + 1
    endif

    do jastart = ja_first, dyn_mat%desca(N_), grid%npcol*dyn_mat%desca(NB_)
        do iastart = ia_first, dyn_mat%desca(M_), grid%nprow*dyn_mat%desca(MB_)
            iaend = min(dyn_mat%desca(M_), iastart+dyn_mat%desca(MB_)-1)
            jaend = min(dyn_mat%desca(N_), jastart+dyn_mat%desca(NB_)-1)

            ia = iastart
            ja = jastart

            call infog2l (ia, ja, dyn_mat%desca, grid%nprow, grid%npcol, &
                          grid%myprow, grid%mypcol, lroffset, lcoffset, rsrc, csrc)

            do ja=jastart,jaend
                do ia=iastart,iaend
                    lrindx = lroffset + (ia-iastart)
                    lcindx = lcoffset + (ja-jastart)
                    ipos = lrindx + (lcindx-1)*dyn_mat%desca(LLD_)
                    if (derivative.ne.0) then
                        if (vel_method=='A') then
                            !call analytical_derivative_dynamical_ij(ipos, q_indx, &
                            !                                        ia, ja, derivative)
                            call create_dynamical_ij(ipos,q_file%points(q_indx,:),ia,ja,1, &
                                                 derivative, dij)
                        else 
                            delta=0.0
                            delta(derivative) = max(DEL_INCREMENT* &
                                                sqrt(abs(moire%rec_lat(1,1)* &
                                                moire%rec_lat(2,2) - &
                                                moire%rec_lat(1,2)*  &
                                                moire%rec_lat(2,1))), &
                                                DEL_INCREMENT)
                            call linsolve(delta, moire%rec_lat, delta_crys)
                            q_shifted = q_file%points(q_indx,:) + delta_crys/2
                            call create_dynamical_ij(ipos, q_shifted, ia, ja, 0, &
                                                     derivative, dijplus)
                            q_shifted = q_file%points(q_indx,:) - delta_crys/2
                            call create_dynamical_ij(ipos, q_shifted, ia, ja, 0, &
                                                     derivative, dijminus)
                            dij = (dijplus-dijminus)/(norm2(delta))    
                        end if
                    else
                        call create_dynamical_ij(ipos,q_file%points(q_indx,:),ia,ja,0, &
                                                 derivative, dij)
                    end if
                    dyn_mat%mat(ipos) = dij
                end do
            end do
        end do
    end do

#ifdef __DEBUG
    write(debug_str,'(A)') "The local dynamical matrix elements in each rank are: "
    call debug_output(0)
    do i=1,mpi_global%size_
        write(format_, '(A,I0,A)') '(A,I0,A,',dyn_mat%size_,'(F12.6,F12.6,6X))'
        if (mpi_global%rank==i-1) then
            write(*,format_) "Rank :", mpi_global%rank, &
                             "\r\nDynamical Matrix : \r\n", (real(dyn_mat%mat(j)), &
                             aimag(dyn_mat%mat(j)), j=1,dyn_mat%size_)
        end if
        call mpi_barrier(mpi_global%comm, mpierr)
    end do
    call mpi_barrier(mpi_global%comm, mpierr)
#endif



#ifdef __QPOOL    
    call mpi_barrier(mpi_local%comm, mpierr)
    write(debug_str,'(A,I0,A)') "Dynamical Matrix created for q-pt no.", q_indx ," on "

    call date_time_message_local(trim(debug_str))
#else
    call mpi_barrier(mpi_global%comm, mpierr)
    if (derivative==0) then
        write(debug_str,'(A,I0,A)') "Dynamical Matrix created for q-pt no. ", q_indx, " on "
        call date_time_message(trim(debug_str))
    end if
#endif

end subroutine


subroutine create_dynamical_ij(ipos, q_point, inp_ia, inp_ja, deriv_order, derivative, dij)
  
  ! -------------------------------------------------------------------
  ! Computes the dynamical matrix elements d_ij
  ! -------------------------------------------------------------------  
    
  use global_variables
  implicit none
  integer, intent(in) :: ipos, inp_ia, inp_ja, deriv_order, derivative
  double precision, dimension(3) :: q_point
  double complex, intent(out) :: dij
  double precision, parameter :: PI  = 3.141592653589793
  double precision, parameter :: tol = 1e-8
  double complex, parameter   :: IM  = (0.0000000000000000,1.0000000000000000)
  double precision, dimension(3) :: rj_sup
  integer :: l,i,j, ix, multiplicity, ia, ja
  double precision :: dist, minimum, prephase
  double precision, dimension(9) :: neigh
  integer, allocatable, dimension(:) :: b
  double complex :: phase

  ia = int((inp_ia+2)/3)
  ja = int((inp_ja+2)/3)
  do l= 1,9
    i = (l-1)/3 - 1
    j = l - 2 - 3*(i+1)
    rj_sup(1) = moire%crys(ja,1) + i
    rj_sup(2) = moire%crys(ja,2) + j
    rj_sup(3) = moire%crys(ja,3)
    rj_sup = matmul(transpose(moire%lat),rj_sup)
    dist = sqrt((moire%real_pos(ia,1)-rj_sup(1))**2 + &
                (moire%real_pos(ia,2)-rj_sup(2))**2 + &
                (moire%real_pos(ia,3)-rj_sup(3))**2)
    neigh(l) = dist
  end do

  minimum = minval(neigh)
  b = pack([(ix, ix=1,size(neigh))], neigh<minimum+tol)

  phase = cmplx(0,0)
  multiplicity = 0
  do l=1,size(b)
    i = (b(l)-1)/3 - 1
    j = b(l)-2-3*(i+1)
    rj_sup(1) = moire%crys(ja,1) + i - moire%crys(ia,1)
    rj_sup(2) = moire%crys(ja,2) + j - moire%crys(ia,2)
    rj_sup(3) = moire%crys(ja,3) - moire%crys(ia,3)
    prephase = dot_product(q_point,rj_sup)
    if (deriv_order==0) then
        phase = phase + exp(-2*PI*IM*prephase)
    else if (deriv_order==1) then
        rj_sup = matmul(transpose(moire%lat),rj_sup)
        phase = phase - exp(-2*PI*IM*prephase)*IM*rj_sup(derivative)
    end if
    multiplicity = multiplicity + 1
  end do

  dij = force_const%mat(ipos)*phase/multiplicity
  dij = dij/sqrt(moire%mass(moire%at_types_i(ia))*moire%mass(moire%at_types_i(ja)))

  return

end subroutine


!subroutine analytical_derivative_dynamical_ij(ipos, q_indx, inp_ia, inp_ja, derivative)
!  
!  ! -----------------------------------------------------------------
!  ! Computes the analytical derivative of the dynamical matrix
!  ! -----------------------------------------------------------------
!    
!  use global_variables
!  implicit none
!  integer, intent(in) :: ipos, inp_ia, inp_ja, q_indx, derivative
!  double precision, parameter :: PI  = 3.141592653589793
!  double precision, parameter :: tol = 0.0000001
!  double complex, parameter   :: IM  = (0,1)
!  double precision, dimension(3) :: rj_sup
!  integer :: l,i,j, ix, multiplicity, ia, ja
!  double precision :: dist, minimum, prephase
!  double precision, dimension(9) :: neigh
!  integer, allocatable, dimension(:) :: b
!  double complex :: phase
!
!  ia = int((inp_ia+2)/3)
!  ja = int((inp_ja+2)/3)
!  do l= -1,9
!    i = (l-1)/3 - 1
!    j = l - 2 - 3*(i+1)
!    rj_sup(1) = moire%crys(ja,1) + i
!    rj_sup(2) = moire%crys(ja,2) + j
!    rj_sup(3) = moire%crys(ja,3)
!    rj_sup = matmul(transpose(moire%lat),rj_sup)
!    dist = sqrt((moire%real_pos(ia,1)-rj_sup(1))**2 + &
!                (moire%real_pos(ia,2)-rj_sup(2))**2 + &
!                (moire%real_pos(ia,3)-rj_sup(3))**2)
!    neigh(l) = dist
!  end do
!
!  minimum = minval(neigh)
!  b = pack([(ix, ix=1,size(neigh))], neigh<minimum+tol)
!
!  phase = cmplx(0,0)
!  multiplicity = 0
!  do l=1,size(b)
!    i = (b(l)-1)/3 - 1
!    j = b(l)-2-3*(i+1)
!    rj_sup(1) = moire%crys(ja,1) + i - moire%crys(ia,1)
!    rj_sup(2) = moire%crys(ja,2) + j - moire%crys(ia,2)
!    rj_sup(3) = moire%crys(ja,3) - moire%crys(ia,3)
!    prephase = dot_product(q_file%points(q_indx,:),rj_sup)
!    rj_sup = matmul(transpose(moire%lat),rj_sup)
!    phase = phase - exp(-2*PI*IM*prephase)*IM*rj_sup(derivative)
!    multiplicity = multiplicity + 1
!  end do
!
!  dyn_mat%mat(ipos) = force_const%mat(ipos)*phase/multiplicity
!  dyn_mat%mat(ipos) = dyn_mat%mat(ipos)/sqrt(moire%mass(moire%at_types_i(ia))*&
!                                             moire%mass(moire%at_types_i(ja)))
!  return
!
!end subroutine








