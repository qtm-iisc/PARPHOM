subroutine create_dynamical_matrix(q_indx)

    use mpi
    use global_variables

    implicit none
    integer, intent(in) :: q_indx
    integer :: ia_first, ja_first, iastart, jastart, iaend, jaend , rsrc, csrc
    integer :: lroffset, lcoffset, ia, ja, lrindx, lcindx, ipos
#ifdef __DEBUG
    integer :: i,j
    character(len=50) :: format_
#endif


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
                    call create_dynamical_ij(ipos, q_indx, ia, ja)
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



#ifdef __KPOOL    
    call mpi_barrier(mpi_local%comm, mpierr)
    write(debug_str,'(A,I0,A)') "Dynamical Matrix created for q-pt no.", q_indx ," on "

    call date_time_message_local(trim(debug_str))
#else
    call mpi_barrier(mpi_global%comm, mpierr)

    write(debug_str,'(A,I0,A)') "Dynamical Matrix created for q-pt no. ", q_indx, " on "
    call date_time_message(trim(debug_str))
#endif

end subroutine


subroutine create_dynamical_ij(ipos, q_indx, inp_ia, inp_ja)
  use global_variables
  implicit none
  integer, intent(in) :: ipos, q_indx, inp_ia, inp_ja
  double precision, parameter :: PI  = 3.141592653589793
  double precision, parameter :: tol = 0.0000001
  double complex, parameter   :: IM  = (0,1)
  double precision, dimension(3) :: rj_sup
  integer :: l,i,j, ix, multiplicity, ia, ja
  double precision :: dist, minimum, prephase
  double precision, dimension(9) :: neigh
  integer, allocatable, dimension(:) :: b
  double complex :: phase

  ia = int((inp_ia+2)/3)
  ja = int((inp_ja+2)/3)
  do l= -1,9
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
    prephase = dot_product(q_file%points(q_indx,:),rj_sup)
    phase = phase + exp(-2*PI*IM*prephase)
    multiplicity = multiplicity + 1
  end do

  dyn_mat%mat(ipos) = force_const%mat(ipos)*phase/multiplicity
  dyn_mat%mat(ipos) = dyn_mat%mat(ipos)/sqrt(moire%mass(moire%at_types_i(ia))*&
                                             moire%mass(moire%at_types_i(ja)))

  return

end subroutine
