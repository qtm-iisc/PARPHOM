SUBROUTINE get_array_size(num_ele_row, num_ele_col, size_needed, lld, Locq)

    ! Get array size to be allocated in each group for a 
    ! (num_ele_row*num_ele_col) matrix distributed in block-cyclic fashion

    USE global_variables

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: num_ele_row, num_ele_col
    INTEGER, INTENT(OUT) :: size_needed, lld, Locq
    INTEGER :: rsrc, csrc, Locp
    EXTERNAL :: numroc
    INTEGER :: numroc

    rsrc = 0
    csrc = 0

    Locp = numroc(num_ele_row, pzheevx_vars%nb, blacs_grid%myprow, rsrc, blacs_grid%nprow)
    lld = max(1,Locp)
    Locq = numroc(num_ele_col, pzheevx_vars%mb, blacs_grid%mypcol, csrc, blacs_grid%npcol)
    Locq = max(1,Locq)

    size_needed = 1+lld*Locq

    if (size_needed.le.0) stop "Memory will overflow while reading force constant"

    return

END SUBROUTINE


subroutine first_array_index(descriptor, st_i, st_j)

    ! Computes the first array index on a local processor in a group

    use global_variables

    implicit none
    integer, dimension(DLEN_), intent(in) :: descriptor
    integer, intent(out) :: st_i, st_j

    if (blacs_grid%myprow.ge.descriptor(RSRC_)) THEN
        st_i = (blacs_grid%myprow - descriptor(RSRC_))*descriptor(MB_) + 1
    ELSE
        st_i = (blacs_grid%myprow + (blacs_grid%nprow-descriptor(RSRC_)))*descriptor(MB_) + 1
    ENDIF


    IF (blacs_grid%mypcol.ge.descriptor(CSRC_)) THEN
        st_j = (blacs_grid%mypcol - descriptor(CSRC_))*descriptor(NB_) + 1
    ELSE
        st_j = (blacs_grid%mypcol + (blacs_grid%npcol-descriptor(CSRC_)))*descriptor(NB_) + 1
    ENDIF

    return

end subroutine
