subroutine cross_product(a,b,cross)

    ! --------------    
    ! cross = a X b    
    ! --------------

    implicit none
    double precision, intent(out), dimension(3) :: cross
    double precision, intent(in), dimension(3) :: a,b

    if ((a(1)*b(2) - a(2)*b(1))>=0) then
        cross(1) = a(2)*b(3) - a(3)*b(2)
        cross(2) = a(3)*b(1) - a(1)*b(3)
        cross(3) = a(1)*b(2) - a(2)*b(1)
    else
        cross(1) = a(3)*b(2) - a(2)*b(3)
        cross(2) = a(1)*b(3) - a(3)*b(1)
        cross(3) = a(2)*b(1) - a(1)*b(2)
    end if
end subroutine
