module bz_integration
    
    ! ============================================================== 
    ! 
    ! Module for evaluating integrals of the form:
    !
    ! I(E) = $\int_{BZ} dk f(\epsilon_k) * \delta(E-\epsilon_k) 
    ! 
    ! ==============================================================
    
    contains
    subroutine linear_triangulation(nbands, loc_tri, len_E, Energy, freq, func, &
                                    rank_, dos_loc, ndos_loc)
        implicit none
        
        integer, intent(in) :: nbands, rank_, len_E
        integer, dimension(:,:), intent(in) :: loc_tri
        double precision, dimension(:,:), intent(in) :: freq, func  ! dim = (nkpt, nbands) 
        double precision, dimension(:), intent(in) :: Energy
        double precision, dimension(len_E), intent(out) :: dos_loc, ndos_loc
        
        integer :: len_tri
        double precision :: d, en, nd
        integer :: i,j,k
        double precision, dimension(4,2) :: temp_arr
    
        len_tri = size(loc_tri, DIM=1)  ! Number of triangles
        do i=1,len_E
            do j=1,len_tri
                do k=1,nbands
                    temp_arr(1,1) = Energy(i)
                    temp_arr(1,2) = 0.0
                    temp_arr(2,1) = freq(loc_tri(j,1), k)
                    temp_arr(2,2) = func(loc_tri(j,1), k)  
                    temp_arr(3,1) = freq(loc_tri(j,2), k)
                    temp_arr(3,2) = func(loc_tri(j,2), k)   
                    temp_arr(4,1) = freq(loc_tri(j,3), k)
                    temp_arr(4,2) = func(loc_tri(j,3), k)  
                    call sort_freq(temp_arr)
                    call integrate_triangle(temp_arr, Energy(i), d, nd)
                    dos_loc(i) = dos_loc(i) + d
                    ndos_loc(i) = ndos_loc(i) + nd 
                end do    
            end do
            if (rank_.eq.0) call progress(i,len_E) 
        end do
        return
    end subroutine
    
    subroutine sort_freq(a)
        implicit none
        double precision, intent(inout) :: a(4, 2)
        integer :: i, j
        double precision :: tmp(2)
        do i = 1,4
           do j=i+1,4
               if (a(i,1).gt.a(j,1)) then
                   tmp = a(i,:)
                   a(i,:) = a(j,:)
                   a(j,:) = tmp
               end if
           end do
        end do
        return
    end subroutine sort_freq
    
    
    subroutine integrate_triangle(temp_array, Ef, dos, nd)
    
        implicit none
        double precision, dimension(4,2), intent(in) :: temp_array
        double precision, intent(in) :: Ef
        double precision, intent(out) :: dos, nd
        double precision :: e0, e1, e2, f0, f1, f2
        double precision :: del10, del20, ddel10, ddel20
        double precision :: del02, del12, ddel02, ddel12
    
        if (temp_array(4,1) == Ef) then
            f0 = temp_array(1,2)
            f1 = temp_array(2,2)
            f2 = temp_array(3,2)
            dos = 0.0
            nd = f0+f1+f2
        end if 
        if (temp_array(3,1) == Ef) then
            
            e0 = temp_array(1,1)
            f0 = temp_array(1,2)
            e1 = temp_array(2,1)
            f1 = temp_array(2,2)
            e2 = temp_array(4,1)
            f2 = temp_array(4,2)
    
            del02 = (e2-Ef)/(e2-e0)
            del12 = (e2-Ef)/(e2-e1)
    
            ddel02 = 1/(e0-e2)
            ddel12 = 1/(e1-e2)
    
            nd = (1-del12*del02**2)*f0 + &
                 (1-del02*del12**2)*f1 + &
                 (1-3*del12*del02+del12*del02**2+del02*del12**2)*f2
            dos = (-2*del02*ddel02*del12 - ddel12*del02**2)*f0 + &
                  (-ddel02*del12**2 - 2*del02*del12*ddel12)*f1 + &
                  (-3*del12*ddel02 - 3*del02*ddel12 + 2*del02*del12*ddel02 + &
                    ddel12*del02**2 + 2*del02*del12*ddel12 + ddel02*del12**2)*f2
        end if
        if (temp_array(2,1) == Ef) then
            
            e0 = temp_array(1,1)
            f0 = temp_array(1,2)
            e1 = temp_array(3,1)
            f1 = temp_array(3,2)
            e2 = temp_array(4,1)
            f2 = temp_array(4,2)
    
            del10 = (Ef-e0)/(e1-e0)
            del20 = (Ef-e0)/(e2-e0)
            ddel10 = 1/(e1-e0)
            ddel20 = 1/(e2-e0)
            nd = (3-del10-del20)*f0 + del10*f1 + del20*f2
            dos = ddel10*del20*nd + del10*ddel20*nd + &
                  del10*del20*((-ddel10-ddel20)*f0+ddel10*f1+ddel20*f2)
            nd = nd*del10*del20
        end if
        if (temp_array(1,1)==Ef) then
            dos = 0.0
            nd = 0.0
        end if
    
        return
    
    end subroutine
    
    subroutine progress(i,n)
        implicit none
        integer(kind=4) :: i,n,k, frac,l
        character(len=32)::bar="???% |                         |"
        frac = int(i*100/n)
        write(unit=bar(1:3),fmt="(i3)")frac
        l = int(i*25/n)
        do k=1, l
            bar(6+k:6+k)="="
        enddo 
        write(unit=6,fmt="(a1,a40,a,$)") char(13), bar, 'of calculations done.'
        return
    end subroutine progress
end module              
        

            
