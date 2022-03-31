module neighbor_list

    contains
    subroutine create_neigh_list(natom, real_pos, lat, neigh_list)
        implicit none
        integer, intent(in) :: natom
        double precision, dimension(natom,3), intent(in) :: real_pos
        double precision, dimension(3,3) :: lat
        double precision, dimension(natom,3,6), intent(out) :: neigh_list

        integer :: i, j, p1, p2
        double precision :: X,Y,Z,R, temp

        temp = (lat(1,2) - lat(1,1))*(lat(2,2)-lat(2,1))
        neigh_list = abs(temp)
        do i=1,natom
            do j=1,natom
                do p1 =-1,1
                    do p2 =-1,1
                        X = real_pos(i,1) - (real_pos(j,1)+p1*lat(1,1)+p2*lat(1,2))
                        Y = real_pos(i,2) - (real_pos(j,2)+p1*lat(2,1)+p2*lat(2,2))
                        Z = real_pos(i,3) - real_pos(j,3)
                        R = sqrt(X**2+Y**2+Z**2)
                        if (R.gt.0.01) then
                            if (R.le.neigh_list(i,1,1)) then
                                neigh_list(i,3,:) = neigh_list(i,2,:)
                                neigh_list(i,2,:) = neigh_list(i,1,:)
                                neigh_list(i,1,:) = (/ R,              &
                                                      dble(j),         &
                                                      dble(p1),        &
                                                      dble(p2),        &
                                                      real_pos(i,1)-X, &
                                                      real_pos(i,2)-Y  /)
                            else if (R.le.neigh_list(i,2,1)) then
                                neigh_list(i,3,:) = neigh_list(i,2,:)
                                neigh_list(i,2,:) = (/ R,              &
                                                      dble(j),         &
                                                      dble(p1),        &
                                                      dble(p2),        &
                                                      real_pos(i,1)-X, &
                                                      real_pos(i,2)-Y  /)
                            else if (R.le.neigh_list(i,3,1)) then
                                neigh_list(i,3,:) = (/ R,              &
                                                      dble(j),         &
                                                      dble(p1),        &
                                                      dble(p2),        &
                                                      real_pos(i,1)-X, &
                                                      real_pos(i,2)-Y  /) 
                            end if
                        end if
                    end do
                end do
            end do
        end do

        return
    end subroutine create_neigh_list

    subroutine calculate_strain(natom, neigh_list, equib_len, strain)
        implicit none
        integer, intent(in) :: natom
        double precision, dimension(natom, 3, 6), intent(in) :: neigh_list
        double precision, intent(in) :: equib_len
        double precision, dimension(natom, 3), intent(out) :: strain
        integer :: i


        do i=1,natom
            strain(i,1) = 100*(neigh_list(i,1,1)-equib_len)/equib_len 
            strain(i,2) = 100*(neigh_list(i,2,1)-equib_len)/equib_len 
            strain(i,3) = 100*(neigh_list(i,3,1)-equib_len)/equib_len 
        end do

        return
    end subroutine
end module
