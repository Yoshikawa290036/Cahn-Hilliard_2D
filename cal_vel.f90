subroutine cal_vel(ni, nj, u, v, xl, yl, dx, dy)
    implicit none
    integer :: ni, nj
    double precision :: u(-6:ni+7, -6:nj+7)
    double precision :: v(-6:ni+7, -6:nj+7)
    double precision :: xl, yl, dx, dy
    integer :: i, j
    double precision :: x, y

    ! u = -y
    do j = -6, nj+7
        do i = -6, ni+7
            y = (dble(j)-0.5d0)*dy-(yl*0.5d0)
            u(i, j) = -y/(yl*0.5d0)
            ! u(i, j) = 0.0
        end do
    end do

    ! v = x
    do j = -6, nj+7
        do i = -6, ni+7
            x = (dble(i)-0.5d0)*dx-(xl*0.5d0)
            v(i, j) = x/(xl*0.5d0)
            ! v(i, j) = 0.35
        end do
    end do
end subroutine cal_vel
