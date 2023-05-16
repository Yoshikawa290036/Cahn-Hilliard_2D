subroutine cal_adv(ni, nj, u, v, advx, advy, dxinv, dyinv)
    implicit none
    integer :: ni, nj
    double precision :: u(-6:ni+7, -6:nj+7)
    double precision :: v(-6:ni+7, -6:nj+7)
    double precision :: advx(-6:ni+7, -6:nj+7)
    double precision :: advy(-6:ni+7, -6:nj+7)

    double precision :: dxinv
    double precision :: dyinv
    
    integer :: i, j
    double precision :: inv12

    inv12 = 1.0d0/12.0d0


!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)

    do j = -4, nj+5
        do i = -4, ni+5
            advx(i,j) = + inv12*dxinv*( &
                        & + u(i-2,j)*u(i-2,j)         &
                        & - u(i-1,j)*u(i-1,j) * 8.0d0 &
                        & + u(i+1,j)*u(i+1,j) * 8.0d0 &
                        & - u(i+2,j)*u(i+2,j)         &
                        & ) &
                        + inv12*dyinv*( &
                        & + u(i,j-2)*v(i,j-2)         &
                        & - u(i,j-1)*v(i,j-1) * 8.0d0 &
                        & + u(i,j+1)*v(i,j+1) * 8.0d0 &
                        & - u(i,j+2)*v(i,j+2)         &
                        & )
        end do
    end do
!$OMP  END PARALLEL DO


!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)

    do j = -4, nj+5
        do i = -4, ni+5
            advy(i,j) = + inv12*dxinv*( &
                        & + u(i-2,j)*v(i-2,j)         &
                        & - u(i-1,j)*v(i-1,j) * 8.0d0 &
                        & + u(i+1,j)*v(i+1,j) * 8.0d0 &
                        & - u(i+2,j)*v(i+2,j)         &
                        & ) &
                        + inv12*dyinv*( &
                        & + v(i,j-2)*v(i,j-2)         &
                        & - v(i,j-1)*v(i,j-1) * 8.0d0 &
                        & + v(i,j+1)*v(i,j+1) * 8.0d0 &
                        & - v(i,j+2)*v(i,j+2)         &
                        & )
        end do
    end do

!$OMP  END PARALLEL DO

end subroutine cal_adv
