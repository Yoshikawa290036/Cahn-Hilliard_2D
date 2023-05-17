subroutine cal_Fs(ni, nj, dxinv, dyinv, rho, Fsx, Fsy, kappa_s)
    implicit none
    integer :: ni, nj
    double precision :: kappa_s, dxinv, dyinv
    double precision ::  rho(-6:ni+7, -6:nj+7)
    double precision ::  Fsx(-6:ni+7, -6:nj+7)
    double precision ::  Fsy(-6:ni+7, -6:nj+7)
    double precision ::    drhodx(-6:ni+7, -6:nj+7)
    double precision ::    drhody(-6:ni+7, -6:nj+7)
    double precision ::  ddrhodxx(-6:ni+7, -6:nj+7)
    double precision ::  ddrhodyy(-6:ni+7, -6:nj+7)
    double precision :: ddrhodxdy(-6:ni+7, -6:nj+7)
    double precision :: rhopx(-6:ni+7, -6:nj+7)
    double precision :: rhopy(-6:ni+7, -6:nj+7)
    integer :: i, j
    double precision :: inv24, inv16, inv12
    double precision :: m2, m1, p0, p1, p2
    double precision :: x, y

    inv12 = 1.0d0/12.0d0
    inv24 = 1.0d0/24.0d0
    inv16 = 1.0d0/16.0d0


    ! do j = -6, nj+7
    !     do i = -6, ni+7
    !         x = (dble(i)-0.5d0)/dxinv
    !         y = (dble(j)-0.5d0)/dyinv
    !         rho(i,j) = sin(x/10.0d0)*cos(y/10.0d0)
    !         ! rho(i,j) = sin(y/10.0d0)
    !     end do
    ! end do


! rho interpolation for x direction

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -6, nj+7
        do i = -5, ni+5
            rhopx(i, j) = inv16* &
                    & (-rho(i-1, j)+9.0d0*rho(i, j)+9.0d0*rho(i+1, j)-rho(i+2, j))
        end do
    end do
!$OMP  END PARALLEL DO

! rho interpolation for y direction

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -5, nj+5
        do i = -6, ni+7
            rhopy(i, j) = inv16* &
                    & (-rho(i, j-1)+9.0d0*rho(i, j)+9.0d0*rho(i, j+1)-rho(i, j+2))
        end do
    end do
!$OMP  END PARALLEL DO


! cal drho/dx & drho/dy
! ====================================================================================

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -6, nj+7
        do i = -3, ni+4
            call nabla(dxinv, drhodx(i, j), inv24, &
                      & rhopx(i-2, j), &
                      & rhopx(i-1, j), &
                      & rhopx(i  , j), &
                      & rhopx(i+1, j))
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -3, nj+4
        do i = -6, ni+7
            call nabla(dyinv, drhody(i, j), inv24, &
                      & rhopy(i, j-2), &
                      & rhopy(i, j-1), &
                      & rhopy(i, j  ), &
                      & rhopy(i, j+1))
        end do
    end do
!$OMP  END PARALLEL DO


! cal ddrho/dxx & ddrho/dyy
! ====================================================================================

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, m2, m1, p0, p1, p2)
    do j = -6, nj+7
        do i = -4, ni+5
            m2 = -inv12*dxinv**2*rho(i-2, j)
            m1 =  inv12*dxinv**2*rho(i-1, j)*16.0d0
            p0 = -inv12*dxinv**2*rho(i  , j)*30.0d0
            p1 =  inv12*dxinv**2*rho(i+1, j)*16.0d0
            p2 = -inv12*dxinv**2*rho(i+2, j)
            ddrhodxx(i, j) = m2+m1+p0+p1+p2
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, m2, m1, p0, p1, p2)
    do j = -4, nj+5
        do i = -6, ni+7
            m2 = -inv12*dyinv**2*rho(i, j-2)
            m1 =  inv12*dyinv**2*rho(i, j-1)*16.0d0
            p0 = -inv12*dyinv**2*rho(i, j  )*30.0d0
            p1 =  inv12*dyinv**2*rho(i, j+1)*16.0d0
            p2 = -inv12*dyinv**2*rho(i, j+2)
            ddrhodyy(i, j) = m2+m1+p0+p1+p2
        end do
    end do
!$OMP  END PARALLEL DO


! cal ddrho/dxdy
! ====================================================================================

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, m2, m1, p0, p1, p2)
    do j = -4, nj+5
        do i = -4, ni+5
            m2 =  inv12*dyinv*drhodx(i, j-2)
            m1 = -inv12*dyinv*drhodx(i, j-1)*8.0d0
            p0 =  0.0d0
            p1 =  inv12*dyinv*drhodx(i, j+1)*8.0d0
            p2 = -inv12*dyinv*drhodx(i, j+2)
            ddrhodxdy(i, j)= m2+m1+p0+p1+p2
        end do
    end do
!$OMP  END PARALLEL DO


! cal Fsx, Fsy
! ====================================================================================

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -4, nj+5
        do i = -4, ni+5
            Fsx(i, j) = kappa_s*( &
                      &  drhody(i,j)*ddrhodxdy(i,j) &
                      & -drhodx(i,j)*ddrhodyy(i,j))
            
            Fsy(i, j) = kappa_s*( &
                      &  drhodx(i,j)*ddrhodxdy(i,j) &
                      & -drhody(i,j)*ddrhodxx(i,j)) 
        end do
    end do
!$OMP  END PARALLEL DO

    ! do j = 1, nj
    !     do i = 1, ni
    !         x = (dble(i)-0.5d0)/dxinv
    !         y = (dble(j)-0.5d0)/dyinv
    !         write (22, '(20e20.10)') x, y, Fsx(i,j), Fsy(i,j), rho(i, j), drhodx(i,j), drhody(i,j), cos(x)*cos(y)
    !         ! write (22, '(20e20.10)') x, dble(j), rhopy(i,j)
    !     end do
    !     write (22,*) 
    ! end do

end subroutine cal_Fs
