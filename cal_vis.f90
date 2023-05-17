subroutine cal_vis(ni, nj, visx, visy, u, v, dxinv, dyinv, muL, muG, rhoL, rhoG, rho)
    implicit none
    integer :: ni, nj
    double precision :: visx(-6:ni+7, -6:nj+7)
    double precision :: visy(-6:ni+7, -6:nj+7)
    double precision :: rho(-6:ni+7, -6:nj+7)
    double precision :: u(-6:ni+7, -6:nj+7)
    double precision :: v(-6:ni+7, -6:nj+7)
    double precision :: dxinv, dyinv
    double precision :: muL, muG, mudif
    double precision :: rhoL, rhoG, rhodifinv
    integer :: i, j
    double precision :: mu
    double precision :: rhomid
    double precision :: inv16


    mudif = muL - muG
    rhodifinv = 1d0/(rhoL - rhoG)
    inv16 = 1.0d0/16.0d0

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, mu, rhomid)

    do j = -5, nj+6
        do i = -5, ni+6

            rhomid = inv16* &
                    & (-rho(i-1, j)+9.0d0*rho(i, j)+9.0d0*rho(i+1, j)-rho(i+2, j))
            mu = muG + mudif*rhodifinv*(rho(i,j)-rhoG)
            visx(i, j) = + dxinv**2*(           &
                            + u(i-1, j )        &
                            - u(i  , j )*2.0d0  &
                            + u(i+1, j ))       &
                          + dyinv**2*(          &
                            + u(i  ,j-1)        &
                            - u(i  ,j  )*2.0d0  &
                            + u(i  ,j+1))
        end do
    end do

!$OMP  END PARALLEL DO




!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)

    do j = -5, nj+6
        do i = -5, ni+6
            visy(i, j) = + dxinv**2*(           &
                            + v(i-1,j  )        &
                            - v(i  ,j  )*2.0d0  &
                            + v(i+1,j  ))       &
                          + dyinv**2*(          &
                            + v(i  ,j-1)        &
                            - v(i  ,j  )*2.0d0  &
                            + v(i  ,j+1))
        end do
    end do

!$OMP  END PARALLEL DO



end subroutine cal_vis
