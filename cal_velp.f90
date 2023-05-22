subroutine cal_velp(ni, nj, dt, rho, up, vp, &
                    advx, advy, Fsx, Fsy, visx, visy)
    implicit none

    integer :: ni, nj
    double precision :: dt
    double precision ::       u(-6:ni+7, -6:nj+7) 
    double precision ::       v(-6:ni+7, -6:nj+7) 
    double precision ::      up(-6:ni+7, -6:nj+7) 
    double precision ::      vp(-6:ni+7, -6:nj+7) 
    double precision ::    advx(-6:ni+7, -6:nj+7) 
    double precision ::    advy(-6:ni+7, -6:nj+7) 
    double precision ::    visx(-6:ni+7, -6:nj+7) 
    double precision ::    visy(-6:ni+7, -6:nj+7) 
    double precision ::     rho(-6:ni+7, -6:nj+7) 
    double precision ::     Fsx(-6:ni+7, -6:nj+7) 
    double precision ::     Fsy(-6:ni+7, -6:nj+7) 
    double precision :: gx, gy
    integer :: i, j

    gx = 0.0d0
    gy = -2.0e-3

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)

    do j = -6, nj+7
        do i = -6, ni+7
            up(i,j) = u(i,j) + dt/rho(i,j)*( &
                    - advx(i,j) + Fsx(i,j) + visx(i,j) + gx)
            
            vp(i,j) = v(i,j) + dt/rho(i,j)*( &
                    - advy(i,j) + Fsy(i,j) + visy(i,j) + gy)
        end do
    end do

!$OMP  END PARALLEL DO


end subroutine cal_velp
