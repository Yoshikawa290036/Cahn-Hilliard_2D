subroutine cal_vis(ni, nj, visx, visy, u, v, ua, va, &
                   dxinv, dyinv, muL, muG, rhoL, rhoG, rho)
    implicit none
    integer :: ni, nj
    double precision :: visx(-6:ni+7, -6:nj+7)
    double precision :: visy(-6:ni+7, -6:nj+7)
    double precision ::  rho(-6:ni+7, -6:nj+7)
    double precision ::    u(-6:ni+7, -6:nj+7)
    double precision ::    v(-6:ni+7, -6:nj+7)
    double precision ::   ua(-6:ni+7, -6:nj+7) 
    double precision ::   va(-6:ni+7, -6:nj+7) 
    double precision :: dxinv, dyinv
    double precision :: muL, muG, mudif
    double precision :: rhoL, rhoG, rhodifinv
    integer :: i, j
    double precision :: mu
    double precision :: rhomid
    double precision :: inv16
    double precision :: ddudxx, ddudyy
    double precision :: ddvdxx, ddvdyy
    double precision :: x, y


    mudif = muL - muG
    rhodifinv = 1d0/(rhoL - rhoG)
    inv16 = 1.0d0/16.0d0

    do j = -6, nj+7
        do i = -6, ni+7
            ! cal u ================
            x = (dble(i)-0.5d0)/dxinv
            y = (dble(j)      )/dyinv
            u(i,j) = sin(x*0.1d0)

            ! cal v ================
            x = (dble(i)      )/dxinv
            y = (dble(j)-0.5d0)/dyinv
            v(i,j) = cos(y*0.1d0)
        end do
    end do


!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, mu, ddudxx, ddudyy, ddvdxx, ddvdyy)

    do j = -4, nj+6
        do i = -4, ni+6

            mu = muG + mudif*rhodifinv*(rho(i,j)-rhoG)

            ddudxx = 0.5d0*dxinv**2*(  &
                     + u(i-2, j)       &
                     - u(i-1, j)       &
                     - u(i  , j)       &
                     + u(i+1, j))

            ddudyy = 0.5d0*dyinv**2*(  &
                     + ua(i-2, j)      &
                     - ua(i-1, j)      &
                     - ua(i  , j)      &
                     + ua(i+1, j))

            ddvdxx = 0.5d0*dxinv**2*(   &
                     + va(i-2, j)       &
                     - va(i-1, j)       &
                     - va(i  , j)       &
                     + va(i+1, j))

            ddvdyy = 0.5d0*dyinv**2*(  &
                     + v(i-2, j)       &
                     - v(i-1, j)       &
                     - v(i  , j)       &
                     + v(i+1, j))

            ! visx(i, j) = + dxinv**2*(           &
            !                 + u(i-1, j )        &
            !                 - u(i  , j )*2.0d0  &
            !                 + u(i+1, j ))       &
            !               + dyinv**2*(          &
            !                 + u(i  ,j-1)        &
            !                 - u(i  ,j  )*2.0d0  &
            !                 + u(i  ,j+1))
            visx(i, j) = mu*(ddudxx + ddudyy)
            visy(i, j) = mu*(ddvdxx + ddvdyy)
        end do
    end do

!$OMP  END PARALLEL DO




! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i, j, mu, ddvdxx, ddvdyy)

!     do j = -4, nj+6
!         do i = -4, ni+6

!             mu = muG + mudif*rhodifinv*(rho(i,j)-rhoG)

!             ddvdxx = 0.5d0*dxinv**2*(   &
!                      + va(i-2, j)       &
!                      - va(i-1, j)       &
!                      - va(i  , j)       &
!                      + va(i+1, j))

!             ddvdyy = 0.5d0*dyinv**2*(  &
!                      + v(i-2, j)       &
!                      - v(i-1, j)       &
!                      - v(i  , j)       &
!                      + v(i+1, j))

!             ! visy(i, j) = + dxinv**2*(           &
!             !                 + v(i-1,j  )        &
!             !                 - v(i  ,j  )*2.0d0  &
!             !                 + v(i+1,j  ))       &
!             !               + dyinv**2*(          &
!             !                 + v(i  ,j-1)        &
!             !                 - v(i  ,j  )*2.0d0  &
!             !                 + v(i  ,j+1))
!             visy(i, j) = mu*(ddvdxx + ddvdyy)
!         end do
!     end do

! !$OMP  END PARALLEL DO


    do j = 1, nj
        do i = 1, ni
            x = (dble(i)-0.5d0)/dxinv
            y = (dble(j)-0.5d0)/dyinv
            mu = muG + mudif*rhodifinv*(rho(i,j)-rhoG)
            write (33, '(20e20.10)') x, y, visx(i,j)/mu, visy(i,j)/mu
        end do
        write (33, *)
    end do


end subroutine cal_vis
