subroutine calphi(ni, nj, u, v, eta, dxinv, dyinv, phi, dt, a, b, temperature, kappa)
    implicit none
    integer :: ni, nj
    double precision :: dxinv, dyinv
    double precision ::   phi(-6:ni+7, -6:nj+7)
    double precision ::   eta(-6:ni+7, -6:nj+7)
    double precision ::     u(-6:ni+7, -6:nj+7)
    double precision ::     v(-6:ni+7, -6:nj+7)
    double precision :: phipx(-6:ni+7, -6:nj+7), dphix, dphipx
    double precision :: phipy(-6:ni+7, -6:nj+7), dphiy, dphipy

    double precision :: ddphix, ddphiy, ddphi(-6:ni+7, -6:nj+7)
    ! double precision :: nphi1(-6:ni+7, -6:nj+7), nphi2(-6:ni+7, -6:nj+7)
    double precision :: dt
    double precision :: a, b, temperature, kappa
    integer :: i, j
    double precision :: advx(-6:ni+7, -6:nj+7), chpx(-6:ni+7, -6:nj+7)
    double precision :: advy(-6:ni+7, -6:nj+7), chpy(-6:ni+7, -6:nj+7)
    double precision ::  Jpx(-6:ni+7, -6:nj+7), zetapx
    double precision ::  Jpy(-6:ni+7, -6:nj+7), zetapy
    double precision :: alpha
    double precision :: m1x, m2x, p0x, p1x, p2x
    double precision :: m1y, m2y, p0y, p1y, p2y
    double precision :: dddphipx, dchpx
    double precision :: dddphipy, dchpy
    double precision :: inv12, inv16, inv24
    inv12 = 1.0d0/12.0d0
    inv16 = 1.0d0/16.0d0
    inv24 = 1.0d0/24.0d0
    alpha = 12.0d0

    ! ddphi = nabla nabla phi
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, m1x, m2x, p0x, p1x, p2x, m1y, m2y, p0y, p1y, p2y, ddphix, ddphiy)
    do j = -4, nj+5
        do i = -4, ni+5
            m2x = -inv12*dxinv**2*phi(i-2, j)
            m1x =  inv12*dxinv**2*phi(i-1, j)*16.0d0
            p0x = -inv12*dxinv**2*phi(i  , j)*30.0d0
            p1x =  inv12*dxinv**2*phi(i+1, j)*16.0d0
            p2x = -inv12*dxinv**2*phi(i+2, j)
            ddphix = m2x+m1x+p0x+p1x+p2x

            m2y = -inv12*dxinv**2*phi(i, j-2)
            m1y =  inv12*dxinv**2*phi(i, j-1)*16.0d0
            p0y = -inv12*dxinv**2*phi(i, j  )*30.0d0
            p1y =  inv12*dxinv**2*phi(i, j+1)*16.0d0
            p2y = -inv12*dxinv**2*phi(i, j+2)
            ddphiy = m2y+m1y+p0y+p1y+p2y
            ddphi(i,j) = ddphix + ddphiy
        end do
    end do
!$OMP  END PARALLEL DO


! cal x direction of advection and chemical potential term
! ========================================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -6, nj+7
        do i = -5, ni+5
            phipx(i, j) = inv16* &
                    & (-phi(i-1, j)+9.0d0*phi(i, j)+9.0d0*phi(i+1, j)-phi(i+2, j))
            ! phipx(i, j) = 0.5d0*(phi(i, j)+phi(i+1, j))
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, zetapx, dphipx, dddphipx)
    do j = -4, nj+5
        do i = -2, ni+3
            zetapx = (temperature/((1.0d0-b*phipx(i, j))**2))-2.0d0*a*phipx(i, j)
            call nabla(dxinv, dphipx  , inv24, phi(i-1,j)  , phi(i,j)  , phi(i+1,j)  , phi(i+2,j)  )
            call nabla(dxinv, dddphipx, inv24, ddphi(i-1,j), ddphi(i,j), ddphi(i+1,j), ddphi(i+2,j))

            Jpx(i, j) = -zetapx*dphipx+kappa*phipx(i, j)*dddphipx
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, dphix, dchpx)
    do j = -1, nj+3
        do i = 0, ni+2
            call nabla(dxinv, advx(i,j), inv24, &
                      & 0.5d0*(u(i-2,j)+u(i-1,j))*phipx(i-2,j), &
                      & 0.5d0*(u(i-1,j)+u(i  ,j))*phipx(i-1,j), &
                      & 0.5d0*(u(i  ,j)+u(i+1,j))*phipx(i  ,j), &
                      & 0.5d0*(u(i+1,j)+u(i+2,j))*phipx(i+1,j))
            ! advx(i, j) = u*dphix

            call nabla(dxinv, dchpx, inv24, &
                      & phipx(i-2,j)*Jpx(i-2,j), &
                      & phipx(i-1,j)*Jpx(i-1,j), &
                      & phipx(i  ,j)*Jpx(i  ,j), &
                      & phipx(i+1,j)*Jpx(i+1,j))
            chpx(i,j) = alpha*dchpx

            ! nphi1(i,j) = phi(i,j)-dt*(advx(i,j)+chpx(i,j))
        end do
    end do
!$OMP  END PARALLEL DO

! cal y direction of advection and chemical potential term
! ========================================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)
    do j = -5, nj+5
        do i = -6, ni+7
            phipy(i, j) = inv16* &
                    & (-phi(i, j-1)+9.0d0*phi(i, j)+9.0d0*phi(i, j+1)-phi(i, j+2))
            ! phipy(i, j)= 0.5d0*(phi(i, j)+phi(i, j+1))
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, zetapy, dphipy, dddphipy)
    do j = -2, nj+3
        do i = -4, ni+5
            zetapy = (temperature/((1.0d0-b*phipy(i, j))**2))-2.0d0*a*phipy(i, j)
            call nabla(dyinv, dphipy  , inv24, phi(i,j-1)  , phi(i,j)  , phi(i,j+1)  , phi(i,j+2)  )
            call nabla(dyinv, dddphipy, inv24, ddphi(i,j-1), ddphi(i,j), ddphi(i,j+1), ddphi(i,j+2))

            Jpy(i, j) = -zetapy*dphipy+kappa*phipy(i, j)*dddphipy
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j, dphiy, dchpy)
    do j = 0, nj+2
        do i = -1, ni+3
            call nabla(dyinv, advy(i,j), inv24,  &
                      & 0.5d0*(v(i,j-2)+v(i,j-1))*phipy(i,j-2), &
                      & 0.5d0*(v(i,j-1)+v(i,j  ))*phipy(i,j-1), &
                      & 0.5d0*(v(i,j  )+v(i,j+1))*phipy(i,j  ), &
                      & 0.5d0*(v(i,j+1)+v(i,j+2))*phipy(i,j+1))
            ! advy(i, j) = v*dphiy

            call nabla( dyinv, dchpy, inv24, &
                      & phipy(i,j-2)*Jpy(i,j-2), &
                      & phipy(i,j-1)*Jpy(i,j-1), &
                      & phipy(i,j  )*Jpy(i,j  ), &
                      & phipy(i,j+1)*Jpy(i,j+1))
            chpy(i,j) = alpha*dchpy

            ! nphi1(i,j) = phi(i,j)-dt*(advy(i,j)+chpy(i,j))
        end do
    end do
!$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do j = 0, nj+1
        do i = 0, ni+1
            phi(i,j) = phi(i,j)-dt*(advx(i,j)+advy(i,j)+chpx(i,j)+chpy(i,j))
            ! phi(i,j) = nphi1(i,j)
        end do
    end do
!$OMP  END PARALLEL DO

end subroutine calphi
