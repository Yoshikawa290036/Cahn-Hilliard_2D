subroutine err_gj(step, ni, nj, dx, dy, T, a, b, kappa, phi)
    implicit none
    integer :: ni, nj, step
    double precision :: dx, dy, dxinv, dyinv
    double precision :: T, a, b, kappa
    double precision :: phi(-6:ni+7, -6:nj+7)
    double precision :: phipx(-6:ni+7, -6:nj+7)
    double precision :: phipy(-6:ni+7, -6:nj+7)

    double precision :: psi, inv16, dphidx, dphidy
    integer :: i, j
    double precision :: x, y, GJ
    character(32) fname

    inv16 = 1.0d0/16.0d0
    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy

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
!$OMP& PRIVATE(i, j)
    do j = -5, nj+5
        do i = -6, ni+7
            phipy(i, j) = inv16* &
                    & (-phi(i, j-1)+9.0d0*phi(i, j)+9.0d0*phi(i, j+1)-phi(i, j+2))
            ! phipy(i, j)= 0.5d0*(phi(i, j)+phi(i, j+1))
        end do
    end do
!$OMP  END PARALLEL DO


    write (fname, '("err_gj",i7.7)') step
    open (11, file=fname)

    do j = -2, nj+3
        do i = -2, ni+3
            x = (dble(i)-0.5d0)*dx
            y = (dble(j)-0.5d0)*dy
            psi = phi(i, j)*(T*log(phi(i, j)/(1.0d0-b*phi(i, j)))-a*phi(i, j))
            call nabla(dxinv, dphidx, &
                     & phipx(i-2,j),  &
                     & phipx(i-1,j),  &
                     & phipx(i  ,j),  &
                     & phipx(i+1,j))
            call nabla(dyinv, dphidy, &
                     & phipy(i,j-2),  &
                     & phipy(i,j-1),  &
                     & phipy(i,j  ),  &
                     & phipy(i,j+1))
            GJ = -psi + 0.5d0*kappa*(dphidx**2+dphidy**2)
            write (11, *) x, y, GJ, -psi
        end do
        write (11, *)
    end do

    close(11)
end subroutine err_gj
