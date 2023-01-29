subroutine init(ni, nj, dx, dy, phi, phimin, phimax, R)
    implicit none

    integer :: ni, nj
    double precision :: dx, dy
    double precision :: phimin, phimax, R
    double precision :: phi(-6:ni+7, -6:nj+7)

    integer :: i, j
    double precision :: x, y, xo, yo
    integer :: width

    width = 10
    xo = (dble(ni/2))*dx
    yo = (dble(nj/2)+25.0d0)*dy

    ! write (*, *) midx, midy

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,x,y)
    do j = -6, nj+7
        do i = -6, ni+7
            x = (dble(i)-0.5d0)*dx
            y = (dble(j)-0.5d0)*dy
            if ((x-xo)**2+(y-yo)**2 <= R**2) then
                if (y < 80.0d0 .and. abs(x-xo) < 2.5d0) then
                    phi(i, j) = phimin
                else
                    phi(i, j) = phimax
                end if
            else
                phi(i, j) = phimin
            end if
        end do
    end do
!$OMP  END PARALLEL DO

! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i,j,x,y)
!     do j = -6, nj+7
!         do i = -6, ni+7
!             x = (dble(i)-0.5d0)*dx
!             y = (dble(j)-0.5d0)*dy
!             ! if ((x-midx)**2+(y-midy)**2 <= R**2 .and. (j >= nj/2)) then
!             if ((x-midx)**2+(y-midy)**2 <= R**2) then
!                 phi(i, j) = phimin
!             else
!                 phi(i, j) = phimax
!             end if
!         end do
!     end do
! !$OMP  END PARALLEL DO

end subroutine init
