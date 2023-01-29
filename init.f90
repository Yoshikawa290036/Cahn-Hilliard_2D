subroutine init(ni, nj, dx, dy, phi, phimin, phimax, R)
    implicit none

    integer :: ni, nj
    double precision :: dx, dy
    double precision :: phimin, phimax, R
    double precision :: phi(-6:ni+7, -6:nj+7)

    integer :: i, j
    double precision :: x, y, midx, midy
    integer :: width

    width = 3
    midx = (dble(ni/2)-0.5d0)*dx
    midy = (dble(nj/2)-0.5d0)*dy

    ! write (*, *) midx, midy

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,x,y)

    ! do j = -6, nj+7
    !     do i = -6, ni+7
    !         x = (dble(i)-0.5d0)*dx
    !         y = (dble(j)-0.5d0)*dy
    !         if ((x-midx)**2+(y-midy)**2 <= R**2 .and. (j >= nj/2)) then
    !         ! if ((x-midx)**2+(y-midy)**2 <= R**2) then
    !             phi(i, j) = phimax
    !         else if (((x-midx)**2+(y-midy)**2 <= R**2) .and. (i >= ni/2+width .or. i <= ni/2-width)) then
    !             phi(i, j) = phimax
    !         else
    !             phi(i, j) = phimin
    !         end if
    !     end do
    ! end do

    ! air bubbles ------------
    do j = -6, nj+7
        do i = -6, ni+7
            x = (dble(i)-0.5d0)*dx
            y = (dble(j)-0.5d0)*dy
            if ((x-midx)**2+(y-midy)**2 <= R**2) then
                phi(i, j) = phimin
            else
                phi(i, j) = phimax
            end if
        end do
    end do

    ! droplet ----------------
    ! do j = -6, nj+7
    !     do i = -6, ni+7
    !         x = (dble(i)-0.5d0)*dx
    !         y = (dble(j)-0.5d0)*dy
    !         if ((x-midx)**2+(y-midy)**2 <= R**2) then
    !             phi(i, j) = phimax
    !         else
    !             phi(i, j) = phimin
    !         end if
    !     end do
    ! end do

    !$OMP  END PARALLEL DO

end subroutine init
