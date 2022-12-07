subroutine init(ni, nj, dx, dy, phi, phimin, phimax, R)
    implicit none

    integer :: ni, nj
    double precision :: dx, dy
    double precision :: phimin, phimax, R
    double precision :: phi(-6:ni+7, -6:nj+7)

    integer :: i, j
    double precision :: x, y, midx, midy

    midx = (dble(ni/2)-0.5d0)*dx
    midy = (dble(nj/2)-0.5d0)*dy

    ! write (*, *) midx, midy
    do j = -6, nj+7
        do i = -6, ni+7
            x = (dble(i)-0.5d0)*dx
            y = (dble(j)-0.5d0)*dy
            if ((x-midx)**2+(y-midy)**2 <= R**2) then
                phi(i, j) = phimax
            else
                phi(i, j) = phimin
            end if
        end do
    end do
end subroutine init
