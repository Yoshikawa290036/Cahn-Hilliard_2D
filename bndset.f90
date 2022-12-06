
! periodic boundary conditions

subroutine bundset(ni, nj, phi)
    implicit none

    integer :: ni, nj
    double precision :: phi(-6:ni+7, -6:nj+7)
    integer :: i, j

    ! left wall
    do j = 1, nj
        phi(0, j) = phi(ni, j)
        phi(-1, j) = phi(ni-1, j)
        phi(-2, j) = phi(ni-2, j)
        phi(-3, j) = phi(ni-3, j)
        phi(-4, j) = phi(ni-4, j)
        phi(-5, j) = phi(ni-5, j)
        phi(-6, j) = phi(ni-6, j)
    end do

    ! right wall
    do j = 1, nj
        phi(ni+1, j) = phi(1, j)
        phi(ni+2, j) = phi(2, j)
        phi(ni+3, j) = phi(3, j)
        phi(ni+4, j) = phi(4, j)
        phi(ni+5, j) = phi(5, j)
        phi(ni+6, j) = phi(6, j)
        phi(ni+7, j) = phi(7, j)
    end do

    ! bottom wall
    do i = 0, ni
        phi(i, 0) = phi(i, nj)
        phi(i, -1) = phi(i, nj-1)
        phi(i, -2) = phi(i, nj-2)
        phi(i, -3) = phi(i, nj-3)
        phi(i, -4) = phi(i, nj-4)
        phi(i, -5) = phi(i, nj-5)
        phi(i, -6) = phi(i, nj-6)
    end do

    ! top wall
    do i = 0, ni
        phi(i, nj+1) = phi(i, 1)
        phi(i, nj+2) = phi(i, 2)
        phi(i, nj+3) = phi(i, 3)
        phi(i, nj+4) = phi(i, 4)
        phi(i, nj+5) = phi(i, 5)
        phi(i, nj+6) = phi(i, 6)
        phi(i, nj+7) = phi(i, 7)
    end do

end subroutine bundset
