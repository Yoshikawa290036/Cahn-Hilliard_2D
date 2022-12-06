
! periodic boundary conditions

subroutine bundset(ni, phi)
    implicit none

    integer :: ni
    double precision :: phi(-6:ni+7)

    phi(0) = phi(ni)
    phi(-1) = phi(ni-1)
    phi(-2) = phi(ni-2)
    phi(-3) = phi(ni-3)
    phi(-4) = phi(ni-4)
    phi(-5) = phi(ni-5)
    phi(-6) = phi(ni-6)

    phi(ni+1) = phi(1)
    phi(ni+2) = phi(2)
    phi(ni+3) = phi(3)
    phi(ni+4) = phi(4)
    phi(ni+5) = phi(5)
    phi(ni+6) = phi(6)
    phi(ni+7) = phi(7)

end subroutine bundset
