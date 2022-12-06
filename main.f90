program main
    implicit none

    integer :: ni
    integer :: i
    double precision :: x
    double precision :: phimin, phimax
    double precision :: dx, dxinv
    double precision :: xl
    double precision :: t, dt
    double precision :: a, b, temperature, kappa
    double precision :: u
    double precision, dimension(:), allocatable :: phi
    integer :: maxstep, step
    integer :: dataou
    character(32) fname

    dataou = 1000
    maxstep = 80000
    ni = 128
    a = 1.0d0
    b = 1.0d0
    kappa = 0.1d0
    temperature = 0.293d0

    dt = 2.5e-2
    dx = 1.0d0

    phimin = 0.265d0
    phimax = 0.405d0

    xl = dx*dble(ni)

    u = 0.5d0
    step = 0
    include'allocate.h'
    dxinv = 1.0d0/dx

    write (*, '("Courant Number      ",20e20.10)') abs(u*dt/dx)
    call init(ni, phi, phimin, phimax, 32)

    include'mkphi.h'

    do step = 1, maxstep
        call bundset(ni, phi)
        ! call calphi(ni, u, dxinv, phi, dt, a, b, temperature, kappa)
        call calphi2(ni, u, dxinv, phi, dt, a, b, temperature, kappa)

        if (mod(step,dataou) == 0) then
            include'mkphi.h'
        end if

    end do

end program main
