program main
    implicit none

    integer :: ni, nj
    integer :: i, j, cnt, initcnt
    double precision :: x, y, sum, initsum
    double precision :: phimin, phimax, phimid
    double precision :: dx, dxinv, dy, dyinv
    double precision :: xl, yl
    double precision :: t, dt
    double precision :: a, b, temperature, kappa
    double precision, dimension(:, :), allocatable :: phi, u, v
    integer :: maxstep, step
    integer :: dataou
    character(32) fname

    dataou = 100
    maxstep = 10000

    ni = 64
    nj = 64

    a = 1.0d0
    b = 1.0d0
    kappa = 0.1d0
    temperature = 0.293d0

    dt = 2.5e-2
    dx = 1.0d0
    dy = 1.0d0

    phimin = 0.265d0
    phimax = 0.405d0
    phimid = 0.5d0*(phimin+phimax)

    xl = dx*dble(ni)
    yl = dy*dble(nj)
    dxinv = 1.0d0/dx
    dyinv = 1.0d0/dy

    include'allocate.h'
    ! write (*, '("Courant Number      ",20e20.10)') abs(u*dt/dx)
    call init(ni, nj, dx, dy, phi, phimin, phimax, 8.0d0)
    call bndset(ni, nj, phi)
    call cal_vel(ni, nj, u, v, xl, yl, dx, dy)
    step = 0
    include'mkphi.h'
    include'cal_erea.h'
    include'mkuvphi.h'
    call err_gj(step, ni, nj, dx, dy, temperature, a, b, kappa, phi)
    initsum = sum
    initcnt = cnt

    do step = 1, maxstep
        call calphi(ni, nj, u, v, dxinv, dyinv, phi, dt, a, b, temperature, kappa)

        call bndset(ni, nj, phi)
        if (mod(step, dataou) == 0) then
            include'mkphi.h'
            include'cal_erea.h'
            call err_gj(step, ni, nj, dx, dy, temperature, a, b, kappa, phi)
        end if
    end do

end program main
