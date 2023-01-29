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
    double precision :: R
    integer :: maxstep, step
    integer :: dataou, hoge
    character(32) fname

    dataou = 314
    maxstep = 62800
    ! maxstep = 1

    ni = 100
    nj = 100
    R = 15.0d0
    a = 1.0d0
    b = 1.0d0
    kappa = 0.1d0
    temperature = 0.293d0

    dt = 0.01d0
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
    call init(ni, nj, dx, dy, phi, phimin, phimax, R)
    call bndset(ni, nj, phi)
    call cal_vel(ni, nj, u, v, xl, yl, dx, dy)
    ! write (*, '("Courant Number      ",20e20.10)') abs(u(0, nj/2+R*dyinv)*dt/dx)
    step = 0
    include'mkphi.h'
    include'cal_erea.h'
    include'mkuvphi.h'
    initsum = sum
    initcnt = cnt

    do step = 0, maxstep
        call calphi(ni, nj, u, v, dxinv, dyinv, phi, dt, a, b, temperature, kappa)

        call bndset(ni, nj, phi)
        if (mod(step, dataou) == 0) then
            include'mkphi.h'
            include'cal_erea.h'
        end if
    end do

end program main
