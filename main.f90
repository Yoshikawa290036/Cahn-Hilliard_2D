program main
    implicit none

    integer :: ni, nj
    integer :: i, j, cnt, initcnt
    double precision :: x, y, sum, initsum
    double precision :: phimin, phimax, phimid
    double precision :: dx, dxinv, dy, dyinv
    double precision :: xl, yl
    double precision :: t, dt
    double precision :: a, b, temperature, kappa_phi, kappa_s
    double precision :: rhoL, rhoG
    double precision :: muL, muG
    double precision, dimension(:, :), allocatable :: phi, u, v
    double precision, dimension(:, :), allocatable :: ua, va
    double precision, dimension(:, :), allocatable :: rho, up, vp
    double precision, dimension(:, :), allocatable :: eta, Fsx, Fsy
    double precision, dimension(:, :), allocatable :: advx, advy
    double precision, dimension(:, :), allocatable :: visx, visy
    double precision :: R
    integer :: maxstep, step
    integer :: dataou, hoge
    character(32) fname

    dataou = 314
    maxstep = 1
    ! maxstep = dataou*5
    
    ni = 64
    nj = 64
    R = 15.0d0
    a = 1.0d0
    b = 1.0d0

    kappa_phi = 0.1d0
    kappa_s = 1.71e3
    rhoL = 1.0e-3
    rhoG = 1.25e-6
    muL = 1.06e-5
    muG = 1.43e-7

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
    call init(ni, nj, dx, dy, phi, phimin, phimax, R)
    call bndset(ni, nj, phi)
    call cal_vel(ni, nj, u, v, xl, yl, dx, dy)
    step = 0
    include'mkphi.h'
    include'cal_erea.h'
    include'mkuvphi.h'
    initsum = sum
    initcnt = cnt

    do step = 0, maxstep
        call cal_phi(ni, nj, u, v, eta, dxinv, dyinv, phi, dt, a, b, temperature, kappa_phi)
        call bndset(ni, nj, phi)        
        call cal_rho(ni, nj, rhoL, rhoG, phimin, phimax, phi, rho)
        call cal_Fs(ni, nj, dxinv, dyinv, rho, Fsx, Fsy, kappa_s)
        call cal_vel_ave(ni, nj, u, v, ua, va)
        call cal_adv(ni, nj, u, v, ua, va, advx, advy, dxinv, dyinv)
        call cal_vis(ni, nj, visx, visy, u, v, ua, va, &
                     dxinv, dyinv, muL, muG, rhoL, rhoG, rho)
        call cal_velp(ni, nj, dt, rho, up, vp, &
                    advx, advy, Fsx, Fsy, visx, visy)
        
        if (mod(step, dataou) == 0) then
            include'mkphi.h'
            include'mkrho.h'
            include'cal_erea.h'
        end if
        flush(6)
        
    end do
end program main
