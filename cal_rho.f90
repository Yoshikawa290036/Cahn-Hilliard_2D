subroutine cal_rho(ni, nj, rhoL, rhoG, phimin, phimax, phi, rho)
    implicit none
    integer :: ni, nj
    double precision :: rhoL, rhoG, phimin, phimax
    double precision :: phi(-6:ni+7, -6:nj+7)
    double precision :: rho(-6:ni+7, -6:nj+7)
    integer :: i, j
    double precision :: phimid, rhomid, rhomrg, PI

    PI = atan(1.0d0)*4.0d0
    phimid = 0.5d0*(phimin+phimax)
    rhomid = 0.5d0*(rhoL+rhoG)
    rhomrg = 0.5d0*(rhoL-rhoG)

    do j = -6, nj+7
        do i = -6, ni+7
            rho(i, j) = rhomid+rhomrg*sin(PI*(phimid-phi(i,j))/(phimax-phimin))
        end do
    end do
end subroutine cal_rho
