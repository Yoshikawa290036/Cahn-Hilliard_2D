subroutine calphi(ni, u, dxinv, phi, dt, a, b, temperature, kappa)
    implicit none
    integer :: ni
    double precision :: u
    double precision :: phi(-6:ni+7), phip(-6:ni+7), dphi, dphip, ddphi(-6:ni+7)
    double precision :: nphi1(-6:ni+7), nphi2(-6:ni+7)
    double precision :: dt
    double precision :: a, b, temperature, kappa
    double precision :: dxinv
    double precision :: inv3
    integer :: i
    double precision :: adv(-6:ni+7), chp(-6:ni+7)
    double precision :: Jp(-6:ni+7), zetap
    double precision :: Gamma
    double precision :: m1, m2, p0, p1, p2
    double precision :: dddphip, dchp

    Gamma = 12.0d0

    do i = -5, ni+5
        phip(i) = 1.0d0/16.0d0* &
                & (-phi(i-1)+9.0d0*phi(i)+9.0d0*phi(i+1)-phi(i+2))
    end do

    ! ddphi = nabla nabla phi
    do i = -4, ni+5
        m2 = -1.0d0/24.0d0*dxinv**2*phi(i-2)
        m1 =  1.0d0/24.0d0*dxinv**2*phi(i-1)*16.0d0
        p0 = -1.0d0/24.0d0*dxinv**2*phi(i  )*30.0d0
        p1 =  1.0d0/24.0d0*dxinv**2*phi(i+1)*16.0d0
        p2 = -1.0d0/24.0d0*dxinv**2*phi(i+2)
        ddphi(i) = m2 + m1 + p0 + p1 + p2
    end do

    do i = -2, ni+3
        zetap = (temperature/((1.0d0-b*phip(i))**2))-2.0d0*a*phip(i)
        call nabla(dxinv, phi(i-1), phi(i), phi(i+1), phi(i+2), dphip)
        call nabla(dxinv, ddphi(i-1), ddphi(i), ddphi(i+1), ddphi(i+2), dddphip)

        ! Jp(i) = -zetap*dphip + kappa*phip(i)*dddphip
        Jp(i) = -zetap*dphip + kappa*phip(i)*dddphip

    end do

    do i = 0, ni+2
        call nabla(dxinv, phip(i-2), phip(i-1), phip(i), phip(i+1), dphi)
        adv(i) = u*dphi

        call nabla(dxinv, phip(i-2)*Jp(i-2), phip(i-1)*Jp(i-1), phip(i)*Jp(i),phip(i+1)*Jp(i+1), dchp)
        chp(i) = gamma*dchp

        nphi1(i) = phi(i) - dt*(adv(i)+chp(i))
        ! nphi1(i) = phi(i) - dt*(chp(i))

    end do

    do i = 0, ni+1
        phi(i) = nphi1(i)
    end do
end subroutine calphi
