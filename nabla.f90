! ans(i+1/2)

subroutine nabla(dxinv, ans, inv24, m2, m1, p1, p2)
    implicit none
    double precision :: dxinv, m2, m1, p1, p2, ans
    double precision :: mm2, mm1, pp1, pp2
    double precision :: inv24

    ! inv24 = 1.0d0/24.0d0

    mm2 = inv24*dxinv*m2
    mm1 = inv24*dxinv*m1*27.0d0
    pp1 = inv24*dxinv*p1*27.0d0
    pp2 = inv24*dxinv*p2

    ans = mm2 - mm1 + pp1 - pp2
end subroutine nabla
