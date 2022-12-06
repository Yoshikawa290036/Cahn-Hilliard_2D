subroutine init(ni, phi, phimin, phimax, width)
    implicit none

    integer :: ni, width
    double precision :: phimin, phimax
    double precision :: phi(-6:ni+7)

    integer :: i

    do i = -3, ni+4
        if (i < ni/2-width/2 .or. i > ni/2+width/2) then
            phi(i) = phimax
        else
            phi(i) = phimin
        end if
    end do
end subroutine init
