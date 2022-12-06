write (fname, '("xphi",i7.7)') step

open (10, file=fname)

do i = 0, ni+1
    x = (dble(i)-0.5d0)*dx
    write (10, '(20e20.10)') x, phi(i)
end do
