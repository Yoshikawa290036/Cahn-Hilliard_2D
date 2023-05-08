write (fname, '("xyrho",i7.7)') step

open (13, file=fname)

do j = -6, nj+7
    do i = -6, ni+7
        x = (dble(i)-0.5d0)*dx
        y = (dble(j)-0.5d0)*dy
        write (13, '(20e20.10)') x, y, rho(i,j)
    end do
    write (13, *)
end do
