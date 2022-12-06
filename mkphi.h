write (fname, '("xyphi",i7.7)') step

open (10, file=fname)

do j = 0, nj+1
    do i = 0, ni+1
        x = (dble(i)-0.5d0)*dx
        y = (dble(j)-0.5d0)*dy
        write (10, '(20e20.10)') x, y, phi(i,j)
    end do
    write (10, *)
end do
