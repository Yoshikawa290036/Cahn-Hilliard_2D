sum = 0.0d0
cnt = 0
do j = 1, nj-1
    do i = 1, ni-1
        sum = sum + phi(i,j)
        if ( phi(i,j) < phimid ) then
            cnt = cnt + 1
        end if
    end do
end do
write (*, *) step, dble(step)*dt, abs(sum/initsum-1.0d0), dble(cnt)/dble(initcnt)
