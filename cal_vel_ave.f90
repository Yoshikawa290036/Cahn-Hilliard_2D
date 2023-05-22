
! =====================================================

    ! ua(i,j) => average of u(i,j) at position v(i,j)
    ! va(i,j) => average of v(i,j) at position u(i,j)

! =====================================================



subroutine cal_vel_ave(ni, nj, u, v, ua, va)
    implicit none
    integer :: ni, nj
    double precision ::  u(-6:ni+7, -6:nj+7)
    double precision ::  v(-6:ni+7, -6:nj+7)
    double precision ::  ua(-6:ni+7, -6:nj+7)
    double precision ::  va(-6:ni+7, -6:nj+7)

    integer :: i, j



!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)

    do j = -6, nj+6
        do i = -5, ni+7
            ua(i, j) = 0.25d0*( + u(i-1,j  )  &
                                + u(i  ,j  )  &
                                + u(i-1,j+1)  &
                                + u(i  ,j+1))
        end do
    end do

!$OMP  END PARALLEL DO


!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, j)

    do j = -5, nj+7
        do i = -6, ni+6
            va(i, j) = 0.25d0*( + v(i  ,j-1)  &
                                + v(i+1,j-1)  &
                                + v(i  ,j  )  &
                                + v(i+1,j  ))
        end do
    end do

!$OMP  END PARALLEL DO



end subroutine cal_vel_ave
