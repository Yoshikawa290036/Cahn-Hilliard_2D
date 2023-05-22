subroutine cal_adv(ni, nj, u, v, ua, va, advx, advy, dxinv, dyinv)
    implicit none
    integer :: ni, nj
    double precision :: u(-6:ni+7, -6:nj+7)
    double precision :: v(-6:ni+7, -6:nj+7)
    double precision :: ua(-6:ni+7, -6:nj+7)
    double precision :: va(-6:ni+7, -6:nj+7)
    double precision :: advx(-6:ni+7, -6:nj+7)
    double precision :: advy(-6:ni+7, -6:nj+7)

    double precision :: dxinv
    double precision :: dyinv
    double precision :: duudx,duvdx,duvdy,dvvdy

    integer :: i, j
    double precision :: inv12, inv24, x, y

    inv12 = 1.0d0/12.0d0
    inv24 = 1.0d0/24.0d0


! =====================================================

    ! ua(i,j) => average of u(i,j) at position v(i,j)
    ! va(i,j) => average of v(i,j) at position u(i,j)

! =====================================================


    ! do j = -6, nj+7
    !     do i = -6, ni+7
    !         ! cal u ================
    !         x = (dble(i)-0.5d0)/dxinv
    !         y = (dble(j)      )/dyinv
    !         u(i,j) = sin(x*0.1d0)

    !         ! cal v ================
    !         x = (dble(i)      )/dxinv
    !         y = (dble(j)-0.5d0)/dyinv
    !         v(i,j) = cos(y*0.1d0)
    !     end do
    ! end do


! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i, j)

!     do j = -6, nj+6
!         do i = -5, ni+7
!             ua(i, j) = 0.25d0*( + u(i-1,j  )  &
!                                 + u(i  ,j  )  &
!                                 + u(i-1,j+1)  &
!                                 + u(i  ,j+1))
!         end do
!     end do

! !$OMP  END PARALLEL DO


! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i, j)

!     do j = -5, nj+7
!         do i = -6, ni+6
!             va(i, j) = 0.25d0*( + v(i  ,j-1)  &
!                                 + v(i+1,j-1)  &
!                                 + v(i  ,j  )  &
!                                 + v(i+1,j  ))
!         end do
!     end do

! !$OMP  END PARALLEL DO




!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j, duudx, duvdy)

    do j = -4, nj+5
        do i = -4, ni+5

            call nabla(dxinv, duudx, inv24,      &
                         u(i-2, j) * u(i-2, j),  &
                         u(i-1, j) * u(i-1, j),  &
                         u(i  , j) * u(i  , j),  &
                         u(i+1, j) * u(i+1, j))
            
            call nabla(dyinv, duvdy, inv24,      &
                        ua(i, j-2) * v(i, j-2),  &
                        ua(i, j-1) * v(i, j-1),  &
                        ua(i, j  ) * v(i, j  ),  &
                        ua(i, j+1) * v(i, j+1))
            
            advx(i,j) = duudx + duvdy
        end do
    end do
!$OMP  END PARALLEL DO


!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j, duvdx, dvvdy)

    do j = -4, nj+5
        do i = -4, ni+5

            call nabla(dxinv, duvdx, inv24,      &
                        u(i-2, j) * va(i-2, j),  &
                        u(i-1, j) * va(i-1, j),  &
                        u(i  , j) * va(i  , j),  &
                        u(i+1, j) * va(i+1, j))
            
            call nabla(dyinv, dvvdy, inv24,      &
                        v(i, j-2) *  v(i, j-2),  &
                        v(i, j-1) *  v(i, j-1),  &
                        v(i, j  ) *  v(i, j  ),  &
                        v(i, j+1) *  v(i, j+1))
            
            advy(i,j) = duvdx + dvvdy
        end do
    end do

!$OMP  END PARALLEL DO


! do j = 1, nj
!     do i = 1, ni
!         x = (dble(i)-0.5d0)/dxinv
!         y = (dble(j)-0.5d0)/dyinv
!         write (33, '(20e20.10)') x, y, advx(i,j),advy(i,j)
!     end do
!     write (33, *)
! end do


end subroutine cal_adv
