
subroutine bndset_phi(ni, nj, phi)
    implicit none
    
    integer :: ni, nj
    double precision :: phi(-6:ni+7, -6:nj+7)
    integer :: i, j
    double precision :: setv


! left wall =================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j, setv)
    do j = 1, nj
        setv = phi(ni, j)

        phi(ni+1, j) = setv
        phi(ni+2, j) = setv
        phi(ni+3, j) = setv
        phi(ni+4, j) = setv
        phi(ni+5, j) = setv
        phi(ni+6, j) = setv
        phi(ni+7, j) = setv
    end do
!$OMP  END PARALLEL DO

! right wall =================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(j, setv)
    do j = 1, nj
        setv = phi(1, j)

        phi( 0, j) = setv
        phi(-1, j) = setv
        phi(-2, j) = setv
        phi(-3, j) = setv
        phi(-4, j) = setv
        phi(-5, j) = setv
        phi(-6, j) = setv
    end do
!$OMP  END PARALLEL DO

! bottom wall =================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, setv)
    do i = 0, ni
        setv = phi(i, 1)

        phi(i,  0) = setv
        phi(i, -1) = setv
        phi(i, -2) = setv
        phi(i, -3) = setv
        phi(i, -4) = setv
        phi(i, -5) = setv
        phi(i, -6) = setv
    end do
!$OMP  END PARALLEL DO

! top wall =================================
!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i, setv)
    do i = 0, ni
        setv = phi(i, nj)

        phi(i, nj+1) = setv
        phi(i, nj+2) = setv
        phi(i, nj+3) = setv
        phi(i, nj+4) = setv
        phi(i, nj+5) = setv
        phi(i, nj+6) = setv
        phi(i, nj+7) = setv
    end do
!$OMP  END PARALLEL DO

end subroutine bndset_phi
