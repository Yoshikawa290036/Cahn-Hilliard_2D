subroutine cal_vel_prediction(ni,nj,dt,rho,u,v,up,vp)
    implicit none
    integer :: ni,nj
    double precision :: dt
    double precision :: u(-6:ni+7, -6:nj+7)
    double precision :: v(-6:ni+7, -6:nj+7)
    double precision :: up(-6:ni+7, -6:nj+7)
    double precision :: vp(-6:ni+7, -6:nj+7)
    double precision :: rho(-6:ni+7, -6:nj+7)
    

end subroutine cal_vel_prediction
