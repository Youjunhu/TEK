module drift
contains
  pure subroutine compute_drift(ns, x, z, mu, vpar, touch_bdry,&
       & phix, phiy, phiz, ax, ay, az, &
       & xdrift0, zdrift0, ydrift0, mirror_force, zdrift00, xdrift1, zdrift1, ydrift1)
    !output are used by the tajectory pusher and the weight pusher
    use constants,only: p_, one,two,twopi,kev
    use table_in_mc,only: w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, bdgxcgz, bdgycgz, bdgxcgy, b_mc
    use magnetic_coordinates,only: mpol, nrad, zgrid, xgrid
    use interpolate_module,only: linear_2d_interpolate
    use gk_module,only: nm_gk, charge_sign_gk, vn_gk
    use normalizing,only: vu, bn, ln, qu,tu
    implicit none
    integer,  intent(in)  :: ns ! species index
    real(p_), intent(in)  :: x(:), z(:), mu(:), vpar(:)
    logical,  intent(in)  :: touch_bdry(:)
    real(p_), intent(in)  :: phix(:), phiy(:), phiz(:), ax(:), ay(:), az(:)
    real(p_), intent(out) :: xdrift0(:), ydrift0(:), zdrift0(:), zdrift00(:)
    real(p_), intent(out) :: mirror_force(:)
    real(p_), intent(out) :: xdrift1(:), ydrift1(:), zdrift1(:)
    real(p_) :: w1v, w2v, w3v, w4v, w5v, w6v, w7v, w8v, w9v, w10v
    real(p_) :: factor1, cs, norm1, norm2, bval, bsq, bdgxcgy0, bdgxcgz0, bdgycgz0
    integer :: k

    cs = charge_sign_gk(ns)
    norm1 = tu*kev/(qu*vn_gk(ns)*bn*ln)
    norm2 = tu*kev/(qu*vu*bn*ln)
    do k = 1, nm_gk(ns) !number of MC markers
       if( touch_bdry(k) .eqv. .true.) cycle
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w1,z(k),x(k),w1v) 
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w2,z(k),x(k),w2v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w3,z(k),x(k),w3v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w4,z(k),x(k),w4v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w5,z(k),x(k),w5v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w6,z(k),x(k),w6v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w7,z(k),x(k),w7v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w8,z(k),x(k),w8v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w9,z(k),x(k),w9v)
       !call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w10,z(k),x(k),w10v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bdgxcgz,z(k),x(k),bdgxcgz0)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bdgycgz,z(k),x(k),bdgycgz0)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bdgxcgy,z(k),x(k),bdgxcgy0)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,b_mc, z(k),x(k), bval)
       bsq = bval**2
       factor1=one+cs*vpar(k)/twopi*w1v
       !factor1=one !for testing, no difference
       xdrift0(k) = cs*vpar(k)**2/twopi*w3v/factor1 +cs*mu(k)/(twopi*factor1)*w6v 
       ydrift0(k) = cs*vpar(k)**2/(twopi*factor1)*w5v +cs*mu(k)/(twopi*factor1)*w8v
       zdrift0(k) = vpar(k)*(w2v+cs*vpar(k)/twopi*w4v)/factor1 &
            &  + cs*mu(k)/(twopi*factor1)*w7v
       !zdrift00(k) = zdrift0(k) - vpar(k)*w2v !removing the parallel streaming
       zdrift00(k) = vpar(k)*(cs*vpar(k)/twopi*w4v)/factor1 &
            &  + cs*mu(k)/(twopi*factor1)*w7v

       !mirror_force(k)=-mu(k)/factor1*(w9v+cs*vpar(k)/twopi*w10v)
       mirror_force(k) = - mu(k) * w9v !tested, agree with the result using the previous line
       !mirror_force2(k) = mu(k)/twopi*vpar(k)*w10v

       xdrift1(k) =  (-phiy(k)*bdgxcgy0 - phiz(k)*bdgxcgz0)/bsq*norm1 &
            &  + vpar(k)/bsq*(ay(k)*bdgxcgy0 + az(k)*bdgxcgz0)*norm2

       ydrift1(k)  = (phix(k)*bdgxcgy0 - phiz(k)*bdgycgz0)/bsq*norm1 & 
            &  + vpar(k)/bsq*(ax(k)*(-bdgxcgy0) + az(k)*bdgycgz0)*norm2

       zdrift1(k)  = (phix(k)*bdgxcgz0 + phiy(k)*bdgycgz0)/bsq*norm1 &
            &  + vpar(k)/bsq*(ax(k)*(-bdgxcgz0) + ay(k)*(-bdgycgz0))*norm2
    enddo

  end subroutine compute_drift
end module drift
