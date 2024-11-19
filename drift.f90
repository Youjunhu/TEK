module drift
contains
  subroutine compute_drift(ns, radcor_e, theta_e, mu_e, vpar_e, touch_bdry_e,&
       & phix_e, phiy_e, phiz_e, ax_e, ay_e, az_e, &
       & xdrift0, zdrift0, ydrift0, mirror_force, zdrift0m, xdrift1, zdrift1, ydrift1)
    !output are used by both the tajectory pusher and the weight pusher
    use constants,only: p_, one,two,twopi,kev
    use table_in_mc,only: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w14,w15, bdgxcgy, b_mc
    use magnetic_coordinates,only: mpol,nflux,theta_1d_array,radcor_1d_array
    use interpolate_module,only: linear_2d_interpolation
    use gk_module,only: nm_gk, charge_sign_e, vn_e
    use normalizing,only: vu, bn, ln, qu,tu
    implicit none
    integer, intent(in)  :: ns !ns is species index
    real(p_), intent(in) :: radcor_e(:), theta_e(:), mu_e(:), vpar_e(:)
    logical, intent(in)  :: touch_bdry_e(:)
    real(p_), intent(in) :: phix_e(:), phiy_e(:), phiz_e(:), ax_e(:), ay_e(:), az_e(:)
    real(p_), intent(out) :: xdrift0(:), zdrift0(:), zdrift0m(:), ydrift0(:), mirror_force(:)
    real(p_), intent(out) :: xdrift1(:), zdrift1(:), ydrift1(:)
    real(p_) :: w1val, w2val, w3val, w4val, w5val,w6val,w7val,w8val,w9val,w10val,w14val, w15val
    real(p_) :: bpar_star, factor1,cs, norm, norm2, bval0, bdgxcgy0
    integer :: k

!!$omp parallel do private(w1val,w2val,w3val,w4val,w5val,w6val,w7val,w8val,w9val,w10val,&
!!$omp & factor1, bpar_star, x,y, epar,ex,ey, gyro_radius,dpsi1,dalpha1,dalpha2)
    cs=charge_sign_e(ns)
    norm=tu*kev/(qu*vn_e(ns)*bn*ln)
    norm2=norm*vn_e(ns)/vu

    do k=1, nm_gk(ns) !number of MC markers
       if( touch_bdry_e(k).eqv..true.) cycle
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w1,theta_e(k),radcor_e(k),w1val) 
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w2,theta_e(k),radcor_e(k),w2val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w3,theta_e(k),radcor_e(k),w3val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w4,theta_e(k),radcor_e(k),w4val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w5,theta_e(k),radcor_e(k),w5val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w6,theta_e(k),radcor_e(k),w6val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w7,theta_e(k),radcor_e(k),w7val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w8,theta_e(k),radcor_e(k),w8val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w9,theta_e(k),radcor_e(k),w9val)
       !call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w10,theta_e(k),radcor_e(k),w10val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w14,theta_e(k),radcor_e(k),w14val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,w15,theta_e(k),radcor_e(k),w15val)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,bdgxcgy,theta_e(k),radcor_e(k),bdgxcgy0)
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,b_mc, theta_e(k),radcor_e(k),bval0)

       bpar_star=bval0*(one+cs*vpar_e(k)/twopi*w1val)
       factor1=one+cs*vpar_e(k)/twopi*w1val
       !factor1=one !for testing, no difference
       xdrift0(k)=cs*vpar_e(k)**2/twopi*w3val/factor1 +cs*mu_e(k)/(twopi*factor1)*w6val 
       zdrift0(k)=vpar_e(k)*(w2val+cs*vpar_e(k)/twopi*w4val)/factor1 &
            &  + cs*mu_e(k)/(twopi*factor1)*w7val
       !zdrift0m(k) = zdrift0(k) - vpar_e(k)*w2val !removing the parallel streaming
       zdrift0m(k)=vpar_e(k)*(cs*vpar_e(k)/twopi*w4val)/factor1 &
            &  + cs*mu_e(k)/(twopi*factor1)*w7val
       ydrift0(k)=cs*vpar_e(k)**2/(twopi*factor1)*w5val +cs*mu_e(k)/(twopi*factor1)*w8val 
       !mirror_force(k)=-mu_e(k)/factor1*(w9val+cs*vpar_e(k)/twopi*w10val)
       mirror_force(k)=-mu_e(k)*w9val !tested, agree with the result using the previous line
       !xdrift1(k) =  ey_e(k)/GSpsi_prime*vu/vn_e(ns) & !/factor1 !including or not factors gives the same ITG spectrum
       xdrift1(k) =  (-phiy_e(k)/bval0**2*bdgxcgy0 - phiz_e(k)/bval0**2*w14val)*norm &
            &  +vpar_e(k)/bval0**2*(ay_e(k)*bdgxcgy0+az_e(k)*w14val)*norm2
       !ydrift1(k)  = -ex_e(k)/GSpsi_prime*vu/vn_e(ns) & !/factor1
       ydrift1(k)  = (phix_e(k)/bval0**2*bdgxcgy0 - phiz_e(k)/bval0**2*w15val)*norm & 
            &  +vpar_e(k)/bval0**2*(ax_e(k)*(-bdgxcgy0)+az_e(k)*w15val)*norm2
       !zdrift1(k)  = -(ex_e(k)*w14val+ey_e(k)*w15val)/(bval0**2)*vu/vn_e(ns) &
       zdrift1(k)  = (phix_e(k)/bval0**2*w14val + phiy_e(k)/bval0**2*w15val)*norm &
            &  +vpar_e(k)/bval0**2*(ax_e(k)*(-w14val)+ay_e(k)*(-w15val))*norm2
!!$         if(isNan(zdrift1(k))) then
!!$          write(*,*) '***in computing_drift***k=', k, ex_e(k), ey_e(k), bval0, vn_e(ns)
!!$          stop 'theta_drift is nan'
!!$          endif
    enddo
!!$omp end parallel do
  end subroutine compute_drift
end module drift
