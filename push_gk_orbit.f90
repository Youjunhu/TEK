module gk_trajectory_pusher
  use constants,only:p_
  use constants,only: pi,one_half
  implicit none
contains

  subroutine push_gc_half_step(ns,xdrift0, zdrift0, ydrift0, mirror_force,&
       & xdrift1, zdrift1, ydrift1) !Euler step (to estimate the values at the half-time-step)
    use gk_module,only: dtao_e, nm_gk, radcor_e, theta_e, alpha_e, vpar_e, mu_e, &
         & touch_bdry_e, & !input_and_output
         & radcor_e_mid, theta_e_mid, alpha_e_mid, vpar_e_mid !output
    use math,only: shift_to_specified_toroidal_range
    use control_parameters,only: gk_nonlinear
    integer, intent(in) :: ns
    real(p_), intent(in) :: xdrift0(:), zdrift0(:), ydrift0(:), mirror_force(:)
    real(p_), intent(in) :: xdrift1(:), zdrift1(:), ydrift1(:)
    real(p_) :: xdrift, zdrift, ydrift
    integer:: nm, k

    nm = nm_gk(ns)
    !    !$omp parallel do
    do k=1,nm
       if( touch_bdry_e(k,ns).eqv..true.) cycle
       if(gk_nonlinear==1) then
          xdrift = xdrift0(k) + xdrift1(k)
          zdrift = zdrift0(k) + zdrift1(k)
          ydrift = ydrift0(k) + ydrift1(k)
       else
          xdrift = xdrift0(k)
          zdrift = zdrift0(k)
          ydrift = ydrift0(k)
       endif

       radcor_e_mid(k,ns)=radcor_e(k,ns)+ xdrift*dtao_e(ns)*one_half
       theta_e_mid(k,ns) =theta_e(k,ns) + zdrift*dtao_e(ns)*one_half
       alpha_e_mid(k,ns) =alpha_e(k,ns) + ydrift*dtao_e(ns)*one_half
       vpar_e_mid(k,ns)  =vpar_e(k,ns)  + mirror_force(k)*dtao_e(ns)*one_half

       if ((theta_e_mid(k,ns) .ge. pi) .or. (theta_e_mid(k,ns)<-pi)) &
            & call shift_gc_theta_then_alpha(radcor_e_mid(k,ns),theta_e_mid(k,ns),alpha_e_mid(k,ns))
       call shift_to_specified_toroidal_range(alpha_e_mid(k,ns))
    enddo
    !    !$omp end parallel do
    call gc_radial_bdry_condition(nm,radcor_e_mid(:,ns), touch_bdry_e(:,ns))
  end subroutine push_gc_half_step

  subroutine push_gc_full_step(ns,xdrift0, zdrift0, ydrift0, mirror_force, xdrift1, zdrift1, ydrift1)
    !using values of (r,v) at t{n+1/2} to compute the rhs of motion equation to push (r,v) from t_{n} to t_{n+1)
    use gk_module,only: dtao_e, nm_gk, & !input
         & radcor_e, theta_e, alpha_e, vpar_e, & !input and output
         & touch_bdry_e !output
    use control_parameters,only: gk_nonlinear
    use interpolate_module
        use magnetic_coordinates,only: toroidal_range 
    use math,only: shift_to_specified_toroidal_range
    implicit none
    integer, intent(in) :: ns
    real(p_), intent(in) :: xdrift0(:), zdrift0(:), ydrift0(:), mirror_force(:)
    real(p_), intent(in) :: xdrift1(:), zdrift1(:), ydrift1(:)
    real(p_) :: xdrift, zdrift, ydrift
    integer:: nm, k

    nm=nm_gk(ns)
    !    !$omp parallel do
    do k=1,nm
       if((touch_bdry_e(k,ns).eqv..true.)) cycle
       if(gk_nonlinear==1) then
          xdrift = xdrift0(k) + xdrift1(k)
          zdrift = zdrift0(k) + zdrift1(k)
          ydrift = ydrift0(k) + ydrift1(k)
       else
          xdrift = xdrift0(k)
          zdrift = zdrift0(k)
          ydrift = ydrift0(k)
       endif
       radcor_e(k,ns)=radcor_e(k,ns)+ xdrift*dtao_e(ns)
       theta_e(k,ns) =theta_e(k,ns) + zdrift*dtao_e(ns)
       alpha_e(k,ns)=alpha_e(k,ns)  + ydrift*dtao_e(ns)
       vpar_e(k,ns) =vpar_e(k,ns)   + mirror_force(k)*dtao_e(ns)

       if(theta_e(k,ns) .ge. pi .or. theta_e(k,ns)<-pi) &
            & call shift_gc_theta_then_alpha(radcor_e(k,ns),theta_e(k,ns),alpha_e(k,ns))
       
       call shift_to_specified_toroidal_range(alpha_e(k,ns))
    enddo

!!$        do k=1,nm_gk(ns)
!!$           if(touch_bdry_e(k,ns).eqv. .true.) cycle
!!$           if(alpha_e(k,ns)>toroidal_range) print *, 'exceed upper******, ', alpha_e(k,ns)
!!$           if(alpha_e(k,ns)<0) print *, 'exceed lower*****in push*, ',  alpha_e(k,ns)
!!$        enddo
    
    !    !$omp end parallel do
    call gc_radial_bdry_condition(nm,radcor_e(:,ns),touch_bdry_e(:,ns))
  end subroutine push_gc_full_step

pure  subroutine shift_gc_theta_then_alpha(radcor,theta,alpha)
    !shift theta, and then shift alpha so that phi is not changed, i.e., keep the particle in the same sptial location
    use constants,only:p_
    use constants,only: twopi,pi
    use magnetic_field, only : qfunc
    implicit none
    real(p_),intent(in):: radcor
    real(p_),intent(inout):: theta,alpha
    integer:: ishift

    ishift=floor(theta/twopi)
    theta=theta-ishift*twopi
    if(theta>pi) then
       theta=theta-twopi
       ishift=ishift+1
    endif

    alpha=alpha+ishift*twopi*qfunc(radcor) !this shift is needed to make the cylindrical toroidal angle phi be the same as the value before the theta-shifting, i.e., make the particle lie at the same sptaial location when doing the theta-shifting
  end subroutine shift_gc_theta_then_alpha


pure  subroutine gc_radial_bdry_condition(nm,radcor,touch_bdry)
    use constants,only:p_
    use magnetic_coordinates,only: radcor_low2,radcor_upp2
    implicit none
    integer,intent(in) :: nm
    real(p_),intent(in) ::  radcor(:)
    logical,intent(out) :: touch_bdry(:)
    integer :: k

    do k=1,nm
       if((radcor(k).le.radcor_upp2) .and. (radcor(k).ge.radcor_low2)) then
          touch_bdry(k)=.false.
       else
          touch_bdry(k)=.true.
       endif
    enddo
  end subroutine gc_radial_bdry_condition


  subroutine count_lost_markers_gk(ns)  !diagnosis
    use gk_module, only : nm_gk, touch_bdry_e
    use domain_decomposition,only: myid
    implicit none
    integer, intent(in) :: ns
    integer:: k, nlost
    nlost=0
    do k=1,nm_gk(ns)
       if(touch_bdry_e(k,ns).eqv..true.) then
          nlost=nlost+1
       endif
    enddo
    write(*,*) 'number of lost gk markers=', nlost, 'myid=',myid, 'species=',ns, ',total number=', nm_gk(ns)
  end subroutine count_lost_markers_gk

end module gk_trajectory_pusher
