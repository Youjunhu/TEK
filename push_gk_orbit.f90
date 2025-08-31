module gk_trajectory_pusher
contains

  pure subroutine push_gc(ns, dt, xdrift0, zdrift0, ydrift0, mirror_force,&
       & xdrift1, zdrift1, ydrift1, x_old, z_old, y_old, vpar_old, x_new, z_new, y_new, vpar_new, touch_bdry)
    use constants, only: p_, pi, one_half
    use gk_module, only: nm_gk, gk_nonlinear
    use math, only: shift_toroidal
    use magnetic_coordinates, only : toroidal_range
    implicit none
    integer, intent(in) :: ns
    real(p_), intent(in) :: dt
    real(p_), intent(in) :: xdrift0(:), zdrift0(:), ydrift0(:), mirror_force(:)
    real(p_), intent(in) :: xdrift1(:), zdrift1(:), ydrift1(:)
    real(p_), intent(in) :: x_old(:), y_old(:), z_old(:), vpar_old(:)
    real(p_), intent(out) :: x_new(:), y_new(:), z_new(:), vpar_new(:)
    logical, intent(inout) :: touch_bdry(:)
    real(p_) :: xdrift, zdrift, ydrift
    integer :: nm, k

    nm = nm_gk(ns)

    do k =1, nm
       if( touch_bdry(k).eqv..true.) cycle
       if(gk_nonlinear==1) then
          xdrift = xdrift0(k) + xdrift1(k)
          zdrift = zdrift0(k) + zdrift1(k)
          ydrift = ydrift0(k) + ydrift1(k)
       else
          xdrift = xdrift0(k)
          zdrift = zdrift0(k)
          ydrift = ydrift0(k)
       endif

       x_new(k) = x_old(k)+ xdrift*dt
       z_new(k) = z_old(k)+ zdrift*dt
       y_new(k) = y_old(k)+ ydrift*dt
       vpar_new(k)   = vpar_old(k)  + mirror_force(k)*dt

       if ((z_new(k) .ge. pi) .or. (z_new(k) < -pi)) &
            & call shift_gc_theta_then_alpha(x_new(k), z_new(k), y_new(k))
       call shift_toroidal(y_new(k),toroidal_range)
    enddo

    call gc_radial_bdry_condition(nm, x_new(:), touch_bdry(:))
  end subroutine push_gc

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
    use constants, only:p_
    use magnetic_coordinates, only: xlow,xupp
    implicit none
    integer, intent(in) :: nm
    real(p_), intent(in) ::  radcor(:)
    logical, intent(out) :: touch_bdry(:)
    integer :: k

    do k=1,nm
       if((radcor(k).le.xupp) .and. (radcor(k).ge.xlow)) then
          touch_bdry(k)=.false.
       else
          touch_bdry(k)=.true.
       endif
    enddo
  end subroutine gc_radial_bdry_condition


  subroutine count_lost_markers_gk(ns)  !diagnosis
    use gk_module, only : nm_gk, touch_bdry_gc
    use domain_decomposition,only: myid
    implicit none
    integer, intent(in) :: ns
    integer:: k, nlost
    nlost=0
    do k=1,nm_gk(ns)
       if(touch_bdry_gc(k,ns).eqv..true.) then
          nlost=nlost+1
       endif
    enddo
    write(*,*) 'number of lost gk markers=', nlost, 'myid=',myid, 'species=',ns, ',total number=', nm_gk(ns)
  end subroutine count_lost_markers_gk

end module gk_trajectory_pusher
