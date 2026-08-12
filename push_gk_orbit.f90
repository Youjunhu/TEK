module gk_trajectory_pusher
contains
  pure subroutine push_gc(ns, lost, dt, xdrift0, zdrift0, ydrift0, mirror_force,&
       & xdrift1, zdrift1, ydrift1, x_old, z_old, y_old, vpar_old, weight, x_new, z_new, y_new, vpar_new)
    use constants, only: p_, pi, one_half
    use gk_module, only: nm_gk, gk_nonlinear
    use math, only: shift_toroidal
    use magnetic_coordinates, only : toroidal_range, xlow, xupp
    use magnetic_field, only: qfunc
    implicit none
    integer, intent(in) :: ns
    logical, intent(in) :: lost(:)
    real(p_), intent(in) :: dt
    real(p_), intent(in) :: xdrift0(:), zdrift0(:), ydrift0(:), mirror_force(:)
    real(p_), intent(in) :: xdrift1(:), zdrift1(:), ydrift1(:)
    real(p_), intent(in) :: x_old(:), y_old(:), z_old(:), vpar_old(:)
    real(p_), intent(inout) :: weight(:)
    real(p_), intent(out) :: x_new(:), y_new(:), z_new(:), vpar_new(:)
    real(p_) :: xdrift, zdrift, ydrift
    integer :: nm, k

    nm = nm_gk(ns)
    do k = 1, nm
       if(lost(k) .eqv. .true.) then
          z_new(k) = z_old(k)
          x_new(k) = x_old(k)
          cycle
       endif

       if(gk_nonlinear(ns)==1) then
          xdrift = xdrift0(k) + xdrift1(k)
          zdrift = zdrift0(k) + zdrift1(k)
          ydrift = ydrift0(k) + ydrift1(k)
       else
          xdrift = xdrift0(k)
          zdrift = zdrift0(k)
          ydrift = ydrift0(k)
       endif

       x_new(k) = x_old(k) + xdrift*dt
       z_new(k) = z_old(k) + zdrift*dt
       y_new(k) = y_old(k) + ydrift*dt
       vpar_new(k)  = vpar_old(k) + mirror_force(k)*dt

       ! Comment out the following if-structure if you do not want particle refilling:
!!$       if ((x_new(k) > xupp) .or. (x_new(k) < xlow) ) then
!!$          x_new(k) = x_old(k)
!!$          z_new(k) = - z_old(k)
!!$          y_new(k) = y_old(k) + 2*qfunc(x_old(k))*z_old(k)
!!$          weight(k) = 0
!!$        endif

       if ((z_new(k) .ge. pi) .or. (z_new(k) < -pi)) &
            & call shift_gc_theta_then_alpha(x_new(k), z_new(k), y_new(k))
       call shift_toroidal(y_new(k), toroidal_range)
    enddo

  end subroutine push_gc

  pure  subroutine shift_gc_theta_then_alpha(radcor, theta, alpha)
    !shift theta, and then shift alpha so that phi is not changed, i.e., keep the particle in the same sptial location
    use constants, only: p_, twopi, pi
    use magnetic_field, only: qfunc
    implicit none
    real(p_), intent(in) :: radcor
    real(p_), intent(inout) :: theta, alpha
    integer :: ishift

    ishift = floor(theta/twopi)
    theta = theta - ishift * twopi
    if(theta .ge. pi) then
       theta = theta - twopi
       ishift = ishift + 1
    endif

    !alpha shift is needed to make the cylindrical toroidal angle phi be the same as the value before the theta-shifting,
    !i.e., make the particle lie at the same sptaial location when doing the theta-shifting
    alpha = alpha + ishift*twopi*qfunc(radcor) 
  end subroutine shift_gc_theta_then_alpha

end module gk_trajectory_pusher
