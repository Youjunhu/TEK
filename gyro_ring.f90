module gyro_ring_mod
  use constants, only : p_, two, twopi
  real(p_), allocatable :: cos_gyro(:,:), sin_gyro(:,:)

contains
  subroutine set_gyro_phase()
    use domain_decomposition, only: myid
    use gk_module, only: gyro_npt, nmmax
    use math, only: sub_random_yj
    implicit none
    integer :: i, j, iseed, next_seed
    real(p_) :: gyro_angle0, gyro_angle, tmp

    allocate(cos_gyro(gyro_npt,nmmax))
    allocate(sin_gyro(gyro_npt,nmmax))

    iseed=-(1177+myid*3) !set the iseed in different procs
    call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 
    do i=1, nmmax
       call sub_random_yj(0,next_seed,gyro_angle0) !0 means using last random number as iseed, first gyro-angle
       gyro_angle0 = gyro_angle0*twopi !scaled to [0:twopi]
       !gyro_angle0=0 !use this if we do not want the first gyro-angle to be random
       do j = 1, gyro_npt ! gyro-angles are uniform distributed relatively to the first one.
          gyro_angle =  gyro_angle0 + twopi/gyro_npt*(j-1)
          cos_gyro(j,i) = cos(gyro_angle)
          sin_gyro(j,i) = sin(gyro_angle)
       enddo
    enddo
  end subroutine set_gyro_phase

  subroutine gyro_ring(ns, touch_bdry, radcor, alpha, theta, mu, x_ring, y_ring, z_ring)
    !output are used by gyro_average and deposit_gk
    use gk_module, only:  gk_flr, nm_gk
    implicit none
    integer, intent(in) :: ns
    logical, intent(in) :: touch_bdry(:)
    real(p_), intent(in) :: theta(:), radcor(:), alpha(:), mu(:)
    real(p_), intent(out) :: x_ring(:, :), y_ring(:, :), z_ring(:, :)
    integer :: nm, k

    nm= nm_gk(ns)

    if(gk_flr(ns) .eqv. .false.) then
       do k = 1, nm
          if (touch_bdry(k) .eqv. .true.) cycle
          x_ring(1,k) = radcor(k)
          y_ring(1,k) = alpha(k)
          z_ring(1,k) = theta(k)
       enddo
    else
       do k =1, nm
          if (touch_bdry(k) .eqv. .true.) cycle
          call gyro_ring_core(ns, k, radcor(k), alpha(k), theta(k), mu(k), x_ring(:,k), y_ring(:,k), z_ring(:,k))
       enddo
    endif
  end subroutine gyro_ring


  subroutine gyro_ring_core(ns, k, x0, y0, z0, mu, x, y, z)
    use magnetic_coordinates, only: mpol, nrad, zgrid, xgrid, dtheta, dradcor, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, &
         & pfn_inner, pfn_bdry, toroidal_range, GSpsi_prime
    use table_in_mc, only : bdgxcgz, b_mc
    use math,only: shift_toroidal
    use gk_module, only:   mass_gk, charge_gk, vn_gk, gyro_npt
    use domain_decomposition, only: theta_start, dtheta2
    use math, only: shift_toroidal
    use interpolate_module, only : linear_2d_interpolate0, locate
    implicit none
    integer, intent(in) :: ns, k
    real(p_), intent(in) :: x0, y0, z0, mu
    real(p_), intent(out) :: x(:), y(:), z(:)
    real(p_) :: bval0, gx0, gxdgy0, gxdgz0, bdgxcgz0
    real(p_) :: vper, gyro_radius, dy1, dy2, dz1, dz2
    integer :: i, j, kp

    call locate(mpol, zgrid, dtheta, z0, i)
    call locate(nrad, xgrid, dradcor,x0, j)
    call linear_2d_interpolate0(mpol,nrad,zgrid,xgrid,dtheta, dradcor, grad_psi,&
         & z0, x0,i,j, gx0)
    call linear_2d_interpolate0(mpol,nrad,zgrid,xgrid,dtheta, dradcor, grad_psi_dot_grad_alpha,&
         & z0, x0,i,j,gxdgy0)
    call linear_2d_interpolate0(mpol,nrad,zgrid,xgrid,dtheta, dradcor, grad_psi_dot_grad_theta,&
         & z0, x0,i,j, gxdgz0)
    call linear_2d_interpolate0(mpol, nrad, zgrid, xgrid, dtheta, dradcor, bdgxcgz,&
         & z0, x0,i,j, bdgxcgz0)
    call linear_2d_interpolate0(mpol, nrad, zgrid, xgrid, dtheta, dradcor, b_mc,&
         & z0, x0,i,j,bval0)

    vper=sqrt(two*bval0*mu)
    gyro_radius=vper*vn_gk(ns)/(bval0*abs(charge_gk(ns))/mass_gk(ns))

    dy1 = gyro_radius*gxdgy0/gx0
    dy2 = gyro_radius*bval0/(GSpsi_prime*gx0)
    dz1 = gyro_radius*gxdgz0/gx0
    dz2 = gyro_radius*bdgxcgz0/(bval0*gx0)

    do kp = 1, gyro_npt
       x(kp)= x0 + gyro_radius*cos_gyro(kp,k)*gx0
       x(kp) = max(x(kp), pfn_inner) !to prevent exceeding the radial computational box
       x(kp) = min(x(kp), pfn_bdry) 
       y(kp) = y0 + cos_gyro(kp,k)*dy1 + sin_gyro(kp,k)*dy2
       call shift_toroidal(y(kp), toroidal_range)
       z(kp) = z0 + cos_gyro(kp,k)*dz1 + sin_gyro(kp,k)*dz2
       z(kp) = max(z(kp), theta_start)
       z(kp) = min(z(kp), theta_start+dtheta2)
    enddo
  end subroutine gyro_ring_core
end module gyro_ring_mod
