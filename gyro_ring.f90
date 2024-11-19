module gyro_ring_mod
  use constants, only : p_, twopi
  real(p_),allocatable :: cos_gyro(:,:), sin_gyro(:,:)

contains
  subroutine set_gyro_phase()
    use domain_decomposition, only : myid
    use gk_module, only : gyro_npt, nmmax
    use math, only : sub_random_yj
    implicit none
    integer :: i, j, iseed, next_seed
    real(p_) :: gyro_angle0, gyro_angle, tmp

    allocate(cos_gyro(gyro_npt,nmmax))
    allocate(sin_gyro(gyro_npt,nmmax))

    iseed=-(1177+myid*3) !set the iseed in different procs
    call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 
    do i=1, nmmax
       call sub_random_yj(0,next_seed,gyro_angle0) !0 means using last random number as iseed, first gyro-angle
       gyro_angle0=gyro_angle0*twopi !scaled to [0:twopi]
       !gyro_angle0=0 !use this if we do not want random gyro-phase
       do j=1, gyro_npt !first gyro-angle on the gyro-ring may be random, but the remainder are uniform distributed relatively to the first one.      
          gyro_angle =  gyro_angle0 + twopi/gyro_npt*(j-1)  !correct
          cos_gyro(j,i) = cos(gyro_angle)
          sin_gyro(j,i) = sin(gyro_angle)
       enddo
    enddo
  end subroutine set_gyro_phase

  subroutine gyro_ring(ns,touch_bdry_e, theta_e, radcor_e, alpha_e, mu_e, x_ring, y_ring)
    use constants, only : p_
    use interpolate_module, only : linear_2d_interpolation0, locate
    use table_in_mc,only: b_mc
    use magnetic_coordinates,only: mpol,nflux,theta_1d_array,radcor_1d_array, toroidal_range, dtheta, dradcor, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha, GSpsi_prime
    use math,only: shift_to_specified_toroidal_range
    use gk_module,only:  gk_flr, nm_gk
    implicit none
    integer, intent(in) :: ns
    logical, intent(in) :: touch_bdry_e(:)
    real(p_), intent(in) :: theta_e(:), radcor_e(:), alpha_e(:), mu_e(:)
    real(p_), intent(out) :: x_ring(:, :), y_ring(:, :)
    real(p_) :: bval0, grad_psi0, grad_psi_dot_grad_alpha0
    integer :: nm, k, i,j

    nm= nm_gk(ns)

    if(gk_flr(ns) .eqv. .false.) then
       do k=1,nm
       if (touch_bdry_e(k) .eqv. .true.) cycle
          x_ring(1,k) = radcor_e(k)
          y_ring(1,k) = alpha_e(k)
       enddo
       return
    end if

    do k=1,nm
       if (touch_bdry_e(k) .eqv. .true.) cycle
       call locate(mpol,theta_1d_array, dtheta,theta_e(k),i)
       call locate(nflux,radcor_1d_array,dradcor,radcor_e(k),j)
       call linear_2d_interpolation0(mpol,nflux,theta_1d_array,radcor_1d_array,dtheta, dradcor, grad_psi,&
            & theta_e(k),radcor_e(k),i,j,grad_psi0)
       call linear_2d_interpolation0(mpol,nflux,theta_1d_array,radcor_1d_array,dtheta, dradcor, grad_psi_dot_grad_alpha,&
            & theta_e(k),radcor_e(k),i,j,grad_psi_dot_grad_alpha0)
       call linear_2d_interpolation0(mpol,nflux,theta_1d_array,radcor_1d_array,dtheta, dradcor, b_mc,&
            & theta_e(k),radcor_e(k),i,j,bval0)

       call gyro_ring_core(ns,k, radcor_e(k), alpha_e(k), mu_e(k),bval0, grad_psi0, &
            grad_psi_dot_grad_alpha0, GSpsi_prime, x_ring(:,k),y_ring(:,k))       
    enddo
  end subroutine gyro_ring


  subroutine gyro_ring_core(ns,k, x0,y0, mu, bval, grad_psi, grad_psi_dot_grad_alpha, GSpsi_prime,x,y)
    use constants,only: p_, two
    use magnetic_coordinates,only: pfn_inner, pfn_bdry
    use gk_module,only:   mass_e, charge_e, vn_e, gyro_npt
    use math,only: shift_to_specified_toroidal_range
    implicit none
    integer,intent(in) :: ns, k
    real(p_),intent(in) :: x0,y0,mu
    real(p_),intent(in) :: bval, grad_psi,grad_psi_dot_grad_alpha, GSpsi_prime
    real(p_),intent(out) :: x(:), y(:)
    real(p_) :: vper, gyro_radius, dalpha1,dalpha2
    integer :: j

    vper=sqrt(two*bval*mu)
    gyro_radius=vper*vn_e(ns)/(bval*abs(charge_e(ns))/mass_e(ns))
    dalpha1=gyro_radius*grad_psi_dot_grad_alpha/grad_psi
    dalpha2=gyro_radius*bval/(GSpsi_prime*grad_psi)

    do j = 1, gyro_npt
       x(j)= x0 + gyro_radius*cos_gyro(j,k)*grad_psi
       x(j)=min(x(j),pfn_bdry) !to prevent exceeding the radial computational box
       x(j)=max(x(j),pfn_inner)
       y(j)= y0 + cos_gyro(j,k)*dalpha1+sin_gyro(j,k)*dalpha2
       call shift_to_specified_toroidal_range(y(j))
    enddo
  end subroutine gyro_ring_core
end module gyro_ring_mod
