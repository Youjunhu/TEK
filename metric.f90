module calculate_toroidal_shift_module
contains
subroutine calculate_toroidal_shift(psival,x_contour,z_contour,np_lcfs,end_i,r,z,tor_shift) !tor_shift=q*delta. zeta=phi-tor_shift
  use constants,only: p_, one_half
  ! use poloidal_flux_2d,only: nx,nz,xarray,zarray
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_GSpsi_prime, i_theta_zero
  use radial_module,only:psi_axis,psi_lcfs
  use domain_decomposition,only:myid
  use magnetic_field, only : psi_z_func,psi_r_func, g_func !toroidal field function
  implicit none
  real(p_),intent(in):: psival
  integer,intent(in):: np_lcfs,end_i
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  real(p_),intent(in):: r,z
  real(p_),intent(out):: tor_shift

  real(p_):: x_mid,z_mid,gx0,dl, g_value
  real(p_):: pfn,mr,costh,sinth,r0
  integer:: i

  g_value=g_func(psival)

  tor_shift=0._p_
  if (end_i.lt.i_theta_zero) then 
     do i=i_theta_zero-1,end_i+1,-1
        x_mid=(x_contour(i)+x_contour(i-1))*one_half
        z_mid=(z_contour(i)+z_contour(i-1))*one_half
        gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
        dl=-sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
        tor_shift=tor_shift+g_value/(x_mid*gx0)*dl
     enddo
     x_mid=(x_contour(end_i+1)+r)*one_half
     z_mid=(z_contour(end_i+1)+z)*one_half
     gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=-sqrt((x_contour(end_i+1)-r)**2+(z_contour(end_i+1)-z)**2)
     tor_shift=tor_shift+g_value/(x_mid*gx0)*dl

  elseif(end_i.ge.i_theta_zero) then
     x_contour(end_i+1)=r !replace No. end_i+1 point by the given point (r,z)
     z_contour(end_i+1)=z
     do i=i_theta_zero, end_i
        x_mid=(x_contour(i)+x_contour(i+1))*one_half
        z_mid=(z_contour(i)+z_contour(i+1))*one_half
        gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
        dl=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
        tor_shift=tor_shift+g_value/(x_mid*gx0)*dl
     enddo
  endif
  !--for testing----analytical formula, for concentric circular configuration
!!$  r0=1.32_p_
!!$  mr=sqrt((x_contour(end_i+1)-r0)**2+(z_contour(end_i+1)-0._p_)**2)
!!$  costh=(x_contour(end_i+1)-r0)/mr
!!$  sinth=z_contour(end_i+1)/mr
!!$  pfn=(psival-psi_axis)/(psi_lcfs-psi_axis)
!!$  tor_shift=2*qfunc(pfn)*atan((r0-mr)/sqrt(r0**2-mr**2)*sinth/(costh+1))
  !----for tesing----------------------

!!$  do i=2,end_i+1
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     tor_shift=tor_shift+g_value/(x_mid*gx0(i-1))*dl(i-1)
!!$  enddo
  tor_shift=tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
end subroutine calculate_toroidal_shift


pure subroutine calculate_toroidal_shift2(j, g_value, x_contour,z_contour,jacob, m, tor_shift, tor_shift_left_bdry_minus_one)
  !cf. calculate_toroidal_shift, the differences are (1) g_value instead of psival is provided as input;
  !(2) compute a series of tor_shift on a magnetic surface (rather than a single value). The range of theta is [-pi:pi]
  use constants,only: p_, one_half, two
  use magnetic_coordinates,only: dtheta, sign_of_jacobian,sign_of_GSpsi_prime,i_theta_zero, GSpsi_array, GSpsi_prime, qhat
  use domain_decomposition,only: dtheta2, multi_eq_cells
  use magnetic_field, only: psi_z_func, psi_r_func
  implicit none
  integer,intent(in) :: j, m
  real(p_),intent(in) :: g_value
  real(p_),intent(in) :: x_contour(m), z_contour(m), jacob(m)
  real(p_),intent(out) :: tor_shift(m), tor_shift_left_bdry_minus_one
  real(p_) :: gx0, dl, sum, x_mid, z_mid, jac_mid
  integer :: i

!!$  tor_shift(i_theta_zero) = 0._p_
!!$  do i=i_theta_zero+1, m !for theta in (0:pi]
!!$     x_mid=(x_contour(i-1)+x_contour(i))/two
!!$     z_mid=(z_contour(i-1)+z_contour(i))/two
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
!!$     tor_shift(i) = tor_shift(i-1) + g_value/(x_mid*gx0)*dl
!!$  enddo
!!$
!!$  do i=i_theta_zero-1, 1, -1 !for theta in (0:-pi]
!!$     x_mid=(x_contour(i+1)+x_contour(i))/two
!!$     z_mid=(z_contour(i+1)+z_contour(i))/two
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2) 
!!$     dl=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
!!$     tor_shift(i) = tor_shift(i+1) - g_value/(x_mid*gx0)*dl
!!$  enddo
!!$  tor_shift = tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign

  !-------------another method------------:
  tor_shift(i_theta_zero) = 0._p_
  do i=i_theta_zero+1, m !for theta in (0:pi]
     x_mid=(x_contour(i-1)+x_contour(i))/two
     jac_mid=(jacob(i-1)+jacob(i))/two
     !tor_shift(i) = tor_shift(i-1) - g_value/(x_mid**2)*jac_mid/GSpsi_prime*dtheta
     tor_shift(i) = tor_shift(i-1) + qhat(i,j)*dtheta
  enddo

  do i=i_theta_zero-1, 1, -1 !for theta in (0:-pi]
     x_mid=(x_contour(i+1)+x_contour(i))/two
     jac_mid=(jacob(i+1)+jacob(i))/two
     tor_shift(i) = tor_shift(i+1) + g_value/(x_mid**2)*jac_mid/GSpsi_prime*dtheta
  enddo
!---------------------------------------
  !another way of calculating tor_shift
!!$ tor_shift(1)=0._p_
!!$  do i=1,m-1 
!!$     x_mid=(x_contour(i+1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i+1)+z_contour(i))*one_half
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i+1)-x_contour(i))**2+(z_contour(i+1)-z_contour(i))**2)
!!$     tor_shift(i+1)=tor_shift(i)+g_value/(x_mid*gx0)*dl
!!$  enddo
!!$  tor_shift=tor_shift-tor_shift(i_theta_zero)

!!$  tor_shift(1)=0._p_
!!$  do i=2,m
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i-1)+z_contour(i))*one_half
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
!!$     tor_shift(i)=tor_shift(i-1)+g_value/(x_mid*gx0)*dl
!!$  enddo
  !tor_shift = tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign

  tor_shift_left_bdry_minus_one = tor_shift(1) - (tor_shift(m)-tor_shift(m-multi_eq_cells))

end subroutine calculate_toroidal_shift2


subroutine calculate_toroidal_shift_at_theta_cut(psival,x_contour,z_contour,np,tor_shift_a,tor_shift_b)
  use constants,only: p_, one_half
  use magnetic_coordinates,only: sign_of_jacobian, sign_of_GSpsi_prime, i_theta_zero
  use domain_decomposition,only:myid
  use magnetic_field, only : psi_z_func,psi_r_func, g_func
  implicit none
  real(p_),intent(in):: psival
  integer,intent(in):: np
  real(p_),intent(in):: x_contour(np),z_contour(np)
  real(p_),intent(out):: tor_shift_a, tor_shift_b
  real(p_):: x_mid,z_mid,gx0,dl, g_value
  integer:: i

  g_value=g_func(psival)

  tor_shift_b=0._p_
  do i=i_theta_zero+1,np ! for location near the cut above the midplane
     x_mid=(x_contour(i)+x_contour(i-1))*one_half
     z_mid=(z_contour(i)+z_contour(i-1))*one_half
     gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
     tor_shift_b=tor_shift_b+g_value/(x_mid*gx0)*dl
  enddo

  tor_shift_a=0._p_
  do i=i_theta_zero,2,-1 ! for location near the cut below the midplane
     x_mid=(x_contour(i)+x_contour(i-1))*one_half
     z_mid=(z_contour(i)+z_contour(i-1))*one_half
     gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=-sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
     tor_shift_a=tor_shift_a+g_value/(x_mid*gx0)*dl
  enddo

  tor_shift_a=tor_shift_a*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
  tor_shift_b=tor_shift_b*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
end subroutine calculate_toroidal_shift_at_theta_cut
end module calculate_toroidal_shift_module


subroutine calculate_metric()
  use constants, only: p_, twopi
  use magnetic_coordinates,only: mpol,nrad,r_mc,z_mc,dtheta,dradcor, &
       & tor_shift_mc, tor_shift_mc_left_bdry_minus_one, & !output
       & qhat, jacobian, & !as output
       & pfn_inner,pfn_bdry,GSpsi_array, &
       & zgrid,xgrid, &
       & sign_of_jacobian,sign_of_GSpsi_prime, GSpsi_prime, &
       & grad_psi, grad_alpha,grad_theta,&
       & grad_psi_dot_grad_alpha,grad_psi_dot_grad_theta,grad_alpha_dot_grad_theta !output
  use radial_module, only: psi_lcfs,psi_axis
  use magnetic_field, only : g_func, qfunc
  use calculate_toroidal_shift_module, only : calculate_toroidal_shift2
  use domain_decomposition, only: myid
  use control_parameters, only: diagnosis
  implicit none
  integer :: i, j, ierr
  real(p_) :: g

  call calculate_gradients_of_psi_and_theta() !computing gradients at mc gridpoints

  if(myid.eq.0) write(*,*) 'max_jacobian=',maxval(jacobian),'min_jacobian=',minval(jacobian)
  !abs_jacobian_max=maxval(abs(jacobian))
  !abs_jacobian_min=minval(abs(jacobian))
  !write(*,*) 'max_abs(jacobian)=',maxval(abs(jacobian)),'min_abs(jacobian)=',minval(abs(jacobian))

  sign_of_jacobian=sign(1._p_,jacobian(mpol/3,nrad/3))
  sign_of_GSpsi_prime=sign(1._p_,GSpsi_prime) !radial coordinator is assumed to be pfn
  call cal_q_with_correct_sign() ! after this call, function "qfunc" is ready to be used
  if(myid.eq.0) write(*,*) 'q_inner=',qfunc(pfn_inner),'q_bdry=',qfunc(pfn_bdry)

  allocate(tor_shift_mc(mpol,nrad))
  allocate(tor_shift_mc_left_bdry_minus_one(nrad)) !at theta=-pi-dtheta

  ! local safety factor
  allocate(qhat(mpol,nrad))
  do j=1,nrad
     g=g_func(GSpsi_array(j))
     do i =1,mpol
        qhat(i,j) = -g/r_mc(i,j)**2*jacobian(i,j)/GSpsi_prime
     enddo
  enddo
  
  do j=1,nrad
     g=g_func(GSpsi_array(j))
     call calculate_toroidal_shift2(j, g, r_mc(:,j), z_mc(:,j), jacobian(:,j), mpol, &
          & tor_shift_mc(:,j), tor_shift_mc_left_bdry_minus_one(j))

  enddo

!!$ if(myid.eq.0) call draw_grids_on_theta_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc)
!!$ if(myid.eq.0) call draw_alpha_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc)
!!$ if(myid.eq.0) call draw_alpha_contours_on_a_magnetic_surface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !i.e., magnetic field lines

  call calculate_gradients_of_generalized_toroidal_angle()

  if((myid==0).and. (diagnosis .eqv. .true.)) then
     call diagnostic1()
     call diagnostic2()
  endif

contains
  subroutine diagnostic1()
    integer:: i,u
    write(*,*) 'maximum of tor_shift_mc=',maxval(tor_shift_mc),'minimum of tor_shift_mc=',minval(tor_shift_mc)
    open(newunit=u,file='tor_shift_mc.txt')
    do j=1,nrad
       do i=1,mpol
          write(u,*) r_mc(i,j),z_mc(i,j),tor_shift_mc(i,j), qhat(i,j)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)
  end subroutine diagnostic1

  subroutine diagnostic2()
    use magnetic_coordinates, only : pfn, dl_mc
    use func_in_mc, only: minor_r_radcor
    integer:: i,u

    real(p_) :: s, q2
    open(newunit=u,file="mc_derivatives_mc3.txt")
    do i=1,mpol
       do j=1,nrad
          write(u,'(2i8.4,9(1pe14.5))')  i,j,grad_alpha(i,j),grad_psi_dot_grad_alpha(i,j),&
               & grad_psi(i,j),jacobian(i,j),grad_alpha_dot_grad_theta(i,j),&
               & tor_shift_mc(i,j),grad_psi_dot_grad_theta(i,j),zgrid(i),xgrid(j)
       enddo
       write(u,*)
    enddo
    close(u)
!!$     do j=1,nrad
!!$        write(*,*) j,xgrid(j),minor_r_radcor(xgrid(j)),r_mag_surf0(1,j)-r_axis,&
!!$             & minor_r_prime(xgrid(j))
!!$     enddo
    open(newunit=u,file="q2.txt")
    do j=1,nrad
       s=0
       do i=1,mpol-1
          s = s + g_func(GSpsi_array(j))/(grad_psi(i,j)*(psi_lcfs-psi_axis)*r_mc(i,j))*dl_mc(i,j)
       enddo
       q2 = s/twopi
       write(u,'(20(ES16.4E2))') pfn(j), minor_r_radcor(pfn(j)),  q2, qfunc(pfn(j))
    enddo
    close(u)
  end subroutine diagnostic2

end subroutine calculate_metric

subroutine cal_q_with_correct_sign()
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_GSpsi_prime
  use radial_module, only : npsi,sign_bphi, qpsi
  use radial_module, only : q_with_sign !as output
  implicit none
  allocate(q_with_sign(npsi))
  q_with_sign=abs(qpsi)*(-sign_bphi)*sign_of_jacobian*sign_of_GSpsi_prime !include the correct sign
  !write(*,*) 'q_with_sign=', q_with_sign(npsi/2)
end subroutine cal_q_with_correct_sign


subroutine calculate_gradients_of_psi_and_theta()
  use constants, only: p_, zero, one, two
  use magnetic_coordinates, only: m=>mpol, n=>nrad, r=>r_mc, z=>z_mc, dtheta, dradcor,&
   & jacobian, abs_jacobian, av_jacobian, & !as output
   & rth, zth, rpsi, zpsi, grad_psi, grad_theta, grad_psi_dot_grad_theta, & !as output
   & grad_psi_r, grad_psi_z, grad_theta_r, grad_theta_z !as output
  use domain_decomposition, only: myid
  !use magnetic_field, only: minor_r_prime, minor_r_radcor
  implicit none
  integer :: i,j,u
  real(p_) ::  minor_r_prime_val
  
  allocate(jacobian(m,n))
  allocate(abs_jacobian(m,n))
  allocate(av_jacobian(n))
  allocate(grad_psi(m,n))
  allocate(grad_theta(m,n))
  allocate(grad_psi_dot_grad_theta(m,n))
  allocate(grad_theta_r(m,n))
  allocate(grad_theta_z(m,n))
  allocate(grad_psi_r(m,n))
  allocate(grad_psi_z(m,n))
  allocate(rpsi(m,n))
  allocate(zpsi(m,n))
  allocate(rth(m,n))
  allocate(zth(m,n))

  !call partial_derivative_in_mc(m,n,r,z,dtheta,dradcor,rpsi,rth,zpsi,zth,jacobian)
  call partial_derivative_in_mc2(m,n,r,z,dtheta,dradcor,rpsi,rth,zpsi,zth,jacobian) 

  abs_jacobian = abs(jacobian)
  
  grad_psi=r/abs_jacobian*sqrt(zth**2+rth**2)
  grad_theta=r/abs_jacobian*sqrt(zpsi**2+rpsi**2)
  grad_psi_dot_grad_theta=-(r/jacobian)**2*(zth*zpsi+rth*rpsi)

  grad_psi_r = -r/jacobian*zth
  grad_psi_z = r/jacobian*rth
  grad_theta_r = r/jacobian*zpsi
  grad_theta_z = -r/jacobian*rpsi

  do j = 1,n
     av_jacobian(j)=sum(abs_jacobian(1:m-1,j))/(m-1) !uniform theta grid is assumed
  enddo
!  if(myid.eq.0) call plot_psi_r_z_theta_r_z_mc()

contains
!!$  subroutine plot_psi_r_z_theta_r_z_mc()
!!$    use magnetic_coordinates,only: zgrid,xgrid
!!$    integer:: u
!!$    real(p_):: minor_r_val,minor_r_prime_val
!!$
!!$    open(newunit=u, file='mc_derivatives_mc1.txt')
!!$    do i=1,m
!!$       do j=1,n
!!$          write(u,*) r(i,j),z(i,j), jacobian(i,j)
!!$
!!$          !minor_r_prime_val=minor_r_prime(xgrid(j))
!!$          !write(u,*) r(i,j),z(i,j),jacobian(i,j)/minor_r_prime_val !for equal arc-length poloidal angle, jacobian(i,j)/minor_r_prime_val is equal -R*minor_r.
!!$          !write(u,*) r(i,j),z(i,j),grad_psi_r(i,j)*minor_r_prime_val, grad_psi_z(i,j)*minor_r_prime_val,&
!!$          !write(u,*) r(i,j),z(i,j),grad_theta_r(i,j),grad_theta_z(i,j),grad_theta(i,j),jacobian(i,j)
!!$       enddo
!!$       write(u,*)
!!$    enddo
!!$    close(u)
!!$
!!$    open(newunit=u, file='psi_r_z_theta_r_z_mc_analytic.txt') !for concentric-circular magnetic field (equal-arc-length theta)
!!$    do i=1,m
!!$       do j=1,n
!!$          !minor_r_val=minor_r_radcor(xgrid(j))
!!$          write(u,*) r(i,j),z(i,j),cos(zgrid(i)),sin(zgrid(i)),&
!!$               &-sin(zgrid(i))/minor_r_val,cos(zgrid(i))/minor_r_val,-R(i,j)*minor_r_val
!!$       enddo
!!$       write(u,*)
!!$    enddo
!!$    close(u)
!!$  end subroutine plot_psi_r_z_theta_r_z_mc
end subroutine calculate_gradients_of_psi_and_theta


subroutine calculate_gradients_of_generalized_toroidal_angle() !compute the gradient of the generalized toroidal angle alpha
  use constants,only: p_, zero,one,two
  use control_parameters, only : diagnosis
  use magnetic_coordinates,only: m=>mpol, n=>nrad, r=>r_mc, z=>z_mc, &
  & tor_shift_mc,jacobian,dtheta,dradcor, &
  & rpsi,zpsi,rth,zth, &
  & grad_alpha_r, grad_alpha_z, grad_alpha_phi, & !as output 
  & grad_alpha, grad_psi_dot_grad_alpha, grad_alpha_dot_grad_theta, & !as output
  & grad_alpha_r_left_bdry_minus_one, grad_alpha_z_left_bdry_minus_one !as output 
  use domain_decomposition,only: myid, dtheta2
  implicit none
  real(p_) :: tor_shift_psi(m,n), tor_shift_th(m,n)
  real(p_) :: tor_shift_psi_left_bdry_minus_one(n)
  integer :: i,j

  allocate(grad_alpha_r(m,n))
  allocate(grad_alpha_z(m,n))
  allocate(grad_alpha_phi(m,n))
  allocate(grad_alpha(m,n))
  allocate(grad_psi_dot_grad_alpha(m,n))
  allocate(grad_alpha_dot_grad_theta(m,n))
  allocate(grad_alpha_r_left_bdry_minus_one(n))
  allocate(grad_alpha_z_left_bdry_minus_one(n))

  call partial_derivative_of_tor_shift_in_mc2(m,n,tor_shift_mc,dtheta,dradcor,tor_shift_psi,tor_shift_th, &
       & tor_shift_psi_left_bdry_minus_one) 
  grad_alpha_r=tor_shift_psi*r/jacobian*zth-tor_shift_th*r/jacobian*zpsi
  grad_alpha_z=tor_shift_th*r/jacobian*rpsi-tor_shift_psi*r/jacobian*rth
  grad_alpha_phi=one/r
  grad_alpha=sqrt(grad_alpha_phi**2+grad_alpha_r**2+grad_alpha_z**2)
  grad_psi_dot_grad_alpha=-r/jacobian*zth*grad_alpha_r+r/jacobian*rth*grad_alpha_z
  grad_alpha_dot_grad_theta=grad_alpha_r*r/jacobian*zpsi-grad_alpha_z*r/jacobian*rpsi

  !i=m-1 !wrong
  i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
  do j=1,n
     grad_alpha_r_left_bdry_minus_one(j)=tor_shift_psi_left_bdry_minus_one(j)*r(i,j)/jacobian(i,j)*zth(i,j)&
          & -tor_shift_th(i,j)*r(i,j)/jacobian(i,j)*zpsi(i,j)
     grad_alpha_z_left_bdry_minus_one(j)=tor_shift_th(i,j)*r(i,j)/jacobian(i,j)*rpsi(i,j)&
          & -tor_shift_psi_left_bdry_minus_one(j)*r(i,j)/jacobian(i,j)*rth(i,j)
  enddo
  !if(myid==0 .and. (diagnosis .eqv. .true.)) call plot_alpha_r_z_mc()
  if(myid==0 .and. (diagnosis .eqv. .true.)) call verification2()
  
contains
  ! subroutine plot_alpha_r_z_mc()
  !   use magnetic_coordinates,only: xgrid,zgrid
  !   use magnetic_field,only: minor_r_radcor, qfunc

  !   integer :: u
  !   real(p_) :: theta,minor_r,major_r0,major_r,q,local_q,dq_dr,factor
  !   real(p_) :: ddelta_minor_r,alpha_r,alpha_z
  !   real(p_) :: tmp,tmp2,dpsi_dr
  !   real(p_), parameter :: g0=1.32*1.91,psi_axis=0.,  psi_lcfs=  0.146018378_p_
  !   open(newunit=u,file='mc_derivatives_mc2.txt')
  !   do i=1,m
  !      do j=1,n
  !         !write(u,*) r(i,j),z(i,j),grad_alpha_r(i,j),grad_alpha_z(i,j), grad_alpha(i,j)
  !         write(u,*) r(i,j),z(i,j),tor_shift_th(i,j),tor_shift_psi(i,j)
  !         !write(u,*) r(i,j),z(i,j),tor_shift_mc(i,j)
  !      enddo
  !      write(u,*)
  !      write(u,*)
  !   enddo

  !   i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
  !   do j=1,n
  !      write(u,*) r(i,j),z(i,j), tor_shift_th(i,j),tor_shift_psi_left_bdry_minus_one(j)
  !   enddo
  !   close(u)
  !   !write(*,*) 'maxval(grad_alpha)=',maxval(grad_alpha),'minval(grad_alpha)=',minval(grad_alpha)
  !   open(newunit=u,file='alpha_r_z_mc_analytic.txt') !for concentric-circular magnetic field (equal-arc-length theta and geometric minor radius as the flux surface label)
  !   do i=1,m
  !      theta=zgrid(i)
  !      do j=1,n
  !         major_r0=1.32_p_!for DIII-D cyclone base case parameter
  !         minor_r=minor_r_radcor(xgrid(j))
  !         major_r=major_r0+minor_r*cos(theta)
  !         q=qfunc(xgrid(j))
  !         local_q=q*sqrt(major_r0**2-minor_r**2)/major_r
  !         dq_dr=0.78*1.4/0.24_p_ !for DIII-D cyclone base case parameter, a linear q profile is assumed, the same as what Ben did
  !         factor=(major_r0-minor_r)/sqrt(major_r0**2-minor_r**2)*tan(theta/two)
  !         ddelta_minor_r=two*dq_dr*atan(factor)+two*q/(1+factor**2)*tan(theta/2._p_)*&
  !              & (-major_r0)/((major_r0+minor_r)*sqrt(major_r0**2-minor_r**2))
  !         alpha_r=-ddelta_minor_r*cos(theta)+local_q*sin(theta)/minor_r
  !         alpha_z=-local_q*cos(theta)/minor_r-ddelta_minor_r*sin(theta)
  !         !write(u,*) r(i,j),z(i,j), alpha_r,alpha_z
  !         tmp=(major_r0-minor_r)/sqrt(major_r0**2-minor_r**2)
  !         tmp2=tan(theta/2)
  !         !dpsi_dr=one/minor_r_prime(xgrid(j)) 
  !         dpsi_dr=g0*minor_r/(q*sqrt(major_r0**2-minor_r**2)*(psi_lcfs-psi_axis)) !another way to calculate dpsi_dr
  !         !write(u,*) r(i,j),z(i,j), local_q, 2*q*tmp*(tan(theta/2)**2/2+0.5)/(tmp**2*tan(theta/2)**2+1)
  !         write(u,*) r(i,j),z(i,j), local_q, ddelta_minor_r/dpsi_dr !local_q is equal to ddelta_dtheta
  !         !write(u,*) r(i,j),z(i,j), (2*atan(tmp2*(major_r0 - minor_r)/sqrt(major_r0**2 - minor_r**2))*dq_dr&
  !         !     &  + 2*(tmp2*minor_r*(major_r0 - minor_r)/(major_r0**2 - minor_r**2)**(3./2) &
  !         !     & - tmp2/sqrt(major_r0**2 - minor_r**2))*q/(tmp2**2*(major_r0 - minor_r)**2&
  !         !     & /(major_r0**2 - minor_r**2) + 1))/dpsi_dr !ddelta/dradcor=ddelta/dr*(dr/dradcor), where the formula for computing ddelta/dr is obtained by using sympy

  !         !write(u,*) r(i,j),z(i,j), 2*q*atan(tmp*tan(theta/2))
  !      enddo
  !      write(u,*)
  !      write(u,*)
  !   enddo
  !   close(u)
  ! end subroutine plot_alpha_r_z_mc

  subroutine verification2() !verify grad_psi_cross_grad_alpha=B0/GSpsi_prime
    use magnetic_coordinates,only: xgrid, grad_psi_r,grad_psi_z, GSpsi_prime
    use math,only:  cross_product_in_cartesian
    use magnetic_field,only : br,bz,bphi
    real(p_):: ax,ay,az, cx,cy,cz,dx,dy,dz
    real(p_):: brval,bzval,bphival
    integer:: u,i,j
    open(newunit=u,file='cross_product_comparision.txt')

!    i=1
!    i=m

 i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
    do j=1,n
       ax=grad_psi_r(i,j)
       ay=0._p_
       az=grad_psi_z(i,j)

!!$       cx=grad_alpha_r(i,j)
!!$       cy=grad_alpha_phi(i,j)
!!$       cz=grad_alpha_z(i,j)

       cx=grad_alpha_r_left_bdry_minus_one(j)
       cy=grad_alpha_phi(i,j)
       cz=grad_alpha_z_left_bdry_minus_one(j)

       call cross_product_in_cartesian(ax,ay,az,cx,cy,cz,dx,dy,dz) !grad_psi_cross_grad_alpha

        brval=br(r(i,j),z(i,j))
        bzval=bz(r(i,j),z(i,j))
        bphival=bphi(r(i,j),z(i,j))
       write(u,*) dx,dy,dz, brval/GSpsi_prime,bphival/GSpsi_prime,bzval/GSpsi_prime
    enddo

close(u)
  end subroutine verification2
end subroutine calculate_gradients_of_generalized_toroidal_angle


subroutine partial_derivative_in_mc(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)
  !calculate the partial derivative of R and Z with respect to the magnetic cooordinates (psi,theta)
  !jacob is also calculated in this subroutine
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: rpsi(m,n),rth(m,n),zpsi(m,n),zth(m,n),jacob(m,n)
  real(p_):: tmp0
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        rpsi(i,j)=(r(i,j+1)-r(i,j-1))/(two*dpsi)
        zpsi(i,j)=(z(i,j+1)-z(i,j-1))/(two*dpsi)
     enddo

     !use linear interpolation to get the value  j=n
     tmp0=(r(i,n)-r(i,n-1))/dpsi
     rpsi(i,n)=two*tmp0-rpsi(i,n-1)
     tmp0=(z(i,n)-z(i,n-1))/dpsi
     zpsi(i,n)=two*tmp0-zpsi(i,n-1)

     !use linear interpolation to get the value j=1
     tmp0=(r(i,2)-r(i,1))/dpsi
     rpsi(i,1)=two*tmp0-rpsi(i,2)

     tmp0=(z(i,2)-z(i,1))/dpsi
     zpsi(i,1)=two*tmp0-zpsi(i,2)

  enddo

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        rth(i,j)= (r(i+1,j)-r(i-1,j))/(two*dtheta)
        zth(i,j)=(z(i+1,j)-z(i-1,j))/(two*dtheta)
     enddo

     !use peroidic property of r and z to calculate the partial derivative for boundary points at theta=0 and 2pi
     rth(1,j)=(r(2,j)-r(m-1,j))/(two*dtheta)
     zth(1,j)=(z(2,j)-z(m-1,j))/(two*dtheta)
     rth(m,j)=rth(1,j)
     zth(m,j)=zth(1,j)
  enddo

  !calculate the Jacobian:
  do i=1,m
!          do j=2,n !the jacobian at the magnetic axis is zero
     do j=1,n
        jacob(i,j)=r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,fai)
     enddo
     !jacob(i,n)=two*jacob(i,n-1)-jacob(i,n-2)
     !use linear interpolation to get the value of jacobian at the magnetic surface near the magnetic axis
  enddo

  !write(*,*) jacob(1:m,5)
end subroutine partial_derivative_in_mc

subroutine partial_derivative_in_mc2(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)
  !calculate the partial derivative of R and Z with respect to the magnetic cooordinates (psi,theta)
  !jacob is also calculated in this subroutine
  use constants,only: p_
  use magnetic_coordinates, only : radcor=>xgrid, theta=>zgrid
  use splines, only: spline3ders
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: r(m,n), z(m,n)
  real(p_), intent(in) :: dtheta, dpsi
  real(p_), intent(out) :: rpsi(m,n), rth(m,n), zpsi(m,n), zth(m,n), jacob(m,n)
  integer :: i, j

  do i = 1, m
     call spline3ders(radcor, r(i,:), radcor, dynew=rpsi(i,:))
     call spline3ders(radcor, z(i,:), radcor, dynew=zpsi(i,:))
  enddo

  do j=1,n
     call spline3ders(theta, r(:,j), theta, dynew=rth(:,j))
     call spline3ders(theta, z(:,j), theta, dynew=zth(:,j))
  enddo

  do i=1,m
     do j=1,n
        jacob(i,j) = r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,phi)
     enddo
     !jacob(i,1) = 2*jacob(i,2) - jacob(i,3)
  enddo
  
end subroutine partial_derivative_in_mc2


subroutine partial_derivative_of_tor_shift_in_mc(m,n,tor_shift,dtheta,dpsi,tor_shift_psi,tor_shift_th,&
     & tor_shift_psi_left_bdry_minus_one)
  !calculate the partial derivative with respect to the magnetic cooordinates (psi,theta)
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  use magnetic_coordinates,only: tor_shift_mc_left_bdry_minus_one
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: tor_shift(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: tor_shift_psi(m,n),tor_shift_th(m,n)
  real(p_),intent(out):: tor_shift_psi_left_bdry_minus_one(n)
  real(p_):: tmp0,twopi_q
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        tor_shift_psi(i,j)=(tor_shift(i,j+1)-tor_shift(i,j-1))/(two*dpsi)
     enddo
     !use linear interpolation to get the value at j=n
     tmp0=(tor_shift(i,n)-tor_shift(i,n-1))/dpsi
     tor_shift_psi(i,n)=two*tmp0-tor_shift_psi(i,n-1)
     !use linear interpolation to get the value j=1
     tmp0=(tor_shift(i,2)-tor_shift(i,1))/dpsi
     tor_shift_psi(i,1)=two*tmp0-tor_shift_psi(i,2)
  enddo

  do j=2,n-1
     tor_shift_psi_left_bdry_minus_one(j)=&
          & (tor_shift_mc_left_bdry_minus_one(j+1)-tor_shift_mc_left_bdry_minus_one(j-1))/(two*dpsi)
  enddo
  !use linear interpolation to get the value at j=n
  tmp0=(tor_shift_mc_left_bdry_minus_one(n)-tor_shift_mc_left_bdry_minus_one(n-1))/dpsi
  tor_shift_psi_left_bdry_minus_one(n)=two*tmp0-tor_shift_psi_left_bdry_minus_one(n-1)
  !use linear interpolation to get the value j=1
  tmp0=(tor_shift_mc_left_bdry_minus_one(2)-tor_shift_mc_left_bdry_minus_one(1))/dpsi
  tor_shift_psi_left_bdry_minus_one(1)=two*tmp0-tor_shift_psi_left_bdry_minus_one(2)

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        tor_shift_th(i,j)= (tor_shift(i+1,j)-tor_shift(i-1,j))/(two*dtheta)
     enddo
     !for boundary points at theta cut
!!$     twopi_q=tor_shift(m,j)
!!$     tor_shift_th(1,j)=(tor_shift(2,j)+twopi_q-tor_shift(m-1,j))/(two*dtheta)
!!$     tor_shift_th(m,j)=tor_shift_th(1,j)
     !tor_shift_th(1,j)=(tor_shift(2,j)-tor_shift_mc_left_bdry_minus_one(j))/(two*dtheta), wrong! since left_bdry_minus_one is defined on the grids for perturbations, which are different from equilibrium grids.
     tor_shift_th(1,j)=two*tor_shift_th(2,j)-tor_shift_th(3,j)   !use linear interpolation to get the value at i=1 and i=m
     tor_shift_th(m,j)=tor_shift_th(1,j) !tor_shift_th is periodic

  enddo
end subroutine partial_derivative_of_tor_shift_in_mc


subroutine partial_derivative_of_tor_shift_in_mc2(m,n,tor_shift,dtheta,dpsi,tor_shift_psi,tor_shift_th,&
     & tor_shift_psi_left_bdry_minus_one)
  !calculate the partial derivative with respect to the magnetic cooordinates (psi,theta)
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  use magnetic_coordinates,only: tor_shift_mc_left_bdry_minus_one, xgrid, zgrid
  use splines, only: spline3ders
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: tor_shift(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: tor_shift_psi(m,n),tor_shift_th(m,n)
  real(p_),intent(out):: tor_shift_psi_left_bdry_minus_one(n)
  integer:: i,j
  real(p_) :: tmp(m,n), tmp0(n)

  do i = 1, m
     call spline3ders(xgrid, tor_shift(i,:), xgrid, tmp(i,:), tor_shift_psi(i,:), tmp(i,:))
  enddo

  do j = 1, n
     call spline3ders(zgrid, tor_shift(:,j), zgrid, tmp(:,j), tor_shift_th(:,j), tmp(:,j))
  enddo

  call spline3ders(xgrid, tor_shift_mc_left_bdry_minus_one(:), xgrid, &
       & tmp0(:), tor_shift_psi_left_bdry_minus_one(:), tmp0(:))


end subroutine partial_derivative_of_tor_shift_in_mc2



function jacobian_func(theta,radcor) result (z) 
  use constants,only: p_
  use magnetic_coordinates,only: mpol,nrad,zgrid,xgrid,jacobian !as input
  use interpolate_module,only: linear_2d_interpolate

  implicit none
  real(p_)::radcor,theta,z
  call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,theta,radcor,z)  !uniform 1darray is assumed
end function jacobian_func
