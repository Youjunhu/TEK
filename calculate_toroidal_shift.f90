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


subroutine calculate_toroidal_shift2(j, g_value, x_contour,z_contour,jacob, m, tor_shift, tor_shift_left_bdry_minus_one)
  !cf. calculate_toroidal_shift, the differences are (1) g_value instead of psival is provided as input;
  !(2) compute a series of tor_shift on a magnetic surface (rather than a single value). The range of theta is [-pi:pi]
  use constants,only: p_, one_half, two
  use magnetic_coordinates,only: dtheta, sign_of_jacobian,sign_of_GSpsi_prime,i_theta_zero, GSpsi_array, GSpsi_prime, qhat
  use domain_decomposition,only: dtheta2, multi_eq_cells
  use magnetic_field, only: psi_z_func, psi_r_func
  implicit none
  integer,intent(in) :: j, m
  real(p_),intent(in) :: g_value, x_contour(m), z_contour(m), jacob(m)
  real(p_),intent(out) :: tor_shift(m)
  real(p_),intent(out) :: tor_shift_left_bdry_minus_one
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
