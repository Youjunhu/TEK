module table_in_mc !used in computing guiding-center drift
  use constants,only:p_
  implicit none
  save
  real(p_), dimension(:,:), allocatable :: br_mc, bz_mc, bphi_mc, &
       & bp_mc, b_mc, bdgxcgy, & !bdgxcgy=B0_dot_grad_x_cross_grad_y
       & w1,w2,w3,w4,w5,w5p,w6,w7,w8,w8p,w9,w10,w12,w13,w14,w15

contains

  subroutine prepare_table_in_mc()
    use constants,only: one
    use domain_decomposition,only: myid
    use magnetic_coordinates,only: mpol, nflux, r_mc,z_mc,jacobian,radcor_1d_array, &
         &  GSpsi_prime, grad_psi, grad_alpha, &
         & grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, grad_alpha_dot_grad_theta, &
         & grad_psi_r, grad_psi_z, grad_theta_r, grad_theta_z, &
         & grad_alpha_r, grad_alpha_z
    use magnetic_field,only: br, bz, bphi, b, &
         & b_r, b_z, unitbr_z, unitbz_r, unitbphi_r,unitbphi_z
    implicit none

    real(p_) :: brval,bzval,bphival,bval
    real(p_) :: b_rval,b_zval
    real(p_) :: dradial_dr_func,dradial_dz_func !function names
    real(p_) :: dradial_dr_func2,dradial_dz_func2 !function names
    real(p_) :: dtheta_dr_func,dtheta_dz_func !function names
    real(p_) :: ddelta_dr_func,ddelta_dz_func,ddelta_dr_lsf_midplane_twopi,ddelta_dz_lsf_midplane_twopi !function names
    real(p_) :: dradial_dr_val,dradial_dz_val
    real(p_) :: dtheta_dr_val,dtheta_dz_val
    real(p_) :: ddelta_dr_val,ddelta_dz_val
    real(p_) :: unitbr,unitbz,unitbphi
    real(p_) :: curl_unitb_rcomp,curl_unitb_zcomp,curl_unitb_phicomp
    real(p_) :: unitb_dot_curl_unitb
    real(p_) :: r,z,radcor
    real(p_) :: grad_psi_val,grad_alpha_val, grad_psi_dot_grad_alpha_val
    real(p_) :: grad_psi_dot_grad_theta_val,grad_alpha_dot_grad_theta_val
    real(p_) :: dalpha_dr_val,dalpha_dz_val,dalpha_dphi_val
    integer :: i,j

    allocate(w1(mpol,nflux))
    allocate(w2(mpol,nflux))
    allocate(w3(mpol,nflux))
    allocate(w4(mpol,nflux))
    allocate(w5(mpol,nflux))
    allocate(w5p(mpol,nflux))
    allocate(w6(mpol,nflux))
    allocate(w7(mpol,nflux))
    allocate(w8(mpol,nflux))
    allocate(w8p(mpol,nflux))
    allocate(w9(mpol,nflux))
    allocate(w10(mpol,nflux))
    allocate(w12(mpol,nflux))
    allocate(w13(mpol,nflux))
    allocate(w14(mpol,nflux))
    allocate(w15(mpol,nflux))
    allocate(b_mc(mpol,nflux))
    allocate(br_mc(mpol,nflux))
    allocate(bz_mc(mpol,nflux))
    allocate(bphi_mc(mpol,nflux))
    allocate(bp_mc(mpol,nflux))
    allocate(bdgxcgy(mpol,nflux))


    do j=1,nflux
       radcor=radcor_1d_array(j)
       do i=1,mpol
          r=r_mc(i,j)
          z=z_mc(i,j)
          brval=br(r,z)
          bzval=bz(r,z)
          bphival=bphi(r,z)
          bval=sqrt(brval**2+bzval**2+bphival**2)
          b_zval=b_z(r,z)
          b_rval=b_r(r,z)
          unitbr=brval/bval
          unitbz=bzval/bval
          unitbphi=bphival/bval
          curl_unitb_rcomp=-unitbphi_z(r,z)
          curl_unitb_phicomp=unitbr_z(r,z)-unitbz_r(r,z)
          curl_unitb_zcomp=unitbphi_r(r,z)+unitbphi/r

          unitb_dot_curl_unitb=unitbr*curl_unitb_rcomp +unitbphi*curl_unitb_phicomp&
               & +unitbz*curl_unitb_zcomp

          dradial_dr_val=grad_psi_r(i,j)
          dradial_dz_val=grad_psi_z(i,j) 
          dtheta_dr_val=grad_theta_r(i,j)
          dtheta_dz_val=grad_theta_z(i,j)
          ddelta_dr_val = -grad_alpha_r(i,j)
          ddelta_dz_val = -grad_alpha_z(i,j)

          dalpha_dr_val = grad_alpha_r(i,j)
          dalpha_dz_val = grad_alpha_z(i,j)
          dalpha_dphi_val=one

          !        if(i.eq.mpol)    ddelta_dr_val=ddelta_dr_lsf_midplane_twopi(r) !at theta=twopi cut
          !        if(i.eq.mpol)    ddelta_dz_val=ddelta_dz_lsf_midplane_twopi(r) !at theta=twopi cut

!!$!write(*,*) 'i,j=',dradial_dr_val, dradial_dr_func2(r,z),dradial_dz_val, dradial_dz_func2(r,z)
          w1(i,j)=unitb_dot_curl_unitb/bval
          w2(i,j)=unitbr*dtheta_dr_val+unitbz*dtheta_dz_val
          !      w2(i,j)=-(psi_lcfs-psi_axis)/(bval*jacobian(i,j))
          w3(i,j)=curl_unitb_rcomp/bval*dradial_dr_val+curl_unitb_zcomp/bval*dradial_dz_val
          w4(i,j)=curl_unitb_rcomp/bval*dtheta_dr_val+curl_unitb_zcomp/bval*dtheta_dz_val
          w5(i,j)=curl_unitb_phicomp/(bval*r)-curl_unitb_rcomp/bval*ddelta_dr_val-curl_unitb_zcomp/bval*ddelta_dz_val
          w5p(i,j)=curl_unitb_phicomp/(bval*r)
          w6(i,j)=(bphival*b_zval)*dradial_dr_val+(-bphival*b_rval)*dradial_dz_val
          w6(i,j)=w6(i,j)/bval**2
          w7(i,j)=(bphival*b_zval)*dtheta_dr_val+(-bphival*b_rval)*dtheta_dz_val
          w7(i,j)=w7(i,j)/bval**2
          w8(i,j)=(bzval*b_rval-brval*b_zval)/r-(bphival*b_zval)*ddelta_dr_val-(-bphival*b_rval)*ddelta_dz_val
          w8(i,j)=w8(i,j)/bval**2
          w8p(i,j)=(bzval*b_rval-brval*b_zval)/(r*bval**2)
          w9(i,j)=unitbr*b_rval+unitbz*b_zval
          w10(i,j)=(b_rval*curl_unitb_rcomp+b_zval*curl_unitb_zcomp)/bval

!!$        grad_psi(i,j)=sqrt(dradial_dr_val**2+dradial_dz_val**2) !now calculated in another subroutine
!!$        grad_alpha(i,j)=sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)
!!$        grad_psi_dot_grad_alpha(i,j)=dradial_dr_val*ddelta_dr_val+dradial_dz_val*ddelta_dz_val

          b_mc(i,j)=bval
          br_mc(i,j)=brval
          bz_mc(i,j)=bzval
          bphi_mc(i,j)=bphival
          bp_mc(i,j)=sqrt(brval**2+bzval**2)

          grad_psi_val=grad_psi(i,j)
          grad_alpha_val=grad_alpha(i,j)
          grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha(i,j)
          grad_psi_dot_grad_theta_val=grad_psi_dot_grad_theta(i,j)
          grad_alpha_dot_grad_theta_val=grad_alpha_dot_grad_theta(i,j)

          w12(i,j)=GSpsi_prime/bval*(grad_psi_dot_grad_theta_val*grad_psi_dot_grad_alpha_val &
               & -grad_psi_val**2*grad_alpha_dot_grad_theta_val)
          w13(i,j)=GSpsi_prime/bval*(grad_psi_dot_grad_theta_val*grad_alpha_val**2 &
               & -grad_psi_dot_grad_alpha_val*grad_alpha_dot_grad_theta_val)

          w14(i,j) = bphival*(dradial_dz_val*dtheta_dr_val-dradial_dr_val*dtheta_dz_val)
          w15(i,j)=brval*(dtheta_dz_val)/r &
               & + bphival*(dalpha_dz_val*dtheta_dr_val-dalpha_dr_val*dtheta_dz_val) &
               & + bzval*(-dtheta_dr_val)/r
          bdgxcgy(i,j) = bval**2/GSpsi_prime
       enddo

    enddo

    if(myid.eq.0) call diagnostic()
  contains
    subroutine diagnostic()
      character(8)::filename
      integer:: file_unit,u
      open(newunit=u,file='grad_theta_alpha.txt')
      do j=1,nflux
         radcor=radcor_1d_array(j)
         do i=1,mpol
!!$          write(u,*) r_mc(i,j), z_mc(i,j), dtheta_dr_val,dtheta_dz_val,&
!!$               & sqrt(dtheta_dz_val**2+dtheta_dr_val**2), ddelta_dr_val,  ddelta_dz_val,&
!!$               sqrt(ddelta_dr_val**2+ddelta_dz_val**2), sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)
!!$          write(u,*) r_mc(i,j), z_mc(i,j),jacobian(i,j),&
!!$               & -w2(i,j)

            !write(u,*) r_mc(i,j), z_mc(i,j),w12(i,j),w13(i,j),w5(i,j)
            write(u,*) r_mc(i,j), z_mc(i,j), w5(i,j), w8(i,j)

            !               & -b_mc(i,j)/GSpsi_prime*w2(i,j)


         enddo
         write(u,*) 
      enddo
      close(u)
!!$  write(filename,'(a4,i4.4)') 'metr',myid
!!$  file_unit=myid+311
!!$  open(file_unit,file=filename)
!!$  do i=1,mpol
!!$     do j=1,nflux
!!$        write(file_unit,'(2i8.4,3(1pe14.5))')  i,j,grad_alpha(i,j),grad_psi_dot_grad_alpha(i,j), grad_psi(i,j)
!!$     enddo
!!$     write(file_unit,*)
!!$  enddo
!!$  close(file_unit)

      ! write(*,*) 'min max w1-4=', maxval(w1),minval(w1),maxval(w2),minval(w2),maxval(w3),minval(w3),maxval(w4),minval(w4)
      !write(*,*) 'min max w5-8=', maxval(w5),minval(w5),maxval(w6),minval(w6),maxval(w7),minval(w7),maxval(w8),minval(w8)
      !write(*,*) 'min max w9-10=', maxval(w9),minval(w9),maxval(w10),minval(w10)

    end subroutine diagnostic

  end subroutine prepare_table_in_mc

end module table_in_mc



module func_in_mc
contains
  function b_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: b_mc
    use magnetic_coordinates,only:mpol,nflux,theta_1d_array,radcor_1d_array
    use interpolate_module,only: linear_2d_interpolation

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,b_mc,theta,radcor,f)
  end function b_mc_func

  function br_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: br_mc
    use magnetic_coordinates,only:mpol,nflux,theta_1d_array,radcor_1d_array
    use interpolate_module,only: linear_2d_interpolation

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,br_mc,theta,radcor,f)
  end function br_mc_func

  function bz_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: bz_mc
    use magnetic_coordinates,only:mpol,nflux,theta_1d_array,radcor_1d_array
    use interpolate_module,only: linear_2d_interpolation

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,bz_mc,theta,radcor,f)
  end function bz_mc_func


  function bphi_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: bphi_mc
    use magnetic_coordinates,only:mpol,nflux,theta_1d_array,radcor_1d_array
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,bphi_mc,theta,radcor,f)
  end function bphi_mc_func

  function tor_shift_func(theta, radcor) result(z)
    use constants,only:p_, pi
    use magnetic_coordinates,only:mpol,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc
    use interpolate_module, only : linear_2d_interpolation
    implicit none
    real(p_) :: theta, radcor,z
    call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc,theta,radcor,z) !interpolating in magnetic coordinates to get tor_shift
  end function tor_shift_func

end module func_in_mc
