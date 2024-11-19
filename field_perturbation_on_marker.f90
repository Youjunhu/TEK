module force
contains

!!$subroutine field_perturbation_on_marker(radcor,theta,alpha,active,&
!!$     & epar0,ex0,ey0,mf_par0,mf_x0,mf_y0)



!!$subroutine potential_on_marker(radcor,theta,alpha,active,potential_val)
!!$  use constants,only:p_
!!$  use constants,only: one,zero,twopi
!!$  use perturbation_field_matrix,only: potential_left,potential_right
!!$  use magnetic_coordinates,only: mtor,nflux2,tor_1d_array,radcor_1d_array2
!!$  use domain_decomposition,only: dtheta2,theta_start
!!$  use interpolate_module,only: linear_2d_interpolation
!!$  implicit none
!!$  real(p_),intent(in)::radcor,theta,alpha
!!$  logical,intent(in)::active
!!$  real(p_),intent(out)::potential_val
!!$  real(p_):: coeff1,coeff2,tmp1,tmp2
!!$
!!$ if(active.eqv..false.) then !force on particles outside the computational region is set to zero
!!$   potential_val=0._p_
!!$  else
!!$
!!$  coeff1=(theta-theta_start)/dtheta2
!!$  coeff2=one-coeff1
!!$  !if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
!!$  call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,potential_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
!!$  call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,potential_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
!!$  potential_val=tmp1*coeff2+tmp2*coeff1
!!$endif
!!$end subroutine potential_on_marker


!!$  subroutine epar_on_marker(radcor,theta,alpha,epar_val)
!!$    use constants,only:p_
!!$    use constants,only: one,zero,twopi
!!$    use perturbation_field_matrix,only: epar_left,epar_right
!!$    !use perturbation_field_matrix,only: epar_left=>epar_left_half,epar_right=>epar_right_half
!!$    use magnetic_coordinates,only: mtor,nflux2,tor_1d_array,radcor_1d_array2
!!$    use domain_decomposition,only: dtheta2,theta_start
!!$    use interpolate_module,only: linear_2d_interpolation
!!$    implicit none
!!$    real(p_),intent(in)::radcor,theta,alpha
!!$    real(p_),intent(out)::epar_val
!!$    real(p_):: coeff1,coeff2,tmp1,tmp2
!!$
!!$    coeff1=(theta-theta_start)/dtheta2
!!$    coeff2=one-coeff1
!!$
!!$    !if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
!!$    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,epar_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
!!$    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,epar_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
!!$    epar_val=tmp1*coeff2+tmp2*coeff1
!!$
!!$  end subroutine epar_on_marker


  subroutine field_perturbation_on_marker_for_adiabatic_e_model(radcor,theta,alpha,&
       & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val) !for adiabatic electron model
    use constants,only:p_
    use constants,only: one,zero,twopi
    use perturbation_field_matrix,only:ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
    use perturbation_field_matrix,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right !field at nearby (larger) poloidal location
    use magnetic_coordinates,only: mtor,nflux2,tor_1d_array,radcor_1d_array2
    use domain_decomposition,only: dtheta2,theta_start !as input
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_),intent(in)::radcor,theta,alpha
    real(p_),intent(out)::ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
    real(p_):: coeff1,coeff2,tmp1,tmp2

    coeff1=(theta-theta_start)/dtheta2
    coeff2=one-coeff1

    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_r_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_r_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
    ef_cyl_r_val=tmp1*coeff2+tmp2*coeff1

    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_z_left,alpha,radcor,tmp1)  
    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_z_right,alpha,radcor,tmp2)  
    ef_cyl_z_val=tmp1*coeff2+tmp2*coeff1

    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_phi_left,alpha,radcor,tmp1)
    call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_phi_right,alpha,radcor,tmp2)  
    ef_cyl_phi_val=tmp1*coeff2+tmp2*coeff1
  end subroutine field_perturbation_on_marker_for_adiabatic_e_model

  subroutine field_perturbation_on_marker2(radcor,theta,alpha,active,&
       & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val) !using cylindrical components
    use constants,only:p_
    use constants,only: one,zero,twopi
    use perturbation_field_matrix,only:ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
    use perturbation_field_matrix,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right !field at nearby (larger) poloidal location
    use magnetic_coordinates,only: mtor,nflux2,tor_1d_array,radcor_1d_array2
    use domain_decomposition,only: dtheta2,theta_start !as input
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_),intent(in)::radcor,theta,alpha
    logical,intent(in):: active
    real(p_),intent(out)::ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
    real(p_):: coeff1,coeff2,tmp1,tmp2

    if(active.eqv..false.) then !force on particles outside the computational region is set to zero
       ef_cyl_r_val=0._p_
       ef_cyl_z_val=0._p_
       ef_cyl_phi_val=0._p_
    else

       coeff1=(theta-theta_start)/dtheta2
       coeff2=one-coeff1

       call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_r_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
       call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_r_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
       ef_cyl_r_val=tmp1*coeff2+tmp2*coeff1

       call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_z_left,alpha,radcor,tmp1)  
       call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_z_right,alpha,radcor,tmp2)  
       ef_cyl_z_val=tmp1*coeff2+tmp2*coeff1

       call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_phi_left,alpha,radcor,tmp1)
       call linear_2d_interpolation(mtor+1,nflux2,tor_1d_array,radcor_1d_array2,ef_cyl_phi_right,alpha,radcor,tmp2)  
       ef_cyl_phi_val=tmp1*coeff2+tmp2*coeff1
    endif
  end subroutine field_perturbation_on_marker2
end module force


!!$subroutine average_em_field()
!!$  use constants,only:one,two
!!$  use perturbation_field_matrix
!!$  implicit none
!!$
!!$  epar_left_half= (epar_left+epar_left_old)/two
!!$  ex_left_half=  (ex_left+ex_left_old)/two
!!$  ey_left_half=  (ey_left+ey_left_old)/two
!!$
!!$  mf_par_left_half=  (mf_par_left+mf_par_left_old)/two
!!$  mf_x_left_half=  (mf_x_left+mf_x_left_old)/two
!!$  mf_y_left_half=  (mf_y_left+mf_y_left_old)/two
!!$
!!$  epar_right_half= (epar_right+epar_right_old)/two
!!$  ex_right_half=  (ex_right+ex_right_old)/two
!!$  ey_right_half=  (ey_right+ey_right_old)/two
!!$
!!$  mf_par_right_half=  (mf_par_right+mf_par_right_old)/two
!!$  mf_x_right_half=  (mf_x_right+mf_x_right_old)/two
!!$  mf_y_right_half=  (mf_y_right+mf_y_right_old)/two
!!$
!!$!  if (isnan(maxval(epar_left_half))) write(*,'(a40,30(1pe12.1))') 'warning**nan_in average_em_field', maxval(epar_left_half)
!!$!write(*,*) 'maxval(epar_left_half)=',maxval(epar_left_half)
!!$end subroutine average_em_field


!!$subroutine store_old_em_field()
!!$  use constants,only:two
!!$  use perturbation_field_matrix !,only:epar_left,ex_left,ey_left,epar_left_old,ex_left_old,ey_left_old,&
!!$!       &   mf_par_left,mf_x_left,mf_y_left,mf_par_left_old,mf_x_left_old,mf_y_left_old
!!$  implicit none
!!$
!!$potential_left_old=potential_left
!!$potential_right_old=potential_right
!!$
!!$  epar_left_old=epar_left
!!$  ex_left_old=ex_left
!!$  ey_left_old=ey_left
!!$  mf_par_left_old=mf_par_left
!!$  mf_x_left_old=mf_x_left
!!$  mf_y_left_old=mf_y_left
!!$
!!$
!!$  epar_right_old=epar_right
!!$  ex_right_old=ex_right
!!$  ey_right_old=ey_right
!!$  mf_par_right_old=mf_par_right
!!$  mf_x_right_old=mf_x_right
!!$  mf_y_right_old=mf_y_right
!!$
!!$!ppar_e_left_old=ppar_e_left
!!$end subroutine store_old_em_field
