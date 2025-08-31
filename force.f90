module force
contains

!!$subroutine field_perturbation_on_marker(radcor,theta,alpha,active,&
!!$     & epar0,ex0,ey0,mf_par0,mf_x0,mf_y0)



!!$subroutine potential_on_marker(radcor,theta,alpha,active,potential_val)
!!$  use constants,only:p_
!!$  use constants,only: one,zero,twopi
!!$  use perturbation_field,only: potential_left,potential_right
!!$  use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
!!$  use domain_decomposition,only: dtheta2,theta_start
!!$  use interpolate_module,only: linear_2d_interpolate
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
!!$  call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,potential_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
!!$  call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,potential_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
!!$  potential_val=tmp1*coeff2+tmp2*coeff1
!!$endif
!!$end subroutine potential_on_marker


!!$  subroutine epar_on_marker(radcor,theta,alpha,epar_val)
!!$    use constants,only:p_
!!$    use constants,only: one,zero,twopi
!!$    use perturbation_field,only: epar_left,epar_right
!!$    !use perturbation_field,only: epar_left=>epar_left_half,epar_right=>epar_right_half
!!$    use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
!!$    use domain_decomposition,only: dtheta2,theta_start
!!$    use interpolate_module,only: linear_2d_interpolate
!!$    implicit none
!!$    real(p_),intent(in)::radcor,theta,alpha
!!$    real(p_),intent(out)::epar_val
!!$    real(p_):: coeff1,coeff2,tmp1,tmp2
!!$
!!$    coeff1=(theta-theta_start)/dtheta2
!!$    coeff2=one-coeff1
!!$
!!$    !if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
!!$    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,epar_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
!!$    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,epar_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
!!$    epar_val=tmp1*coeff2+tmp2*coeff1
!!$
!!$  end subroutine epar_on_marker


  subroutine field_perturbation_on_marker_for_adiabatic_e_model(radcor,theta,alpha,&
       & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val) !for adiabatic electron model
    use constants,only:p_
    use constants,only: one,zero,twopi
    use perturbation_field,only:ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
    use perturbation_field,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right !field at nearby (larger) poloidal location
    use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
    use domain_decomposition,only: dtheta2,theta_start !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_),intent(in)::radcor,theta,alpha
    real(p_),intent(out)::ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
    real(p_):: coeff1,coeff2,tmp1,tmp2

    coeff1=(theta-theta_start)/dtheta2
    coeff2=one-coeff1

    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
    ef_cyl_r_val=tmp1*coeff2+tmp2*coeff1

    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_left,alpha,radcor,tmp1)  
    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_right,alpha,radcor,tmp2)  
    ef_cyl_z_val=tmp1*coeff2+tmp2*coeff1

    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_left,alpha,radcor,tmp1)
    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_right,alpha,radcor,tmp2)  
    ef_cyl_phi_val=tmp1*coeff2+tmp2*coeff1
  end subroutine field_perturbation_on_marker_for_adiabatic_e_model

  subroutine field_perturbation_on_marker2(radcor,theta,alpha,active,&
       & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val) !using cylindrical components
    use constants,only:p_
    use constants,only: one,zero,twopi
    use perturbation_field,only:ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
    use perturbation_field,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right !field at nearby (larger) poloidal location
    use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
    use domain_decomposition,only: dtheta2,theta_start !as input
    use interpolate_module,only: linear_2d_interpolate
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

       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
       ef_cyl_r_val=tmp1*coeff2+tmp2*coeff1

       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_left,alpha,radcor,tmp1)  
       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_right,alpha,radcor,tmp2)  
       ef_cyl_z_val=tmp1*coeff2+tmp2*coeff1

       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_left,alpha,radcor,tmp1)
       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_right,alpha,radcor,tmp2)  
       ef_cyl_phi_val=tmp1*coeff2+tmp2*coeff1
    endif
  end subroutine field_perturbation_on_marker2
end module force



subroutine push_electron_cylindrical_single_particle(dtao,mu,r,z,phi,vpar)
  !input: initial condition of the orbit: mu,r,z,phi,vpar
  !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
  use constants,only:p_, zero,one,two,one_half,three,six,twopi
  implicit none
  real(p_),intent(in):: dtao,mu
  real(p_),intent(inout):: r,z,phi,vpar  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvpar
  real(p_):: r_e_dot,z_e_dot,phi_e_dot,vpar_e_dot !equations of motion
  real(p_):: kr1,kz1,kvpar1,kphi1,kr2,kz2,kvpar2,kphi2,kr3,kz3,kvpar3,kphi3,kr4,kz4,kvpar4,kphi4 !Runge-Kutta steps

  !4nd order Rung-Kutta method
  kr1=dtao*r_e_dot(r,z,phi,vpar,mu)
  kz1=dtao*z_e_dot(r,z,phi,vpar,mu)
  kvpar1=dtao*vpar_e_dot(r,z,phi,vpar,mu)
  kphi1=dtao*phi_e_dot(r,z,phi,vpar,mu)

  kr2=dtao*r_e_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kz2=dtao*z_e_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kvpar2=dtao*vpar_e_dot(r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kphi2=dtao*phi_e_dot  (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)

  kr3=dtao*r_e_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kz3=dtao*z_e_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kvpar3=dtao*vpar_e_dot(r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kphi3=dtao*phi_e_dot  (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)

  kr4=dtao*r_e_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kz4=dtao*z_e_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kvpar4=dtao*vpar_e_dot(r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kphi4=dtao*phi_e_dot  (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)

  dr=kr1/six+kr2/three+kr3/three+kr4/six
  dz=kz1/six+kz2/three+kz3/three+kz4/six
  dphi=kphi1/six+kphi2/three+kphi3/three+kphi4/six
  dvpar=kvpar1/six+kvpar2/three+kvpar3/three+kvpar4/six
  !update
  r=r+      dr
  z=z+      dz
  vpar=vpar+dvpar
  phi=phi+  dphi
  !write(*,*) 'dvpar=',dvpar
end subroutine push_electron_cylindrical_single_particle



function r_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_, twopi
    use gk_module,only: charge_sign_gk  
use magnetic_field, only : b,bphi,b_z,bz,b_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_r,bstar_e_parallel,bstar_e_parallelval
  
  

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=bstar_e_r(r,z,phi,vpar)/bstar_e_parallelval*vpar +&
       & charge_sign_gk(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*&
       & (bphi(r,z)*b_z(r,z)-bz(r,z)*b_phi(r,z)/r)

end function r_e_dot


function z_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_, twopi
    use gk_module,only: charge_sign_gk
use magnetic_field, only : b,bphi,b_z,bz,b_phi, br,b_r
    implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_z,bstar_e_parallel,bstar_e_parallelval
  
 

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval= bstar_e_z(r,z,phi,vpar)/bstar_e_parallelval*vpar + &
       & charge_sign_gk(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*(br(r,z)*b_phi(r,z)/r-bphi(r,z)*b_r(r,z))
end function z_e_dot

function phi_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_
  use constants,only:twopi
  use gk_module,only: charge_sign_gk
  use magnetic_field, only : b,bphi,b_z,bz,b_phi,br,b_r
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_phi,bstar_e_parallel,bstar_e_parallelval
   

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=bstar_e_phi(r,z,phi,vpar)/bstar_e_parallelval*vpar+&
       & charge_sign_gk(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*(bz(r,z)*b_r(r,z)-br(r,z)*b_z(r,z))
  funcval=funcval/r
end function phi_e_dot

function vpar_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_
  use constants,only:twopi
use magnetic_field, only : b,bphi,b_z,bz,b_phi, b_r
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_r,bstar_e_z,bstar_e_phi,bstar_e_parallel,bstar_e_parallelval
  

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=-mu*(bstar_e_r(r,z,phi,vpar)/bstar_e_parallelval*b_r(r,z)+bstar_e_z(r,z,phi,vpar)/bstar_e_parallelval*b_z(r,z) &
       & +bstar_e_phi(r,z,phi,vpar)/bstar_e_parallelval*b_phi(r,z)/r )
end function vpar_e_dot


function bstar_e_r(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_gk
use magnetic_field, only : br,unitbphi_z,unitbz_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  

  funcval=br(r,z)+charge_sign_gk(1)*vpar/twopi*(unitbz_phi(r,z)/r-unitbphi_z(r,z))

end function bstar_e_r


function bstar_e_z(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only:charge_sign_gk
use magnetic_field, only : b,bphi,bz,unitbphi_r,unitbr_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_):: unitbphi
  

  unitbphi=bphi(r,z)/b(r,z)
  funcval=bz(r,z)+charge_sign_gk(1)*vpar/twopi*(unitbphi_r(r,z)+unitbphi/r-unitbr_phi(r,z)/r)

end function bstar_e_z

function bstar_e_phi(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_gk
use magnetic_field, only : bphi,unitbr_z,unitbz_r
    implicit none
  real(p_)::funcval,r,z,phi,vpar
  

  funcval=bphi(r,z)+charge_sign_gk(1)*vpar/twopi*(unitbr_z(r,z)-unitbz_r(r,z))

end function bstar_e_phi


function bstar_e_parallel(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_gk
use magnetic_field, only : b,br,bz,bphi, unitbr_z,unitbz_r,unitbphi_r,unitbphi_z, unitbr_phi,unitbz_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_)::unitb_dot_curl_unitb
  real(p_):: unitbr,unitbphi,unitbz,bval
  
  bval=b(r,z)
  unitbr=br(r,z)/bval
  unitbz=bz(r,z)/bval
  unitbphi=bphi(r,z)/bval

  unitb_dot_curl_unitb=unitbr*(unitbz_phi(r,z)/r-unitbphi_z(r,z))&
       & +unitbphi*(unitbr_z(r,z)-unitbz_r(r,z))&
       & +unitbz*(unitbphi_r(r,z)+unitbphi/r-unitbr_phi(r,z)/r)
  funcval=bval+charge_sign_gk(1)*vpar/twopi*unitb_dot_curl_unitb
end function bstar_e_parallel

