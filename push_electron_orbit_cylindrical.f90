
subroutine push_electron_cylindrical_single_particle(dtao,mu,r,z,phi,vpar)
  !input: initial condition of the orbit: mu,r,z,phi,vpar
  !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
  use constants,only:p_
  use constants,only:zero,one,two,one_half,three,six,twopi
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
  use constants,only:p_
  use constants,only:twopi
    use gk_module,only: charge_sign_e  
use magnetic_field, only : b,bphi,b_z,bz,b_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_r,bstar_e_parallel,bstar_e_parallelval
  
  

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=bstar_e_r(r,z,phi,vpar)/bstar_e_parallelval*vpar +&
       & charge_sign_e(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*&
       & (bphi(r,z)*b_z(r,z)-bz(r,z)*b_phi(r,z)/r)

end function r_e_dot


function z_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_
  use constants,only:twopi
    use gk_module,only: charge_sign_e
use magnetic_field, only : b,bphi,b_z,bz,b_phi, br,b_r
    implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_z,bstar_e_parallel,bstar_e_parallelval
  
 

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval= bstar_e_z(r,z,phi,vpar)/bstar_e_parallelval*vpar + &
       & charge_sign_e(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*(br(r,z)*b_phi(r,z)/r-bphi(r,z)*b_r(r,z))
end function z_e_dot

function phi_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_
  use constants,only:twopi
  use gk_module,only: charge_sign_e
  use magnetic_field, only : b,bphi,b_z,bz,b_phi,br,b_r
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_phi,bstar_e_parallel,bstar_e_parallelval
   

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=bstar_e_phi(r,z,phi,vpar)/bstar_e_parallelval*vpar+&
       & charge_sign_e(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*(bz(r,z)*b_r(r,z)-br(r,z)*b_z(r,z))
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
    use gk_module,only: charge_sign_e
use magnetic_field, only : br,unitbphi_z,unitbz_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  

  funcval=br(r,z)+charge_sign_e(1)*vpar/twopi*(unitbz_phi(r,z)/r-unitbphi_z(r,z))

end function bstar_e_r


function bstar_e_z(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only:charge_sign_e
use magnetic_field, only : b,bphi,bz,unitbphi_r,unitbr_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_):: unitbphi
  

  unitbphi=bphi(r,z)/b(r,z)
  funcval=bz(r,z)+charge_sign_e(1)*vpar/twopi*(unitbphi_r(r,z)+unitbphi/r-unitbr_phi(r,z)/r)

end function bstar_e_z

function bstar_e_phi(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_e
use magnetic_field, only : bphi,unitbr_z,unitbz_r
    implicit none
  real(p_)::funcval,r,z,phi,vpar
  

  funcval=bphi(r,z)+charge_sign_e(1)*vpar/twopi*(unitbr_z(r,z)-unitbz_r(r,z))

end function bstar_e_phi


function bstar_e_parallel(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_e
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
  funcval=bval+charge_sign_e(1)*vpar/twopi*unitb_dot_curl_unitb
end function bstar_e_parallel

