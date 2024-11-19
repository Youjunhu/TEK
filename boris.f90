module boris
contains

subroutine push_full_orbit_cylindrical_boris(charge_sign,dtao,r,phi,z,vr,vphi,vz)
  !input: initial condition of the orbit: r,phi,z,vr,vphi,vz
  !Output: the instanteous value of the orbit after dtao: r,phi,z,vr,vphi,vz
!in essence, the algorithm uses Cartesian coordinates, and then takes into account the rotation of the the basis vectors due to the change of particle's location
  use constants,only:p_
  use constants,only:zero,one,two,twopi
  use math,only: cross_product_in_cartesian
  use magnetic_field, only : br,bz,bphi 
  implicit none
  real(p_),intent(in):: charge_sign, dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of orbit

!  real(p_):: er,ez,ephi !function names
  real(p_):: tx,ty,tz,t,factor,sx,sy,sz
  real(p_):: vx_minus,vy_minus,vz_minus
  real(p_):: vx_prime,vy_prime,vz_prime
  real(p_):: vx_plus,vy_plus,vz_plus
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: alpha

  vx_minus=vr  +  er  (r,z,phi)*dtao/two*(twopi*charge_sign)
  vy_minus=vphi+  ephi(r,z,phi)*dtao/two*(twopi*charge_sign)
  vz_minus=vz  +  ez  (r,z,phi)*dtao/two*(twopi*charge_sign)

  tx=br(r,z)*dtao/two*(twopi*charge_sign)
  ty=bphi(r,z)*dtao/two*(twopi*charge_sign)
  tz=bz(r,z)*dtao/two*(twopi*charge_sign)

  call cross_product_in_cartesian(vx_minus,vy_minus,vz_minus,tx,ty,tz,cx,cy,cz)

  vx_prime=  vx_minus+ cx
  vy_prime=  vy_minus+ cy
  vz_prime=  vz_minus+ cz

  t=sqrt(tx*tx+tz*tz+ty*ty)
  factor=two/(one+t*t)
  sx=tx*factor
  sy=ty*factor
  sz=tz*factor

  call cross_product_in_cartesian(vx_prime,vy_prime,vz_prime,sx,sy,sz,cx,cy,cz)

  vx_plus =  vx_minus+ cx
  vy_plus =  vy_minus+ cy
  vz_plus =  vz_minus+ cz

  vx=vx_plus+  er  (r,z,phi)*dtao/two*(twopi*charge_sign)
  vy=vy_plus+  ephi(r,z,phi)*dtao/two*(twopi*charge_sign)
  vz=vz_plus+  ez  (r,z,phi)*dtao/two*(twopi*charge_sign)


  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  alpha=asin(y/r)

  phi=phi+alpha

  vr=cos(alpha)*vx+sin(alpha)*vy
  vphi=-sin(alpha)*vx+cos(alpha)*vy
end subroutine push_full_orbit_cylindrical_boris


subroutine push_full_orbit_cylindrical_boris_with_additional_input_output(dtao,radcor,theta,alpha, active,&
     & r,phi,z,vr,vphi,vz,vr_mid,vphi_mid,vz_mid)
  !derived from "push_full_orbit_cylindrical_boris", the modification is that: additional input phi_mid, which is value of phi at t_{n+1/2}, and additional output: vr_mid,vphi_mid,vz_mid, which are projection of velocity at t_{n+1/2} onto the basis vector at t_{n+1/2}, instead of t_{n+1}. Here, for general cases, n can also be an half-integer
  use constants,only:p_
  use constants,only:zero,one,two,twopi
  use control_parameters,only: fk_nonlinear
  use magnetic_field, only : br,bz,bphi 
  use math,only: cross_product_in_cartesian
  use force
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of (r,phi,z)  at t_{n} (as input) and t_{n+1} (as output), (vr,vphi,vz) is the projection of velocity at t_{n-1/2} at basis vector at t_{n} (as input); and the projection of velocity at t_{n+1/2} on the basis vector at t_{n+1} (as output}
  !  real(p_),intent(in):: phi_mid !value of phi at t_{n+1/2}, not used now, instead this value is calculated direclty in this subroutine
  real(p_),intent(in):: radcor,theta,alpha !perturbed field is interpolated in magnetic coordinates, so we need the magnetic coordinates of particles
  logical,intent(in):: active
  real(p_),intent(out):: vr_mid,vphi_mid,vz_mid !projection velocity at t_{n+1/2} on the basis vector at t_{n+1/2}
  !real(p_):: er,ez,ephi !function names
  real(p_):: er_val,ez_val,ephi_val
  real(p_):: tx,ty,tz,t,factor,sx,sy,sz
  real(p_):: vx_minus,vy_minus,vz_minus
  real(p_):: vx_prime,vy_prime,vz_prime
  real(p_):: vx_plus,vy_plus,vz_plus
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: dphi,phi0


  if(fk_nonlinear.eq.1) then
     call  field_perturbation_on_marker2(radcor,theta,alpha, active, er_val,ez_val,ephi_val) 
  else
     er_val=0._p_
     ez_val=0._p_
     ephi_val=0._p_
  endif
  phi0=phi !record the old value of phi
  vx_minus=vr   +  er_val*dtao/two*(twopi)
  vy_minus=vphi +  ephi_val*dtao/two*(twopi)
  vz_minus=vz   +  ez_val*dtao/two*(twopi)

  tx=br(r,z)*dtao/two*(twopi)
  ty=bphi(r,z)*dtao/two*(twopi)
  tz=bz(r,z)*dtao/two*(twopi)

  call cross_product_in_cartesian(vx_minus,vy_minus,vz_minus,tx,ty,tz,cx,cy,cz)

  vx_prime=  vx_minus+ cx
  vy_prime=  vy_minus+ cy
  vz_prime=  vz_minus+ cz

  t=sqrt(tx*tx+tz*tz+ty*ty)
  factor=two/(one+t*t)
  sx=tx*factor
  sy=ty*factor
  sz=tz*factor

  call cross_product_in_cartesian(vx_prime,vy_prime,vz_prime,sx,sy,sz,cx,cy,cz)

  vx_plus =  vx_minus+ cx
  vy_plus =  vy_minus+ cy
  vz_plus =  vz_minus+ cz

  vx=vx_plus +  er_val*dtao/two*(twopi)
  vy=vy_plus +  ephi_val*dtao/two*(twopi)
  vz=vz_plus +  ez_val*dtao/two*(twopi)

  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  dphi=asin(y/r)

  phi=phi+dphi

  vr=cos(dphi)*vx+sin(dphi)*vy !projected to basis vectors at the t_{n+1} particle location
  vphi=-sin(dphi)*vx+cos(dphi)*vy !projected to basis vectors at the t_{n+1} particle location

  !  dphi=phi_mid-phi0 !turn out not reliable because the toroidal angular momentum conservation is not as good as the following simple treatment
  dphi=(phi-phi0)/two !equivalent to dphi=dphi/two
  vr_mid=cos(dphi)*vx+sin(dphi)*vy !the projection of v_{n+1/2} on basis vectors at t_{n+1/2}
  vphi_mid=-sin(dphi)*vx+cos(dphi)*vy
  vz_mid=vz

end subroutine push_full_orbit_cylindrical_boris_with_additional_input_output


function er(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: er,r,z,phi

  er=0._p_

end function er


function ez(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: ez,r,z,phi

  ez=0._p_

end function ez


function ephi(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: ephi,r,z,phi

  ephi=0._p_

end function ephi





subroutine normalize_full_orbit_variables(nmarker,r,phi,z,vr,vphi,vz)
 use constants,only:p_
  use normalizing,only: Ln,vn_i

  implicit none
  integer,intent(in):: nmarker
  real(p_),intent(inout):: r(nmarker),phi(nmarker),z(nmarker)
  real(p_),intent(inout):: vr(nmarker),vphi(nmarker),vz(nmarker)

 !normalization
  r=r/Ln !convert to unit Ln
  z=z/Ln !convert to unit Ln
!  phi=phi
  vr=vr/vn_i !normalized by vn_i
  vphi=vphi/vn_i
  vz=vz/vn_i
end subroutine normalize_full_orbit_variables


subroutine particle_variables_to_guiding_center_variables(mass, charge, vn, r,phi,z,vr0,vphi0,vz0,rg,phig,zg, mu, vpar)
  use constants,only:p_,two
  use normalizing,only: bn
  use magnetic_field, only : br,bz,bphi 
  implicit none
  real(p_),intent(in):: mass, charge, vn
  real(p_),intent(in):: r,phi,z,vr0,vphi0,vz0 !particle variables, vr0, vphi0, vz0 in unit of vn
  real(p_),intent(out):: rg,phig,zg, mu, vpar !guiding-center variables

  real(p_):: brval,bphival,bzval,bval,factor
  real(p_):: vr,vz,vphi,v !in SI unit
  
  brval=br(r,z)
  bzval=bz(r,z)
  bphival=bphi(r,z)
  bval=sqrt(brval**2+bphival**2+bzval**2)
  vr=vr0*vn !to SI unit
  vz=vz0*vn !to SI unit
  vphi=vphi0*vn !to SI unit

  v=sqrt(vr**2+vz**2+vphi**2)

factor=mass/(bval**2*charge)

  !  rg=r+mass/(bval**2*charge)*(-vz*bphival+vphi*bzval) ! wrong!
  rg=sqrt((r+factor*(vphi*bzval-vz*bphival))**2+(factor*(vz*brval-vr*bzval))**2)
  !  phig=phi+atan(mass/(bval**2*charge)*(-vr*bzval+vz*brval)/r) !wrong!
  phig=phi+asin(factor*(-vr*bzval+vz*brval)/rg) !range of asin function is [-pi/2,pi/2]
  !zg=z+mass/(bval**2*charge)*(-vr*bphival-vphi*brval) !wrong, pointed by Yingfeng Xu
  zg=z+factor*(vr*bphival-vphi*brval) !corrected.

vpar=(vr*brval+vz*bzval+vphi*bphival)/bval
mu=(v**2-vpar**2)/(two*bval)
vpar=vpar/vn
mu=mu/(vn**2/bn)
end subroutine particle_variables_to_guiding_center_variables


subroutine guiding_center_variables_to_particle_variables(mass, charge, r_gc,z_gc,phi_gc, mu, vpar, r,z,phi, vr,vz,vphi) !see my notes for the formula
    use constants,only:p_
    use constants,only: twopi,two,one,four
    use normalizing,only:ln,bn
    use math,only: cross_product_in_cartesian
    use domain_decomposition,only: myid
    use magnetic_field, only : br,bz,bphi
    implicit none
    real(p_),intent(in):: mass, charge
    real(p_),intent(in):: r_gc,z_gc,phi_gc, mu, vpar
    real(p_),intent(out):: r,z,phi, vr,vz,vphi
    real(p_):: vn
    real(p_)::br_val,bz_val,bphi_val,bval,omega, gyro_radius
    real(p_):: mu_si, vpar_si,vper_si  !*_si indicates quanties in S.I. units
    real(p_):: root1,root2,a,b,c !coefficients of the quadratic equation
    real(p_):: vx,vy, vz2
    real(p_):: shiftx, shifty, shiftz
    integer,parameter:: niteration=2 !one iteration is usually enought
    integer:: k
    vn=ln/(twopi/(bn*charge/mass))
    mu_si=mu*(mass*vn**2/bn)
    vpar_si=vpar*vn

    r=r_gc !initial guess
    z=z_gc
    phi=phi_gc 
    do k=1, niteration
       br_val=br(r,z)
       bz_val=bz(r,z)
       bphi_val=bphi(r,z)
       bval=sqrt(br_val**2+bz_val**2+bphi_val**2)

       omega=bval*abs(charge)/mass
       vper_si=sqrt(mu_si*two*bval/mass)
!!$    gyro_radius=vper_si/omega
!!$    r=r_gc+gyro_radius
!!$    z=z_gc
!!$    phi=phi_gc

       vr=0._p_
       a=(bz_val/bphi_val)**2+one
       b=-two*bval*bz_val/bphi_val**2*vpar_si
       c=(bval*vpar_si/bphi_val)**2-two*bval*mu_si/mass-vpar_si**2
       root1=(-b+sqrt(b**2-four*a*c))/(two*a)
       root2=(-b-sqrt(b**2-four*a*c))/(two*a)
!!$       if((bphi_val*charge).gt.0) then
!!$          vz=root1/vn
!!$       else
!!$          vz=root2/vn
!!$       endif


       !          vz=root1/vn
       vz=root2/vn

       vphi=(bval*vpar-vz*bz_val)/bphi_val

       vx=vr*vn
       vy=vphi*vn
       vz2=vz*vn

       call cross_product_in_cartesian(vx,vy,vz2,br_val,bphi_val,bz_val,shiftx,shifty,shiftz)

       shiftx=-shiftx/(bval**2*charge/mass)
       shifty=-shifty/(bval**2*charge/mass)
       shiftz=-shiftz/(bval**2*charge/mass)

       r=r_gc + shiftx
       z=z_gc + shiftz
       phi=phi_gc+asin(shifty/r)
       !if(myid.eq.0)write(*,*) 'iteration=', k, 'r, z, phi=', r, z, phi
    enddo

!!$    write(*,*) 'v_squared=',vphi**2+vz**2, vpar**2+mu_si*two*bval/mass/vn**2 !to very they agree with each other
!!$    yj: block
!!$      real(p_):: rg2,phig2,zg2, mu2, vpar2
!!$      call particle_variables_to_guiding_center_variables(mass, charge, vn, r,phi,z,vr,vphi,vz,rg2,phig2,zg2, mu2, vpar2)
!!$     write(*,*) 'rg,phig,zg, mu, vpar=', rg2,r_gc, phig2, phi_gc, zg2, z_gc, mu2, mu, vpar2, vpar !to very they agree with each other
!!$    end block yj
  end subroutine guiding_center_variables_to_particle_variables


end module boris



module initial_half_step_for_boris
implicit none
contains
subroutine backward_half_step_for_boris(charge_sign,dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm
  !input:  r,z,phi,vr,vz,vphi, at the same time t_{0}
  !Output: the projection of velocity at t_{-1/2} onto the basis vectors at t_{0}
  use constants,only:p_
  use constants,only:one_half
  use constants,only:zero,one,two,one_half,three,six,twopi,kev
  implicit none
  real(p_),intent(in):: charge_sign,dtao
  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvr,dvz,dvphi
!  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot
!  real(p_):: vz_fo_dot !equations of motion
!  real(p_):: vr_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
  real(p_):: step,vr_new,vphi_new,r0,z0,phi0,dt
  integer,parameter:: m=100
  integer:: i

!write(*,*)  "energy_before evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
r0=r
phi0=phi
z0=z

step=-0.5_p_*dtao
dt=step/m
do i=1,m
  kr1=    one_half*dt*r_fo_dot(r,z,phi,vr,vz,vphi)
  kz1=    one_half*dt*z_fo_dot(r,z,phi,vr,vz,vphi)
  kphi1=  one_half*dt*phi_fo_dot(r,z,phi,vr,vz,vphi)
  kvr1=   one_half*dt*vr_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  kvphi1= one_half*dt*vphi_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  kvz1=   one_half*dt*vz_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)

  dr=   dt*r_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dz=   dt*z_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dphi= dt*phi_fo_dot  (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvr=  dt*vr_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvphi=dt*vphi_fo_dot (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvz=  dt*vz_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)

  z=z+dz
  r=r+dr
  phi=phi+dphi
  vz=vz+dvz
  vr=vr+dvr
  vphi=vphi+dvphi
enddo

dphi=phi-phi0
vr_new=vr
vphi_new=vphi

!write(*,*)  "energy_before projection=", 0.5_p_*mass*(vr_new**2+vz**2+vphi_new**2)*vn**2/kev
!projection of the new velocity onto the old basis vectors
vr=  -vphi_new*sin(dphi)+vr_new*cos(dphi) 
vphi=vphi_new*cos(dphi)+vr_new*sin(dphi) 
!write(*,*)  "energy_after projection=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
r=r0 !go back to the original value
z=z0 !go back to the original value
phi=phi0 !go back to the original value
!write(*,*) 'dvphi=',dvphi
end subroutine backward_half_step_for_boris


subroutine forward_half_step_for_boris(charge_sign,dtao,r0,z0,phi0,vr0,vz0,vphi0,r1,z1,phi1,vr1,vz1,vphi1) 
  !input: initial condition of the orbit: r0,z0,phi0,vr0,vz0,vphi0, at the same time t_{0}
  !Output: the instanteous value of (r,z,phi) after half_dtao, i.e., t_{1/2}, and the projection of velocity at t_{0} onto the basis vectors at t_{1/2}
  use constants,only:p_
  use constants,only:zero,one,two,one_half,twopi
    use math,only: shift_to_specified_toroidal_range
  implicit none
  real(p_),intent(in):: charge_sign,dtao
  real(p_),intent(in):: r0,z0,phi0,vr0,vz0,vphi0  !instantaneous value of orbit at t0
  real(p_),intent(out):: r1,z1,phi1 !location at t0+half_dtao
  real(p_),intent(out):: vr1,vz1,vphi1   !projection of the old velocity (t=t0) on the new basis vector (determined by the new (r,z,phi))
  real(p_):: step,dt,dr,dz,dphi,dvr,dvz,dvphi
!  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot,vr_fo_dot,vz_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
!  real(p_):: kr2,kz2,kphi2,kvr2,kvz2,kvphi2 !Runge-Kutta steps
  integer,parameter::m=100 ! if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple steps
  integer:: k
  real(p_):: r,z,phi,vr,vz,vphi !working variables

  r=r0
  z=z0
  phi=phi0
  vr=vr0
  vz=vz0
  vphi=vphi0

  step=0.5_p_*dtao
  dt=step/m
  do k=1,m
     !2nd order Rung-Kuta method
     kr1=     one_half*dt*r_fo_dot(r,z,phi,vr,vz,vphi)
     kz1=     one_half*dt*z_fo_dot(r,z,phi,vr,vz,vphi)
     kphi1=   one_half*dt*phi_fo_dot(r,z,phi,vr,vz,vphi)
     kvr1=    one_half*dt*vr_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
     kvz1=    one_half*dt*vz_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
     kvphi1=  one_half*dt*vphi_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)

     dr=    dt*r_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dz=    dt*z_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dphi=  dt*phi_fo_dot  (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvr=   dt*vr_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvz=   dt*vz_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvphi= dt*vphi_fo_dot (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)

     !update
     r=r+      dr
     z=z+      dz
     phi=phi+  dphi
     vr=vr+dvr
     vz=vz+dvz
     vphi=vphi+dvphi
     !write(*,*) 'dvphi=',dvphi
  enddo

r1=r
z1=z
phi1=phi

dphi=phi1-phi0

vz1=vz0
vr1=  vphi0*sin(dphi)+vr0*cos(dphi) !projection of the old velocity on the new basis vector 
vphi1=vphi0*cos(dphi)-vr0*sin(dphi) !projection of the old velocity on the new basis vector 

!!$phi1=phi1-int(phi1/twopi)*twopi !shift into the range [0:twopi]
!!$if(phi1.lt.0) phi1=phi1+twopi !shift into the range [0:twopi]
!   call shift_to_zero_twopi_range(phi1)
   call shift_to_specified_toroidal_range(phi1)
end subroutine forward_half_step_for_boris



function r_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: r_fo_dot,r,z,phi,vr,vz,vphi

  r_fo_dot=vr
end function r_fo_dot

function phi_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: phi_fo_dot,r,z,phi,vr,vz,vphi

  phi_fo_dot=vphi/r
end function phi_fo_dot

function z_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: z_fo_dot,r,z,phi,vr,vz,vphi

  z_fo_dot=vz
end function z_fo_dot

function vr_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  use magnetic_field, only : bphi,bz
  implicit none
  real(p_):: charge_sign,vr_fo_dot,r,z,phi,vr,vz,vphi
 
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi)) !wrong, term vphi**2/r**2 is missing
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi))+ vphi**2/r**2 !2017-Aug.12, a bug found, where r**2 should be replaced by r. kinetic energy is not well conserved by the buggy code (compared with results given by the Cartesian version of orbit integrator), which forced me to examine the codes to find possible bugs, and I finally found this bug. After correcting this code, the conservation of kinetic energy given by this code is as good as that of the Cartesian version.
  vr_fo_dot=twopi*charge_sign*(-vz*bphi(r,z)+vphi*bz(r,z))+ vphi**2/r !correct
end function vr_fo_dot



function vphi_fo_dot(charge_sign,r,z,phi,vr,vz,vphi) 
  use constants,only:p_
  use constants,only:twopi
  use magnetic_field, only : br,bz
  implicit none
  real(p_):: charge_sign,vphi_fo_dot,r,z,phi,vr,vz,vphi

  vphi_fo_dot=twopi*charge_sign*(-vr*bz(r,z)+vz*br(r,z))-vphi*vr/r !the last term is inertial force

end function vphi_fo_dot

function vz_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  use constants,only:p_,twopi
  use magnetic_field, only : br,bphi
  implicit none
  real(p_):: charge_sign,vz_fo_dot,r,z,phi,vr,vz,vphi
  vz_fo_dot=twopi*charge_sign*(vr*bphi(r,z)-vphi*br(r,z))
end function vz_fo_dot
end module initial_half_step_for_boris
