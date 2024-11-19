module push_ion_weight_module
contains
subroutine push_ion_weight(dtao,nmarker_i,active_i,radcor_i,theta_i,alpha_i, &
     & grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i,&
     & v_i, vpar_i,vx_i,vy_i,w_i,w_i_star)
  use constants,only:p_
  use constants,only: two,twopi,one_half,kev
  use normalizing,only: vn_i,Ln
  use fk_module,only: ni0,mass_i,ps_vol_i,normalizing_factor
!  use force,only: field_perturbation_on_marker
  use magnetic_field,only: minor_r_prime
  use magnetic_coordinates, only : GSpsi_prime
  use density_temperature_profile_mod,only: kappa_ti_func, kappa_ni_func, ti_func
  implicit none
  integer,intent(in):: nmarker_i
  logical,intent(in):: active_i(nmarker_i)
  real(p_),intent(in):: dtao, radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  real(p_),intent(in):: grad_psi_i(nmarker_i),grad_alpha_i(nmarker_i),grad_psi_dot_grad_alpha_i(nmarker_i),&
       & bval_i(nmarker_i), v_i(nmarker_i), vpar_i(nmarker_i),&
       & vx_i(nmarker_i),vy_i(nmarker_i),w_i(nmarker_i)
  real(p_),intent(out):: w_i_star(nmarker_i)
  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: tmp,gradient_term
  integer:: k

  !tmp=ti0*kev/(mass_i*vn_i**2)
  !eq_particle_number=ni0*sqrt((mass_i/(twopi*ti0*kev))**3)*Ln**3*vn_i**3/normalizing_factor !explicitly cancel the exp(-v^2/vt^2) dependence, valid only for the case that both physical equilibrium distribution and marker distribution are Maxwellian.

!$omp parallel do private(gradient_term,ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val,wiprime)
  do k=1,nmarker_i
     if( active_i(k).eqv..false.) then
        w_i_star(k)=0._p_
     else
        tmp=ti_func(radcor_i(k))*kev/(mass_i*vn_i**2)
        eq_particle_number=ps_vol_i(k)*fi0(radcor_i(k),v_i(k))
        gradient_term=kappa_ni_func(radcor_i(k))+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti_func(radcor_i(k))
 !       call field_perturbation_on_marker(radcor_i(k),theta_i(k),alpha_i(k),active_i(k),&
        !            & ef_par_val,ef_x_val,ef_y_val)
        ef_x_val=0 !testing
        ef_y_val=0
        ef_par_val=0
        wiprime=eq_particle_number*twopi/tmp*(ef_par_val*vpar_i(k)+ef_x_val*vx_i(k)+ef_y_val*vy_i(k) ) & 
             &+eq_particle_number*gradient_term*GSpsi_prime*minor_r_prime(radcor_i(k))/bval_i(k)**2 &
             !& *(grad_psi_i(k)**2*grad_alpha_i(k)**2-grad_psi_dot_grad_alpha_i(k)**2)*ef_y_val &
             & *(bval_i(k)**2/GSpsi_prime**2)*ef_y_val !&
!em             &-eq_particle_number*gradient_term*minor_r_prime(radcor_i(k))/bval_i(k)*(vx_i(k)*mf_par_val&
!em             & -vpar_i(k)*(mf_x_val*grad_psi_i(k)**2+mf_y_val*grad_psi_dot_grad_alpha_i(k)))

        w_i_star(k)=w_i(k)+wiprime*dtao
        !if (isnan(wiprime)) write(*,'(a20,30(1pe12.1))') 'warning**nan appear', w_i_star(k),wiprime
          !write(*,*) ef_par_val
     endif
  enddo
!$omp end parallel do
!w_i_star=0._p_ !for testing
        !write(*,*) 'max wi=',maxval(abs(w_i)),'max wi_star', maxval(abs(w_i_star))
end subroutine push_ion_weight

subroutine push_ion_weight_using_electric_cylindrical_components(dtao,radcor_i,theta_i,alpha_i,&
     & r_i,z_i,phi_i,vr_i,vz_i,vphi_i,active_i,w_i,w_i_star,nmarker_i)
  use constants,only:p_
  use constants,only: two,twopi,one_half,kev
  use normalizing,only: vn_i,ln
  use fk_module,only:ps_vol_i,mass_i,v_i,ni0 ,normalizing_factor
  use magnetic_coordinates,only:mpol,nflux,theta_1d_array,radcor_1d_array, &
       & grad_psi_r,grad_psi_z
  use magnetic_field, only : br,bz,bphi
  use force
  use interpolate_module
  use magnetic_field,only: minor_r_prime
  use density_temperature_profile_mod,only: kappa_ti_func, kappa_ni_func, ti_func
  implicit none
  real(p_),intent(in):: dtao
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i),w_i(nmarker_i)
  real(p_),intent(in):: r_i(nmarker_i),z_i(nmarker_i),phi_i(nmarker_i)
  real(p_),intent(in):: vr_i(nmarker_i),vz_i(nmarker_i),vphi_i(nmarker_i)
  logical,intent(in):: active_i(nmarker_i)
  real(p_),intent(out):: w_i_star(nmarker_i)

  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
  real(p_):: tmp,tmp2

  real(p_):: br_val,bz_val,bphi_val, grad_psi_r0,grad_psi_z0,bval
  integer:: k,ierr

  do k=1,nmarker_i
     v_i(k)=sqrt(vr_i(k)**2+vz_i(k)**2+vphi_i(k)**2)
  enddo
  !tmp=ti0*kev/(mass_i*vn_i**2)

  !eq_particle_number=ni0*sqrt((mass_i/(twopi*ti0*kev))**3)*ln**3*vn_i**3/normalizing_factor !explicitly cancel the exp(-v^2/vt^2) dependence, valid only for the case that both physical equilibrium distribution and marker distribution are Maxwellian.
  do k=1,nmarker_i
     if( active_i(k).eqv..false.) then
        w_i_star(k)=0._p_
     else
        tmp=ti_func(radcor_i(k))*kev/(mass_i*vn_i**2)
        eq_particle_number=ps_vol_i(k)*fi0(radcor_i(k),v_i(k)) !general case
        tmp2=kappa_ni_func(radcor_i(k))+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti_func(radcor_i(k))

        bz_val=bz(r_i(k),z_i(k))
        br_val=br(r_i(k),z_i(k))
        bphi_val=bphi(r_i(k),z_i(k))
        bval=sqrt(bz_val*bz_val+br_val*br_val+bphi_val*bphi_val)

        call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
             & grad_psi_r,theta_i(k),radcor_i(k),grad_psi_r0)  !in magnetic coordinates
        call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
             & grad_psi_z, theta_i(k),radcor_i(k),grad_psi_z0)  !in magnetic coordinates
        call field_perturbation_on_marker2(radcor_i(k),theta_i(k),alpha_i(k), active_i(k),&
             & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val)

        wiprime=eq_particle_number*twopi/tmp*(ef_cyl_r_val*vr_i(k)+ef_cyl_z_val*vz_i(k)+ef_cyl_phi_val*vphi_i(k)) &
             &+eq_particle_number*tmp2*minor_r_prime(radcor_i(k))/bval**2*&
             (grad_psi_z0*ef_cyl_r_val*bphi_val-grad_psi_r0*ef_cyl_z_val*bphi_val&
             & -grad_psi_z0*ef_cyl_phi_val*br_val+grad_psi_r0*ef_cyl_phi_val*bz_val)

        w_i_star(k)=w_i(k)+wiprime*dtao
     endif
  enddo
end subroutine push_ion_weight_using_electric_cylindrical_components


subroutine ion_velocity_and_metric_in_mc(nmarker_i, r_i , z_i,  &
     & radcor_i , theta_i , active_i,  &
     & vr_i , vphi_i , vz_i,  &
     & grad_psi_i , grad_alpha_i , grad_psi_dot_grad_alpha_i , bval_i,  & 
     & v_i, vpar_i , vx_i , vy_i ) !evaluate variables appearing on the rhs of the ion weight evolution equation, for linear case, these variables are independent of perturbations, and thus can be calculated once and stored, to be used multiple times in the iterations involved in the implicit delta-f method
  use constants,only: p_
  use constants,only: one
  use magnetic_coordinates,only: mpol,nflux,radcor_1d_array,theta_1d_array,r_mc,z_mc, &
  & grad_psi,grad_alpha_r,grad_alpha_z,grad_alpha,grad_psi_dot_grad_alpha, &
  & grad_psi_r,grad_psi_z
use magnetic_field, only : br,bphi,bz
  use interpolate_module
  
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: r_i (nmarker_i), z_i (nmarker_i)
  real(p_),intent(in):: radcor_i (nmarker_i), theta_i (nmarker_i)
  logical,intent(in)::  active_i (nmarker_i)
  real(p_),intent(in):: vr_i (nmarker_i), vphi_i (nmarker_i), vz_i (nmarker_i)
  real(p_),intent(out):: grad_psi_i (nmarker_i), grad_alpha_i (nmarker_i), &
       & grad_psi_dot_grad_alpha_i (nmarker_i), bval_i (nmarker_i)
  real(p_),intent(out):: v_i (nmarker_i), vpar_i (nmarker_i), vx_i (nmarker_i), vy_i (nmarker_i)


  real(p_):: br_val,bphi_val,bz_val
  real(p_):: grad_psi_r_val,grad_psi_z_val,grad_alpha_r_val,grad_alpha_z_val
  integer:: k

  v_i (1:nmarker_i)=sqrt(vr_i (1:nmarker_i)**2+vz_i (1:nmarker_i)**2+vphi_i (1:nmarker_i)**2)

  do k=1, nmarker_i
     if(active_i (k).eqv..false.) cycle
     br_val=br(r_i (k),z_i (k))
     bphi_val=bphi(r_i (k),z_i (k))
     bz_val=bz(r_i (k),z_i (k))
     bval_i (k)=sqrt(br_val**2+bphi_val**2+bz_val**2)

     vpar_i (k)=(br_val*vr_i (k)+bz_val*vz_i (k)+bphi_val*vphi_i (k))/bval_i (k)

     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi_r,theta_i (k),radcor_i (k),grad_psi_r_val)  !in magnetic coordinates

     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi_z,theta_i (k),radcor_i (k),grad_psi_z_val)  !in magnetic coordinates

     vx_i (k)=(vr_i (k)*grad_psi_r_val+vz_i (k)*grad_psi_z_val)

     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
          & grad_alpha_r,theta_i (k),radcor_i (k),grad_alpha_r_val)  !in magnetic coordinates
     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
          & grad_alpha_z,theta_i (k),radcor_i (k),grad_alpha_z_val)  !in magnetic coordinates

     vy_i (k)=vphi_i (k)/r_i (k)+vr_i (k)*grad_alpha_r_val+vz_i (k)*grad_alpha_z_val

     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi,theta_i (k),radcor_i (k),grad_psi_i (k))  !in magnetic coordinates

     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
          & grad_alpha,theta_i (k),radcor_i (k),grad_alpha_i (k))  !in magnetic coordinates

     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,&
          & grad_psi_dot_grad_alpha,theta_i (k),radcor_i (k),grad_psi_dot_grad_alpha_i (k))  !in magnetic coordinates
  enddo

  !for testing
!!$  do k=1,nmarker_i
!!$     if(active_i (k).eqv..true.) then
!!$        cos_theta=grad_psi_dot_grad_alpha_i(k)/(grad_psi_i(k)*grad_alpha_i(k))
!!$        sin_theta=sqrt(1-cos_theta**2) !sign?
!!$        theta1=atan((cos_theta-(vy_i(k)/vx_i(k)*grad_psi_i(k)/grad_alpha_i(k)))/sin_theta)
!!$        vper_i(k)=vx_i(k)/grad_psi_i(k)/cos(theta1)
!!$      write(*,*) sqrt(vpar_i(k)**2+vper_i(k)**2),v_i(k)  !tested, good agreement
!!$    endif
!!$  enddo

end subroutine ion_velocity_and_metric_in_mc

subroutine fk_push_first_step()
  use constants,only:two, one_half
  use fk_module,only: nmarker_i,dtao_i, w_i, w_i_mid
  use fk_module,only: r_i, z_i, phi_i, r_i_old,z_i_old,phi_i_old, r_i_mid, z_i_mid, phi_i_mid
  use fk_module,only: vr_i_old,vz_i_old,vphi_i_old
  use fk_module,only: vr_i_integer, vz_i_integer, vphi_i_integer
  use fk_module,only: vr_i_mid, vz_i_mid, vphi_i_mid
  use fk_module,only: touch_bdry_i,active_i,radcor_i, theta_i, alpha_i
  use fk_module,only: touch_bdry_i_mid,active_i_mid,radcor_i_mid, theta_i_mid, alpha_i_mid
  use fk_module,only: grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i, v_i,vpar_i,vx_i,vy_i
  use fk_particle_coordinates_transform_module,only:compute_particle_magnetic_coordinates
  use sort_ions, only:sort_ions_according_to_poloidal_location
  implicit none
  r_i_old=r_i; z_i_old=z_i; phi_i_old=phi_i
  vr_i_old=vr_i_mid; vz_i_old=vz_i_mid;vphi_i_old=vphi_i_mid
  call push_ion_orbit_first_step(dtao_i)
  vr_i_integer=(vr_i_old+vr_i_mid)/two; vz_i_integer=(vz_i_old+vz_i_mid)/two; vphi_i_integer=(vphi_i_old+vphi_i_mid)/two
!!$  call ion_velocity_and_metric_in_mc(nmarker_i, r_i_old, z_i_old,  &
!!$       & radcor_i, theta_i, active_i,  vr_i_integer, vphi_i_integer, vz_i_integer,  &
!!$       & grad_psi_i , grad_alpha_i , grad_psi_dot_grad_alpha_i , bval_i,  & 
!!$       & v_i, vpar_i , vx_i , vy_i) !output is used by ion weight pusher
!!$
!!$  call push_ion_weight(one_half*dtao_i,nmarker_i,active_i,radcor_i,theta_i,alpha_i, &
!!$       & grad_psi_i , grad_alpha_i, grad_psi_dot_grad_alpha_i , bval_i,  & 
!!$       & v_i, vpar_i , vx_i, vy_i, w_i, w_i_mid) 
    call push_ion_weight_using_electric_cylindrical_components(one_half*dtao_i,radcor_i,theta_i,alpha_i,&
             r_i_old,z_i_old,phi_i_old,vr_i_integer,vz_i_integer,vphi_i_integer,active_i,w_i,w_i_mid,nmarker_i)

  r_i_mid=(r_i+r_i_old)/two;  z_i_mid=(z_i+z_i_old)/two; phi_i_mid=(phi_i+phi_i_old)/two
  call compute_particle_magnetic_coordinates(nmarker_i,r_i_mid,phi_i_mid,z_i_mid,&
       & radcor_i_mid,active_i_mid,touch_bdry_i_mid,theta_i_mid,alpha_i_mid) !prepare ion magnetic coordinates so that markers can be sorted and deposited in the magnetic coordinates
  call sort_ions_according_to_poloidal_location(theta_i_mid) !prepare for ion weight pusher and deposition
end subroutine fk_push_first_step


subroutine fk_push_second_step()
  use fk_module,only: dtao_i, nmarker_i, r_i_mid, z_i_mid, phi_i_mid,&
       & radcor_i_mid , theta_i_mid , alpha_i_mid, active_i_mid,  vr_i_mid , vphi_i_mid , vz_i_mid,  &
       & grad_psi_i_mid , grad_alpha_i_mid , grad_psi_dot_grad_alpha_i_mid , bval_i_mid,  & 
       & v_i_mid, vpar_i_mid , vx_i_mid , vy_i_mid
  use fk_module,only: w_i, w_i_star !output
  use fk_module,only: r_i,phi_i,z_i, radcor_i,theta_i,alpha_i,active_i,touch_bdry_i
  use fk_particle_coordinates_transform_module,only:compute_particle_magnetic_coordinates
  use sort_ions, only:sort_ions_according_to_poloidal_location
  implicit none

!!$  call ion_velocity_and_metric_in_mc(nmarker_i, r_i_mid, z_i_mid,  &
!!$       & radcor_i_mid , theta_i_mid , active_i_mid,  vr_i_mid , vphi_i_mid , vz_i_mid,  &
!!$       & grad_psi_i_mid , grad_alpha_i_mid , grad_psi_dot_grad_alpha_i_mid , bval_i_mid,  & 
!!$       & v_i_mid, vpar_i_mid , vx_i_mid , vy_i_mid) !output is used by ion weight pusher
!!$  call push_ion_weight(dtao_i,nmarker_i,active_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid, &
!!$       & grad_psi_i_mid , grad_alpha_i_mid, grad_psi_dot_grad_alpha_i_mid , bval_i_mid,  & 
!!$       & v_i_mid, vpar_i_mid , vx_i_mid , vy_i_mid,w_i,w_i_star) 
   call push_ion_weight_using_electric_cylindrical_components(dtao_i,radcor_i_mid,theta_i_mid,alpha_i_mid,&
             r_i_mid,z_i_mid,phi_i_mid,vr_i_mid,vz_i_mid,vphi_i_mid,active_i_mid,w_i,w_i_star,nmarker_i)

  w_i(1:nmarker_i)=w_i_star(1:nmarker_i)
  call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,&
       & radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !prepare ion magnetic coordinates so that markers can be sorted and deposited in the magnetic coordinates
  call sort_ions_according_to_poloidal_location(theta_i) !prepare for ion weight pusher and deposition

end subroutine fk_push_second_step

subroutine push_ion_orbit_first_step(dtao_i) !wrapper of "push_full_orbit_cylindrical_boris_with_additional_output" subroutine
  use constants,only:p_
  use constants,only: twopi
  use fk_module,only: nmarker_i,touch_bdry_i
  use fk_module,only: r_i,phi_i,z_i,vr_i,vphi_i,vz_i !as input and output
!  use fk_module,only: phi_i_mid
  use fk_module,only: vr_i_mid,vphi_i_mid,vz_i_mid !as output
  use fk_module,only: radcor_i,theta_i,alpha_i,active_i
  use boris
  implicit none
  real(p_),intent(in):: dtao_i
  integer:: k

  !$omp parallel do
  do k=1,nmarker_i
     if(touch_bdry_i(k).eqv..true.) cycle
     !call push_full_orbit_cylindrical_boris(dtao_i,r_i(k),phi_i(k),z_i(k),vr_i(k),vphi_i(k),vz_i(k))
     call push_full_orbit_cylindrical_boris_with_additional_input_output(dtao_i,radcor_i(k),theta_i(k),alpha_i(k), active_i(k), &
          & r_i(k),phi_i(k),z_i(k),vr_i(k),vphi_i(k),vz_i(k),&
          & vr_i_mid(k),vphi_i_mid(k),vz_i_mid(k)) !obtain velocity at t_{n+1/2}, outputting the projections of this velocity onto both the basis vectors at t_{n+1/2} and those at t_{n+1}
!!$     !{r,z,phi}_i_mid is known before entering this subroutine ({r,z,phi}_i_mid is ethier the initial codintion of the second-pusher or the output of it.)
  enddo
  !$omp end parallel do
end subroutine push_ion_orbit_first_step

!!$subroutine push_ion_orbit_second_step(dtao_i) 
!!$  !calculate the velocity at t_{n+1}, giving the projections of this velocity onto (1) the basis vectors at t_{n+1} (to prepare for the deposition process) and (2) the basis vectors at t_{n+3/2} (to prepare input for the next secod-pusher)
!!$  !output the location at t_{n+3/2} (to prepare input for the next secod-pusher)
!!$  !(the spatial location at t_{n+1} is already computed by the first boris-pusher)
!!$  use constants,only:p_
!!$  use constants,only: twopi
!!$  use fk_module,only: nmarker_i,touch_bdry_i_mid
!!$  use fk_module,only: r_i_mid,phi_i_mid,z_i_mid !as input (t_{n+1/2}) and output (t_{n+3/2})
!!$  use fk_module,only: vr_i_integer_mid,vphi_i_integer_mid,vz_i_integer_mid !input (projection of v at t_{n} onto basis vector at t_{n+1/2}) and output (projection of v at t_{n+1} onto basis vector at t_{n+3/2})
!!$  use fk_module,only: vr_i_integer,vphi_i_integer,vz_i_integer !as output, projection of velocity at t_{n+1} onto the basis vectors at t_{n+1}
!!$  use fk_module,only: phi_i !as input, location at t_{n+1}
!!$  use domain_decomposition,only: myid
!!$  use boris
!!$  implicit none
!!$  real(p_),intent(in):: dtao_i
!!$  integer:: k
!!$
!!$  !$omp parallel do
!!$  do k=1,nmarker_i
!!$     if(touch_bdry_i_mid(k).eqv..true.) cycle
!!$     call push_full_orbit_cylindrical_boris_with_additional_input_output(dtao_i,r_i_mid(k),phi_i_mid(k),z_i_mid(k),&
!!$          & vr_i_integer_mid(k),vphi_i_integer_mid(k),vz_i_integer_mid(k),&
!!$          & phi_i(k),vr_i_integer(k),vphi_i_integer(k),vz_i_integer(k))
!!$ enddo
!!$  !$omp end parallel do
!!$!if(myid.eq.0) write(*,*) z_i_mid(20), dtao_i*vz_i_integer(20)
!!$end subroutine push_ion_orbit_second_step

function fi0(radcor,v) result (z) ! v in unit of vn, fi0 in unit 1/(Ln**3*vn**3)
  use constants,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: vn_i,Ln
  use fk_module,only: mass_i
  use density_temperature_profile_mod,only: ni_func, ti_func
  implicit none
  real(p_),intent(in)  :: radcor,v
  real(p_)  :: z
  real(p_)  :: v_si, ni, ti !local variables
  v_si=v*vn_i
  ti=ti_func(radcor)
  ni=ni_func(radcor)
  z=ni*sqrt((mass_i/(twopi*ti*kev))**3)*exp(-mass_i*v_si**2/(two*ti*kev))
  z=z*(vn_i**3*Ln**3)
end function fi0

end module push_ion_weight_module





