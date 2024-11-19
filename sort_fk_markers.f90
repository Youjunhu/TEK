module sort_ions

contains
subroutine sort_ions_according_to_poloidal_location(theta)
  use constants,only:p_
  use constants,only: twopi
  use pputil
  use fk_module
  implicit none
  real(p_),intent(in):: theta(:)
  integer:: ierr,np_old,np_new
  !assign particles to the different processors according to their theta coordinates, using the subroutines provided in pputil_yj.f90

  np_old=nmarker_i
  call init_pmove(theta(:),np_old,twopi,ierr)

  call pmove(ps_vol_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_i_star(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(r_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

 call pmove(r_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phi_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(radcor_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(theta_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(radcor_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(theta_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(vr_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(vr_i_integer_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i_integer_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i_integer_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vr_i_integer(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i_integer(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i_integer(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(vr_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
!!$
  call pmove(v_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vx_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vy_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_alpha_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_dot_grad_alpha_i(:),    np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(bval_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(v_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vx_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vy_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_alpha_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_dot_grad_alpha_i_mid(:),    np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(bval_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove2(active_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(active_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(touch_bdry_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(touch_bdry_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call end_pmove(ierr)
  nmarker_i=np_new

  !     call check_domain_particles(theta,nmarker_i)
end subroutine sort_ions_according_to_poloidal_location

end module sort_ions

subroutine check_domain_particles(theta,nmarker_i) !pass the test, comfirming domain decomposition is consistent with particles grouping
  use constants,only:p_
  use domain_decomposition,only:theta_start,dtheta2
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: theta(nmarker_i)
integer:: k

do k=1,nmarker_i
if(theta(k)<theta_start .or. theta(k)>theta_start+dtheta2) write(*,*) 'warningg*** particle not in domain'
enddo

end subroutine check_domain_particles
