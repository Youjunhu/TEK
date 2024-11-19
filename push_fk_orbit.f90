module fk_particle_coordinates_transform_module
contains
subroutine compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !get the radial coordinates of markers, relocate markers if they are outside of the computational box, then calculate the (theta,alpha) coordinates
  use constants,only:p_
  use constants,only:pi,twopi
  use poloidal_flux_2d,only:xarray,zarray,nx,nz
  use magnetic_coordinates,only:radcor_low2,radcor_upp2,radcor_low0,radcor_upp0
  use magnetic_coordinates,only:mpol,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc
  use magnetic_field,only:radcor_as_func_of_pfn
  use magnetic_field, only : pfn_func
  use interpolate_module
    use math,only: shift_to_specified_toroidal_range
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(inout):: r_i(nmarker_i),phi_i(nmarker_i),z_i(nmarker_i)
  real(p_),intent(out):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  logical,intent(out):: active_i(nmarker_i),touch_bdry_i(nmarker_i)
  integer:: k,flag,next_seed
  real(p_):: theta,rannum1,rannum2
  logical:: outside_box
  real(p_):: tor_shift

  do k=1,nmarker_i
     outside_box=r_i(k).ge.xarray(nx) .or.r_i(k).le.xarray(1) .or. z_i(k).ge.zarray(nz) .or.  z_i(k).le.zarray(1) 
     if(outside_box.eqv..true.) then
        touch_bdry_i(k)=.true. !marker is lost forever, will never be re-introduced
        active_i(k)=.false.
     else
        radcor_i(k)=radcor_as_func_of_pfn(pfn_func(r_i(k),z_i(k))) !calculate the radial coordinate
        if(radcor_i(k).ge.radcor_upp0 .or. radcor_i(k).le.radcor_low0) then 
           touch_bdry_i(k)=.true. !marker is lost forever, will never be re-introduced
           active_i(k)=.false.
        else if(radcor_i(k).lt.radcor_low2 .or. radcor_i(k).gt.radcor_upp2) then
           z_i(k)=-z_i(k) !relocate the marker by up-down reversing
           radcor_i(k)=radcor_as_func_of_pfn(pfn_func(r_i(k),z_i(k))) !re-calculate the radial coordinate, which is different from the original value if the flux-surface is not up-down symmetric
           active_i(k)=.false.
           touch_bdry_i(k)=.false.
        else
           touch_bdry_i(k)=.false.
           active_i(k)=.true.
        endif
     endif
  enddo

!calculate (theta,alpha) coordinates
 do k=1,nmarker_i
     !if(active_i(k).eqv..false.) cycle !wrong, inactive markers' theta must be computed so that they can be sorted by the sorting subroutine
     if(touch_bdry_i(k).eqv..true.) cycle
     call interpolate_from_cylindrical_to_magnetic_coordinates1(r_i(k),z_i(k),theta_i(k))
     call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc,&
          & theta_i(k),radcor_i(k),tor_shift) !interpolating in magnetic coordinates to get tor_shift
     alpha_i(k)=phi_i(k)-tor_shift !generalized toroidal angle
     call shift_to_specified_toroidal_range(alpha_i(k))
     call shift_to_specified_toroidal_range(phi_i(k)) !phi_i is needed in deposition, so should be shifted to the desired range
  enddo

end subroutine compute_particle_magnetic_coordinates


subroutine count_lost_markers_fk() 
  use fk_module
  use domain_decomposition,only: myid
  implicit none
  integer:: k, nlost
  nlost=0
  do k=1,nmarker_i
     if(touch_bdry_i(k).eqv..true.) then
        nlost=nlost+1
     endif
  enddo
write(*,*) 'number of lost fk markers=', nlost, 'myid=',myid, ',total number=', nmarker_i
end subroutine count_lost_markers_fk


subroutine clean_up_lost_markers_fk() !for fully kinetic markers, drop lost markers so that we we have smaller particle arrays
  use fk_module,only:nmarker_i,w_i, ps_vol_i, active_i, touch_bdry_i, radcor_i, theta_i, alpha_i,&
       & r_i,phi_i,z_i, vr_i,vphi_i,vz_i, vr_i_mid, vz_i_mid, vphi_i_mid
  implicit none
  integer:: k, kclean

  kclean=0
  do k=1,nmarker_i
     if(touch_bdry_i(k).eqv..false.) then

        kclean=kclean+1

        r_i(kclean)=r_i(k)
        phi_i(kclean)=phi_i(k)
        z_i(kclean) =z_i (k)
        vr_i(kclean)=vr_i(k)
        vphi_i(kclean)=vphi_i(k)
        vz_i(kclean)=vz_i(k)
        w_i(kclean)=w_i(k)
        ps_vol_i(kclean)=ps_vol_i(k)
        active_i(kclean)=active_i(k)
        touch_bdry_i(kclean)=touch_bdry_i(k)
        radcor_i(kclean)=radcor_i(k)
        theta_i(kclean)=theta_i(k)
        alpha_i(kclean)=alpha_i(k)
        vr_i_mid(kclean)=vr_i_mid(k)
        vz_i_mid(kclean)= vz_i_mid(k)
        vphi_i_mid(kclean)= vphi_i_mid(k)
     endif
  enddo
  nmarker_i=kclean !update the number of markers
end subroutine clean_up_lost_markers_fk


end module fk_particle_coordinates_transform_module
