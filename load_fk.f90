subroutine load_fk()
  use constants,only:p_, one,twopi,pi,two,kev,fourpi
  use control_parameters,only: ion_spatial_loading_scheme  !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
  use control_parameters,only: ion_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
  use normalizing,only: vn_i
  use magnetic_coordinates,only: radcor_low=>radcor_low2,radcor_upp=>radcor_upp2,vol, j_low=>j_low2, j_upp=>j_upp2 !as input
  use magnetic_coordinates,only: mpol,nflux,radcor_1d_array,theta_1d_array,jacobian,toroidal_range,jacobian
  use fk_module,only: total_nmarker_i,mass_i,ti0 !as input
  use fk_module,only: nmarker_i,radcor_i,theta_i,active_i ! as output
  use fk_module,only:active_i_mid,touch_bdry_i_mid,touch_bdry_i
  use fk_module,only: alpha_i,tor_shift_i !only allocate these array, their values are not assigned in the present subroutine
  use fk_module,only: r_i,z_i,phi_i,v_i,vr_i,vz_i,vphi_i !as output, loading using magnetic coordinates, the Boris pusher works in cylindrical coordinates, therefore, we need to transform from mag. cor. to cylin. cor.
  use fk_module,only: r_i_old,z_i_old,phi_i_old
  use fk_module,only: ps_vol_i,normalizing_factor !as output
  use fk_module,only: w_i,w_i_mid,w_i_star !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vpar_i, vx_i,vy_i  !only allocate the array, the value is not set in this subroutine
  use fk_module,only: grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i !only allocate the array, the value is not set in this subroutine
  use fk_module,only: r_i_mid,z_i_mid,phi_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vr_i_old,vz_i_old,vphi_i_old
  use fk_module,only: vr_i_mid,vz_i_mid,vphi_i_mid !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vr_i_integer,vz_i_integer,vphi_i_integer !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vr_i_integer_mid,vz_i_integer_mid,vphi_i_integer_mid
  use fk_module,only: vt_i,vmin_i,vmax_i !as output, vt_i in SI unit
  use pputil !containing subroutines that sort particles into differnt processors
  use domain_decomposition,only: numprocs,myid
  use math, only :  random_yj, sub_random_yj
  use interpolate_module
  implicit none

  integer:: iseed,next_seed
  integer,parameter:: max_try=10000
  real(p_):: radcor_val,theta_val,rannum1,rannum2,rannum3,tmp
  !  real(p_):: random_yj
  integer:: i,ierr,j,file_unit
  !  integer:: status(MPI_STATUS_SIZE)
  character(5):: filename
  real(p_):: pos1,pos2,jacobian_val
  real(p_):: abs_jacobian_func,abs_jacobian_max
  integer:: np_old,np_new
  real(p_):: vt,vmin,vmax,v_val,maxwellian_func_ion
  !real(p_),allocatable:: theta_v(:),phi_v(:)
  !real(p_),allocatable::vx(:),vy(:),vz(:),tmp_array(:)
  real(p_)::vx(nmarker_i),vy(nmarker_i),vz(nmarker_i),tmp_array(3*nmarker_i)
  real(p_)::maxwellian_min,maxwellian_max

  !  allocate(pitch_angle_i(fixed_large_size))
  !  allocate(gyro_angle_i(fixed_large_size))

  !  allocate(theta_v(nmarker_i)) !local array
  !  allocate(phi_v(nmarker_i)) !local array
!!$  allocate(vx(nmarker_i)) !local array
!!$  allocate(vy(nmarker_i)) !local array
!!$  allocate(vz(nmarker_i)) !local array
!!$  allocate(tmp_array(3*nmarker_i))

abs_jacobian_max=maxval(abs(jacobian(:,j_low:j_upp)))
  !  radcor_min=minval(radcor_1d_array)
  !  radcor_max=maxval(radcor_1d_array)


  ! ---random generator, when use MPI_send to generate iseed for other processes, it is actual a sequence generator,instead of parallel generator
!!$  if ( myid .eq. 0 ) then ! master generates random numbers first, others wait in line
!!$     iseed = 0
!!$  else 
!!$     call MPI_Recv(iseed, 1, MPI_INT, myid-1, 1, MPI_COMM_WORLD, status,ierr) !other processes wait to receive the iseed
!!$  endif

  iseed=-(1777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
  !  write(*,*) 'myid=',myid, 'iseed=',iseed

  ! now generate the random numbers
  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 

  if(ion_spatial_loading_scheme.eq.1) then
     do i=1,nmarker_i
        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
        call sub_random_yj(0,next_seed,rannum2) 
        radcor_val=radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
        !theta_val=rannum2*twopi  !theta in [0:2*pi]
        theta_val=(rannum2-0.5_p_)*twopi !theta in [-pi:+pi]
        radcor_i(i)=radcor_val
        theta_i(i)=theta_val
     enddo
  elseif (ion_spatial_loading_scheme.eq.2) then

     do i=1,nmarker_i     
        do j=1,max_try !rejection method to generate nonuniform random numbers
           call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
           call sub_random_yj(0,next_seed,rannum2) !use last random number as iseed
           !call sub_random_yj(0,next_seed,rannum3) !use last random number as iseed
           radcor_val=radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
           theta_val=-pi+rannum2*twopi
           pos1=abs_jacobian_func(theta_val,radcor_val)
           !       write(*,*) 'abs_jacobian_func(theta_val,radcor_val)=',pos1
           call sub_random_yj(0,next_seed,pos2) !use last random number as iseed   
           pos2=pos2*abs_jacobian_max !scaled to the range [0: abs_jacobian_max]
           if(pos1<pos2) then
              cycle
           else
              radcor_i(i)=radcor_val
              theta_i(i)=theta_val
              !phi_i(i)=toroidal_range*rannum3
              exit
           endif
        enddo
        !     if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution"
        !     write(*,*) 'j=',j
     enddo
  else
     stop 'please specify a loading scheme for the spatial distribution ion markers'
  endif

  do i=1,nmarker_i
     call magnetic_coordinates_to_cylindrical_coordinates(theta_i(i),radcor_i(i),r_i(i),z_i(i)) !to get the corresponding (R,Z) coordinates
  enddo

  do i=1,nmarker_i !setting toroidal coordinate of particles
     call sub_random_yj(0,next_seed,rannum3) !use last random number as iseed
     phi_i(i)=toroidal_range*rannum3
  enddo

  !setting velocity
  vt=sqrt(two*ti0*kev/mass_i)
  vmin=-3._p_*vt/vn_i !normalized by vn_i
  vmax=+3._p_*vt/vn_i !normalized by vn_i
  maxwellian_max=maxwellian_func_ion(0._p_)

!if(myid.eq.0) write(*,*)  'vt=',vt, 'vmin=',vmin,'vmax=',vmax, 'vmin*vn_i=',vmin*vn_i, 'vt/vn_i=',vt/vn_i

  if(ion_velocity_loading_scheme.eq.1) then !using uniform loading in v, instead of Gaussian
     do i=1,3*nmarker_i        
        call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
        v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
        tmp_array(i)=v_val
     enddo
  elseif (ion_velocity_loading_scheme.eq.2) then
     do i=1,3*nmarker_i        
        do j=1,max_try !rejection method to generate nonuniform random numbers
           call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
           v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
           pos1=maxwellian_func_ion(v_val*vn_i)
           call sub_random_yj(0,next_seed,pos2) !0 means using last random number as iseed   
           pos2=pos2*maxwellian_max !scaled to [0,maxwellian_max]
           if(pos1<pos2) then
              cycle
           else
              tmp_array(i)=v_val
              exit
           endif
        enddo
!        if(myid.eq.0) write(*,*) 'j=',j
     enddo
  else
     stop 'please specify a loading scheme for the velocity distribution ion markers'
  endif

  do i=1,nmarker_i
     vx(i)=tmp_array(i)
     vy(i)=tmp_array(i+nmarker_i)
     vz(i)=tmp_array(i+2*nmarker_i)
  v_i(i)=sqrt(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
  enddo

!if(myid.eq.0) call calculate_possibility_density(vz,nmarker_i,100,vmin,vmax)

  do i=1,nmarker_i
     vz_i(i)=vz(i)
     vr_i(i)=vx(i)*cos(phi_i(i))+vy(i)*sin(phi_i(i))
     vphi_i(i)=-vx(i)*sin(phi_i(i))+vy(i)*cos(phi_i(i))
  enddo

!  v_i=sqrt(vr_i*vr_i+vz_i*vz_i+vphi_i*vphi_i)
  !if(myid.eq.3) call calculate_possibility_density(v_i,nmarker_i,100,vmin,vmax)

  !  v_i=v_i/vn_i !normalized by vn_i
!!$  do i=1,nmarker_i !setting direction of velocity
!!$     call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
!!$     theta_v(i)=pi*rannum1
!!$     call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed
!!$     phi_v(i)=twopi*rannum1
!!$  enddo
  !transform to components in cylindrical coordinates
!!$  do i=1,nmarker_i !velocity components in a constant Cartesian coordinate system
!!$     vz_i(i)=v_i(i)*cos(theta_v(i))
!!$     vx(i)=v_i(i)*sin(theta_v(i))*cos(phi_v(i))
!!$     vy(i)=v_i(i)*sin(theta_v(i))*sin(phi_v(i))
!!$  enddo

!!$  do i=1,nmarker_i !projected onto the basis vectors of cylindrical coordinates
!!$     vr_i(i)=vx(i)*cos(phi_i(i))+vy(i)*sin(phi_i(i))
!!$     vphi_i(i)=vy(i)*cos(phi_i(i))-vx(i)*sin(phi_i(i))
!!$  enddo

  if(ion_spatial_loading_scheme.eq.1 .and. ion_velocity_loading_scheme.eq.1) then
     !normalizing_factor=total_nmarker_i/(twopi*toroidal_range*(radcor_upp-radcor_low)*twopi*pi*(vmax-vmin))
     normalizing_factor=total_nmarker_i/(twopi*toroidal_range*(radcor_upp-radcor_low)*(vmax-vmin)**3)
     do i=1,nmarker_i
        call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor)
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor)
     enddo
  elseif(ion_spatial_loading_scheme.eq.1 .and. ion_velocity_loading_scheme.eq.2) then
     !normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*twopi*pi*vt/vn_i*sqrt(pi)/two)
     normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*(sqrt(pi)*vt/vn_i)**3)
     do i=1,nmarker_i
        call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i))
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i))
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.1) then
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*(vmax-vmin))
     normalizing_factor=total_nmarker_i/(vol*(vmax-vmin)**3)
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor)
        ps_vol_i(i)=one/(normalizing_factor)
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.2) then
     !  normalizing_factor=total_nmarker_i/(vol*fourpi*vt/vn_i*sqrt(pi)/two) !wrong
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*vt/vn_i*sqrt(pi)/two) !wrong again
     normalizing_factor=total_nmarker_i/(vol*(sqrt(pi)*vt/vn_i)**3) !corrected
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i)) !wrong
        ps_vol_i(i)=one/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_i))
     enddo
  endif


!!$   iseed=next_seed
!!$  if (myid .ne. numprocs-1) then
!!$     call MPI_Send(iseed, 1, MPI_INT, myid+1, 1, MPI_COMM_WORLD,ierr)  !send the iseed to next process
!!$  endif

  !assign the loaded particles to the corresponding processors, using the subroutines provided in pputil_yj.f90
  np_old=nmarker_i
  call init_pmove(theta_i(:),np_old,twopi,ierr)
  call pmove(theta_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(radcor_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(r_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(v_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vr_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(ps_vol_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call end_pmove(ierr)
  nmarker_i=np_new

  !  call some_test3(myid,numprocs)

!!$  write(filename,'(a1,i4.4)') 'i',myid
!!$  open(newunit=file_unit, file=filename)
!!$  do i=1,nmarker_i
!!$     write(file_unit,'(4(1pe14.5),i6,i3)')  radcor_i(i) ,theta_i(i),r_i(i),z_i(i),i,myid
!!$     !     write(file_unit,*)  vr_i(i),vphi_i(i),vz_i(i)
!!$  enddo
!!$  close(file_unit)

  active_i=.true. ! initially, all markers are active, i.e., within the computational region
  touch_bdry_i=.false.
  active_i_mid=.true. ! initially, all markers are active, i.e., within the computational region
  touch_bdry_i_mid=.false.

  vmin_i=vmin
  vmax_i=vmax
  vt_i=vt
end subroutine load_fk


function maxwellian_func_ion(v) result(z)
use constants,only:p_
  use constants,only:two,kev
  use fk_module,only: mass_i,ti0 !as input
implicit none
real(p_):: v,z

z=exp(-mass_i*v*v/(two*ti0*kev)) !the normalizing factor (mi/(twopi*Ti*kev))^(3/2) is not included


end function maxwellian_func_ion






