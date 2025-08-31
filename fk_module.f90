module fk_module !fully-kinetic for ions
  use constants,only:p_, zero
  implicit none
  save
  real(p_) :: mass_i,charge_i,dtao_fk
  real(p_) :: ni0,ti0,kappa_ni,kappa_ti !ti0 is a typical value of ion temperature in the compuational region, in flux tube model, ti0 is the temperature at the reference magnetic surface
  real(p_) :: vt_i,vmin_i,vmax_i
  real(p_) ::normalizing_factor
  real(p_) :: omegan_fk, tn_fk, vn_fk
  
  integer :: ntouch_bdry_i=0, total_ntouch_bdry_i !numbe of markers that touch the boundary in each process and all processes, respectiveyl
  integer :: total_nmarker_i  ! total number of ion markers (including the particles in all the processors).
  integer :: nmarker_i_per_cell  ! total number of ion markers (including the particles in all the processors).
  integer :: nmarker_i !particle number in a single processor, its value will be differnt for differnt processors and at differnt time

  real(p_),allocatable :: ps_vol_i(:) !determined by the initial loading, is constant for each marker over the time-evoultion
  real(p_),allocatable :: w_i(:)  !weight of ion markers
  real(p_),allocatable :: w_i_mid(:)  !weight of ion markers at t_{n+1/2}
  real(p_),allocatable :: w_i_star(:)  !weight of ion markers

  real(p_),allocatable :: r_i(:),z_i(:),phi_i(:) !Cylindrical coordinates at time t_{n},  unit Ln, rad
  real(p_),allocatable :: r_i_old(:),z_i_old(:),phi_i_old(:) !Cylindrical coordinates at time t_{n},  unit Ln, rad, temporary working arrays for computing averaing
  real(p_),allocatable :: r_i_mid(:),z_i_mid(:),phi_i_mid(:) !Cylindrical coordinates at time t_{n+1/2}, unit Ln, rad, 
  real(p_),allocatable :: radcor_i(:),theta_i(:), alpha_i(:),tor_shift_i(:) !magnetic coordinates at integer-time-step, alpha is the generalized toroidal angle
  real(p_),allocatable :: radcor_i_mid(:),theta_i_mid(:),alpha_i_mid(:) !magnetic coordinates at half time-step

  real(p_),allocatable :: vr_i(:),vz_i(:),vphi_i(:) ! projection of velocity at t_{n-1/2} to the basis cylindrical vector at integer-time-step t_{n}, unit vn_fk=ln/tn_fk, 
  real(p_),allocatable :: vr_i_old(:),vz_i_old(:),vphi_i_old(:) !temporary working arrays for computing averaing
  real(p_),allocatable :: vr_i_integer_mid(:),vz_i_integer_mid(:),vphi_i_integer_mid(:) !projection of velocity at t_{n} to the basis cylindrical vector at t_{n+1/2}, initial condition for the second boris pusher, unit vn_fk=ln/tn_fk, 

  real(p_),allocatable :: vr_i_mid(:),vz_i_mid(:),vphi_i_mid(:) !projection of velocity at t_{n+1/2} to the local basis vectors at t_{n+1/2}
  real(p_),allocatable :: vr_i_integer(:),vz_i_integer(:),vphi_i_integer(:) !projection of velocity at t_{n} to the local basis vectors at t_{n}

  real(p_),allocatable :: v_i(:),vpar_i(:),vx_i(:),vy_i(:) !vx is defined by vx=v_dot_grad_x, vy is defined by vy=v_dot_grad_y,note that grad_x and grad_y are not perpendicular to each other
  real(p_),allocatable :: grad_psi_i(:),grad_alpha_i(:), grad_psi_dot_grad_alpha_i(:),bval_i(:)
  real(p_),allocatable :: v_i_mid(:),vpar_i_mid(:),vx_i_mid(:),vy_i_mid(:)
  real(p_),allocatable :: grad_psi_i_mid(:),grad_alpha_i_mid(:), grad_psi_dot_grad_alpha_i_mid(:),bval_i_mid(:)

  logical,allocatable :: touch_bdry_i(:),active_i(:) !indicates whether the orbit of a marker touches the boundary
  logical,allocatable :: touch_bdry_i_mid(:),active_i_mid(:)
  real(p_),dimension(:,:),allocatable :: my_den_i_left, my_den_i_right !fk ion density

  integer  :: ion_spatial_loading_scheme  !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
  integer  :: ion_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
  integer :: fk_nonlinear
contains
  subroutine initialize_fk()
    use domain_decomposition,only: numprocs, myid
    use magnetic_coordinates, only : nrad, mpol2,mtor
    namelist/fk_nmlt/mass_i,charge_i,ni0,ti0,kappa_ni,kappa_ti,nmarker_i_per_cell, &
         & ion_spatial_loading_scheme, ion_velocity_loading_scheme, fk_nonlinear
    integer:: fixed_large_size, u

     open(newunit=u,file='input.nmlt')
     read(u,fk_nmlt)
     close(u)
     if(myid==0)  write(*,fk_nmlt)

    
    total_nmarker_i=nmarker_i_per_cell*nrad*mpol2*mtor
    if(myid.eq.0) write(*,*) 'total number of ions=',     total_nmarker_i
    nmarker_i=total_nmarker_i/numprocs !nmarker_i initially store the number of markers initially loaded per processor (i.e.total_nmarker_i/numprocs), latter actual number of markers per proc will be assigned to nmarker_i, the value of which will be differnt for differnt processors and at differnt time
    fixed_large_size=(total_nmarker_i/numprocs)*3/2 !the number of particle per proc after re-arranging the particles between the processors may exceed the number of original loaded particles per proc (i.e., total_nmarker_i/numprocs), increasing the array length by a factor of 3/2 is needed to make sure that the array is big enough to contain all the particles that belong to the domain for which the processor is responsible.
    !  write(*,*), 'nmarker_i, fixed_large_size=',nmarker_i, fixed_large_size

    allocate(radcor_i(fixed_large_size)) 
    allocate(theta_i(fixed_large_size))
    allocate(alpha_i(fixed_large_size))  
    allocate(tor_shift_i(fixed_large_size))

    allocate(radcor_i_mid(fixed_large_size)) 
    allocate(theta_i_mid(fixed_large_size))
    allocate(alpha_i_mid(fixed_large_size))

    allocate(v_i(fixed_large_size))
    allocate(vr_i(fixed_large_size))
    allocate(vz_i(fixed_large_size))
    allocate(vphi_i(fixed_large_size))

    allocate(r_i(fixed_large_size))
    allocate(z_i(fixed_large_size))
    allocate(phi_i(fixed_large_size))
    allocate(r_i_mid(fixed_large_size)) 
    allocate(z_i_mid(fixed_large_size)) 
    allocate(phi_i_mid(fixed_large_size))

    allocate(w_i(fixed_large_size)) 
    allocate(w_i_mid(fixed_large_size)) 
    allocate(w_i_star(fixed_large_size)) 

    allocate(ps_vol_i(fixed_large_size))
    allocate(active_i(fixed_large_size)) !whether particles are within computational boundary
    allocate(active_i_mid(fixed_large_size)) !whether particles are within computational boundary
    allocate(touch_bdry_i(fixed_large_size)) !whether particles are within computational boundary
    allocate(touch_bdry_i_mid(fixed_large_size)) !whether particles are within computational boundary

    allocate(vr_i_integer(fixed_large_size)) 
    allocate(vz_i_integer(fixed_large_size)) 
    allocate(vphi_i_integer(fixed_large_size)) 

    allocate(vr_i_integer_mid(fixed_large_size)) 
    allocate(vz_i_integer_mid(fixed_large_size)) 
    allocate(vphi_i_integer_mid(fixed_large_size)) 

    allocate(vr_i_mid(fixed_large_size)) 
    allocate(vz_i_mid(fixed_large_size)) 
    allocate(vphi_i_mid(fixed_large_size)) 

    allocate(vpar_i(fixed_large_size)) !velocity components in magnetic coordinates
    allocate(vx_i(fixed_large_size)) 
    allocate(vy_i(fixed_large_size)) 
    allocate(grad_psi_i(fixed_large_size)) 
    allocate(grad_alpha_i(fixed_large_size)) 
    allocate(grad_psi_dot_grad_alpha_i(fixed_large_size)) 
    allocate(bval_i(fixed_large_size)) 

    allocate(r_i_old(fixed_large_size)) 
    allocate(z_i_old(fixed_large_size)) 
    allocate(phi_i_old(fixed_large_size)) 

    allocate(vr_i_old(fixed_large_size)) 
    allocate(vz_i_old(fixed_large_size)) 
    allocate(vphi_i_old(fixed_large_size)) 

    allocate(v_i_mid(fixed_large_size))
    allocate(vpar_i_mid(fixed_large_size))
    allocate(vx_i_mid(fixed_large_size))
    allocate(vy_i_mid(fixed_large_size))

    allocate(grad_psi_i_mid(fixed_large_size))
    allocate(grad_alpha_i_mid(fixed_large_size))
    allocate( grad_psi_dot_grad_alpha_i_mid(fixed_large_size))
    allocate(bval_i_mid(fixed_large_size))

    allocate(my_den_i_left(mtor,nrad), source=zero)
    allocate(my_den_i_right(mtor,nrad), source=zero)
  
  end subroutine initialize_fk

end module fk_module



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


subroutine load_fk()
  !spatial location in magnetic coordinates (psi,theta,phi) and then transform to cylindrical coordinates
  use constants,only:p_, one,twopi,pi,two,kev,fourpi
 
  use magnetic_coordinates,only: radcor_low=>xlow,radcor_upp=>xupp,vol, &
       &  mpol,nrad,xgrid,zgrid,jacobian,toroidal_range,jacobian
  use misc, only: magnetic_coordinates_to_cylindrical_coordinates
  use fk_module,only: total_nmarker_i,mass_i,ti0, vn_fk
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
  use fk_module, only: ion_spatial_loading_scheme, ion_velocity_loading_scheme
  use pputil !containing subroutines that sort particles into differnt processors
  use domain_decomposition,only: numprocs,myid
  use math, only :  random_yj, sub_random_yj
  use misc, only: abs_jacobian_func
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
  real(p_) :: abs_jacobian_max
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

abs_jacobian_max=maxval(abs(jacobian(:,1:nrad)))
  !  radcor_min=minval(xgrid)
  !  radcor_max=maxval(xgrid)


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
  vmin=-3._p_*vt/vn_fk !normalized by vn_fk
  vmax=+3._p_*vt/vn_fk !normalized by vn_fk
  maxwellian_max=maxwellian_func_ion(0._p_)

!if(myid.eq.0) write(*,*)  'vt=',vt, 'vmin=',vmin,'vmax=',vmax, 'vmin*vn_fk=',vmin*vn_fk, 'vt/vn_fk=',vt/vn_fk

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
           pos1=maxwellian_func_ion(v_val*vn_fk)
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

  !  v_i=v_i/vn_fk !normalized by vn_fk
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
        call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor)
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor)
     enddo
  elseif(ion_spatial_loading_scheme.eq.1 .and. ion_velocity_loading_scheme.eq.2) then
     !normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*twopi*pi*vt/vn_fk*sqrt(pi)/two)
     normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*(sqrt(pi)*vt/vn_fk)**3)
     do i=1,nmarker_i
        call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk))
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk))
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.1) then
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*(vmax-vmin))
     normalizing_factor=total_nmarker_i/(vol*(vmax-vmin)**3)
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor)
        ps_vol_i(i)=one/(normalizing_factor)
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.2) then
     !  normalizing_factor=total_nmarker_i/(vol*fourpi*vt/vn_fk*sqrt(pi)/two) !wrong
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*vt/vn_fk*sqrt(pi)/two) !wrong again
     normalizing_factor=total_nmarker_i/(vol*(sqrt(pi)*vt/vn_fk)**3) !corrected
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk)) !wrong
        ps_vol_i(i)=one/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk))
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

