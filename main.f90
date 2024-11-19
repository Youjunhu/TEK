!code name: TEK (Tokamak Electromagnetic Kinetic simulation)
!Particle in cell simulation of low-frequency electromagnetic turbulence in tokamak plasmas
!using gyrokinetic and/or fully kinetic model for ions, and gyrokinetic for electrons
!Main devevloper: Youjun Hu, Email: yjhu@ipp.cas.cn, youjunhu@gmail.com

program main
  use constants,only:p_, twopi,pi,kev,two,mu0, one_half, zero
  use control_parameters,only:kstart,kend,dtao_omega_i_axis,niter,iplot_mode_structure,&
       &  fk_ions_switch, gk_species_switch, diagnosis
  use normalizing,only: bn,ln,tn_i, vn_i, beta_ni,vu, tn_main, dtao_main, omegan_i
  use fk_module,only: mass_i,charge_i,dtao_i, nmarker_i,vt_i,ti0,ni0,ti0,kappa_ti, &
       & nmarker_i_per_cell, r_i,z_i,phi_i, radcor_i,theta_i,alpha_i,tor_shift_i,ps_vol_i,&
       & touch_bdry_i,active_i, touch_bdry_i_mid,active_i_mid, &
       & radcor_i_mid,theta_i_mid,alpha_i_mid, &
       & vr_i,vz_i,vphi_i, ntouch_bdry_i, total_ntouch_bdry_i, &
       & vr_i_integer,vz_i_integer,vphi_i_integer,&
       & vr_i_mid,vphi_i_mid,vz_i_mid, &
       & r_i_mid,z_i_mid,phi_i_mid,vr_i_integer_mid,vz_i_integer_mid,vphi_i_integer_mid, &
       & grad_psi_i_mid,grad_alpha_i_mid,grad_psi_dot_grad_alpha_i_mid,bval_i_mid, &
       & v_i_mid,vpar_i_mid,vx_i_mid,vy_i_mid, &
       & w_i,w_i_star,w_i_mid, allocate_fk_particle_array
  use gk_module,only:nsm, mass_e, charge_e,charge_sign_e,te0,nm_gk, &
       & dtao_e,ne0,kappa_te, rho_e, omegan_e, omega_e_axis,tn_e, vn_e, beta_ne, w_unit, &
       & touch_bdry_e, &
       & radcor_e,theta_e,alpha_e,vpar_e,mu_e, v_e, ptcl_num0_e, w_e, w_e_mid, &
       & radcor_e_mid,theta_e_mid,alpha_e_mid,vpar_e_mid, &
       & equ_radial_drift, equ_theta_drift, equ_alpha_drift,mirror_force, equ_theta_drift0,&
       & perturbed_radial_drift, perturbed_theta_drift, perturbed_alpha_drift, &
       & x_ring, y_ring, phix_e, phiy_e, phiz_e, ax_e, ay_e, az_e, ahx_e, ahy_e, ahz_e, ah_e, &
       & allocate_gk_particle_array, sort_gk_markers
  use load_gk_mod, only: load_gk
  use gyro_ring_mod, only : set_gyro_phase, gyro_ring
  use gyro_average_mod, only : gyro_average
  use drift, only : compute_drift
  use gk_trajectory_pusher, only: push_gc_half_step, push_gc_full_step
  use gk_weight_pusher,only: push_gk_weight
  use magnetic_coordinates,only: mpol,nflux,radcor_low2,radcor_upp2, nsegment, dtheta, vol, grad_psi,&
       & mpol2, nflux2,mtor,tor_1d_array,radcor_1d_array,theta_1d_array,tor_shift_mc, toroidal_range, GSpsi_prime
  use magnetic_field, only: qfunc, pfn_func, minor_r_radcor
  use table_in_mc,only: b_mc, prepare_table_in_mc
  use math,only: shift_to_specified_toroidal_range
  use radial_module,only: baxis, r_axis
  use flux_tube_model,only: radcor_fixed, j_fixed
  use perturbation_field_matrix,only:   allocate_field_matrix, potential, apara_s, apara_s_old,&
       & my_den_e_left, my_den_e_right, my_jpar_e_left, my_jpar_e_right
  use domain_decomposition,only: myid,numprocs,tube_comm,grid_comm,ntube,gclr,tclr,GCLR_cut, GCLR_cut_left,&
       & dtheta2,theta_start,my_right,my_left, my_right2, my_left2, multi_eq_cells, ipol_eq
  use pputil, only: ppinit
  use fk_particle_coordinates_transform_module,only: compute_particle_magnetic_coordinates, clean_up_lost_markers_fk, &
       & count_lost_markers_fk
  use mode_structure
  use diagnostic_mod
  use deposit_fk_module,only : deposit_ions
  use deposit_gk_module,only : deposit_gk
  use FFTW3
  use spectrum_diagnostic
  use report_module,only: report, mode_evolution_analysis2, mode_evolution_analysis5
  use poisson, only: prepare_poisson_matrix,solve_poisson
  use ampere, only : prepare_ampere_matrix, solve_ampere, apara_s_evolution, apara_resplit_and_weight_pullback
  use restart_mod
  use push_ion_weight_module
  use initial_half_step_for_boris
  use initialization_mod
  use interpolate_module, only : linear_2d_interpolation
  use communication_connection
  use mpi
  implicit none

  integer ::   ierr, id_writing_evolution
  real(p_) :: t_omega_i_axis
  integer :: k,kt,i_new,iter,i,j
  real(p_) :: omega_i_axis,rho_i,k_binormal1,k_binormal2, bval0
  real(p_) :: radcor_as_func_of_pfn
  real(p_) :: minor_r_min,minor_r_max,minor_r_width,va0
  real(p_) :: t1,t2, tarray(2) !store the cpu clock time
  INTEGER :: count1,count2, count3, count_rate, count_max
  integer :: file_unit_i,file_unit_e,file_unit
  integer :: u_evolution, ns

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  if(myid==0) write(*,*) 'numprocs=', numprocs, 'myid=', myid
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  call cpu_time(tarray(1))  !cpu_time is a f95 intrinsic subroutine

  call read_parameters()
  call ppinit(ntube, tube_comm, grid_comm)
  t1=mpi_wtime() !measure the wall time
  GCLR=INT(myid/ntube)
  TCLR=MOD(myid,ntube)
  mpol2=numprocs/ntube !poloidal grids for perturbed field
  if(mod(numprocs,ntube) .ne. 0) then
     write(*,*) "mod(numprocs,ntube) must be zero, please adjust numprocs or ntube"
     goto 1234 !end the job
  endif
  if(mod(mpol-1,mpol2).ne.0) then
     write(*,*) '***, mod(mpol-1,mpol2) must be zero, please adjust poloidal gridpoint number'
     goto 1234 !end the job
  endif
  dtheta2=twopi/mpol2
  my_right=GCLR+1
  my_right2=GCLR+2
  if(my_right==mpol2) my_right=0
  if(my_right2==mpol2) my_right2=0
  if(my_right2==mpol2+1) my_right2=1

  my_left=GCLR-1
  my_left2=GCLR-2
  if(my_left==-1) my_left=mpol2-1
  if(my_left2==-1) my_left2=mpol2-1
  if(my_left2==-2) my_left2=mpol2-2
  !Domain decomposition. A mpi process is responsible for the poloidal range between [theta_start : theta_start+dtheta2]
  theta_start = -pi+GCLR*dtheta2 
  GCLR_cut=0
  GCLR_cut_left = mpol2-1

  call construct_numerical_tokamak_equilibrium() !in cylindrical coordinates
  call construct_magnetic_coordinates() 
  !dtheta is the poloidal angle spacing of the equilibrium grids, dtheta2 is the poloidal angle spacing of grids for perturbations
  multi_eq_cells = NINT(dtheta2/dtheta) 
  ipol_eq = 1+nint((theta_start-theta_1d_array(1))/dtheta) !equilibrium poloidal index of the present MPI processes

  call mapping_cylindrical_to_magnetic_coordinates() !for fk species
  !if(myid.eq.0) call field_lines_analyse()
  call prepare_table_in_mc() 

  if(myid==0 .and. diagnosis .eqv. .true.) call visualize_grid()

  call allocate_fk_particle_array()
  call allocate_gk_particle_array()

  omegan_i=bn*charge_i/mass_i !cyclotron angular frequency in Hz
  tn_i=twopi/omegan_i !time unit used in this program
  vn_i=Ln/tn_i !the value of the normalizing velocity in SI unit m/s
  beta_ni=ni0*mass_i*vn_i**2/two/(bn**2/(two*mu0)) !following the notation used in my notes

  omegan_e(:) = bn*abs(charge_e(:))/mass_e(:) !cyclotron angular frequency in Hz
  tn_e(:) = twopi/omegan_e(:)
  vn_e(:) = ln/tn_e(:)
  tn_main = tn_e(1) !choose which species is used as main ion species
  vu = vn_e(1) 
  charge_sign_e(:) = sign(1._p_,charge_e(:))
  beta_ne = ne0*mass_e*vn_e**2/two/(bn**2/(two*mu0)) !following the notation used in my notes

  if(kstart.eq.1) then
     call load_fk() !spatial location in magnetic coordinates (psi,theta,phi) and then transform to cylindrical coordinates
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,radcor_i,active_i,touch_bdry_i,theta_i,alpha_i)

     do ns=1,nsm
        call load_gk(ns, nm_gk(ns), radcor_e(:,ns), theta_e(:,ns), alpha_e(:,ns), vpar_e(:,ns), &
             & mu_e(:,ns),v_e(:,ns), ptcl_num0_e(:,ns), touch_bdry_e(:,ns))
     enddo
     call initialize_gk_weight()     
     w_unit = ne0(1)*vol/nm_gk(1)
     if(myid==0) write(*,*) 'w_unit=', w_unit
     w_e = w_e / w_unit
     ptcl_num0_e = ptcl_num0_e / w_unit
  endif

  call set_gyro_phase()

  call allocate_field_matrix()

  omega_i_axis=abs(baxis)*charge_i/mass_i
  omega_e_axis=abs(baxis)*abs(charge_e)/mass_e
  rho_i=sqrt(ti0*kev/mass_i)/omega_i_axis
  rho_e=sqrt(te0*kev/mass_e)/omega_e_axis

  if(myid.eq.0) write(*,*) 'bulk ion cyclotron angular frequency, omega_i_axis (MHz)=', omega_i_axis/(1.d6)
  if(myid.eq.0) write(*,*) 'kappa_ti*rho_i=', kappa_ti*rho_i, 'te0=',te0

  dtao_i=dtao_omega_i_axis/omega_i_axis/tn_i !time step in unit of tn_i
  dtao_e(:)=dtao_i*(mass_i/mass_e(:)) !time step in unit of the gyro-period, tn_e=twopi/abs(omegan_e)
  dtao_main=dtao_e(1)
  if(myid.eq.0) write(*,*) 'dt (seconds)=',dtao_i*tn_i, 'dtao_i=dt/tn_i',dtao_i
  if(myid.eq.0) write(*,*) 'dtao_e=dt/tn_e',dtao_e
  if(myid.eq.0) write(*,*) 'Total_evolution_time (seconds)=',dtao_i*tn_i*(kend-kstart+1)

  minor_r_min=minor_r_radcor(radcor_low2)
  minor_r_max=minor_r_radcor(radcor_upp2)
  minor_r_width=minor_r_max-minor_r_min
  if(myid.eq.0) write(*,*) 'minor_r_min, minor_r_max, minor_r_width, center (m)=',&
       & minor_r_min, minor_r_max, minor_r_width, (minor_r_min+minor_r_max)/two
  if(myid.eq.0) write(*,*) 'minor_r_width/rho_i=',minor_r_width/rho_i
  if(myid.eq.0) write(*,*) 'minor_r_width/rho_e(:)=',minor_r_width/rho_e

  if(myid.eq.0) write(*,*) 'first radial sine harmonic: kr*rho_i=',pi/minor_r_width*rho_i
  if(myid.eq.0) write(*,*) 'R0/rho_i=',r_axis/rho_i
  !if(myid==0) write(*,*) radcor_minor_r(0.5d0)
  k_binormal1=nsegment*abs(qfunc(radcor_fixed))/minor_r_radcor(radcor_fixed)*rho_i
  k_binormal2=nsegment*abs(baxis)/(GSpsi_prime*&
       & (grad_psi(1,j_fixed)+grad_psi(mpol/2,j_fixed))/two)*rho_i
  if(myid.eq.0) write(*,*) 'k_binorm*rhoi1=',k_binormal1,'k_binorm*rhoi2=',k_binormal2
  if(myid.eq.0) write(*,*) 'number of radial harmonics that should be included (pi*shear0*k_binormal/kr1)=',&
       & pi*0.84*k_binormal1/(pi/minor_r_width*rho_i)

  if(myid.eq.0) write(*,*) 'omega_star1/twopi (kHz)=', k_binormal1*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_
  if(myid.eq.0) write(*,*) 'omega_star2/twopi (kHz)=', k_binormal2*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_
  !if(myid.eq.0) write(*,*) 'Ion beta=',ti0*kev*ni0/(baxis**2/(two*mu0)), 'ti0 (keV)=',ti0,'ni0 (m^-3)=',ni0
  if(myid.eq.0) write(*,*) 'species beta=',te0*kev*ne0/(baxis**2/(two*mu0)), 'te0 (keV)=',ti0,'ne0 (m^-3)=',ne0
  va0=abs(baxis)/sqrt(mu0*mass_i*ni0)
  block
    real(p_), allocatable :: cs(:)
    allocate(cs(nsm))
    cs=sqrt(te0*kev/mass_i)
    if(myid.eq.0) write(*,*) 'VA0 (10^6m/s)=',va0/(1.d6), 'sound speed (10^6m/s)=',cs/(1.d6), 'VA0/Cs=',va0/cs
  end block
  if(myid.eq.0) write(*,*) 'compressional Alfven wave freq/omega_i=k_binorma*va0/omega_axis=',k_binormal1/rho_i*va0/omega_i_axis
  if(myid.eq.0) write(*,*) 'dt*vte/Lte=',dtao_omega_i_axis/omega_i_axis*(sqrt(te0*kev/mass_e)/(1/kappa_te))
  vr_i_integer=vr_i !vr_i_integer is the projection of velocity at t_{n} to the basis vectors at t_{n}
  vz_i_integer=vz_i 
  vphi_i_integer=vphi_i 

  if(kstart.eq.1) then
     do k=1,nmarker_i !vr_i is initially the the projection of velocity at t_{0} to basis vectors at t_{0}, after the backward pushing, vr_i is the projection of velocity at t_{-1/2} to basis vector at t_{0}
        call backward_half_step_for_boris(sign(1._p_,charge_i),dtao_i,r_i(k),z_i(k),phi_i(k),vr_i(k),vz_i(k),vphi_i(k)) !push only velocity, to set initial condition for the first step of boris algorithm, using mutiple steps (2nd RK)
!!$        call forward_half_step_for_boris(sign(1._p_,charge_i),dtao_i,r_i(k),z_i(k),phi_i(k),&
!!$             & vr_i_integer(k),vz_i_integer(k),vphi_i_integer(k),&
!!$             & r_i_mid(k),z_i_mid(k),phi_i_mid(k),vr_i_integer_mid(k),vz_i_integer_mid(k),vphi_i_integer_mid(k)) !push location half-step and then the velocity is projected onto the new local vector basis
     enddo
     call initialize_weight_i(w_i)
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,&
          & radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !so that markers can be sorted and deposited in magnetic coordinates
     call push_ion_orbit_first_step(dtao_i)
  endif

  call initialize_fft()
  if(kstart.gt.1)  then
     call read_data_for_restarting(kstart)
     !call update_efield_at_right_boundary_of_present_cell() !Electric field value at right-boundary of the present cell is needed when pushing particle weight, smoothing and taking z derivatives
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,&
          & radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !theta of lost markers are not calculated
     call clean_up_lost_markers_fk() !for fully kinetic markers, drop lost markers so that we we have smaller particle arrays
  endif

  id_writing_evolution=numprocs/2 !corresponding to theta=0, roughly the low-field-side midplane
  if(myid.eq.id_writing_evolution) open(newunit=u_evolution, file='evolution.txt')

  t2=mpi_wtime(); if(myid.eq.0) write(*,*) 'mpi_wall_time (seconds)=',t2-t1
  do ns=1,nsm
!!$     do k=1, nm_gk(ns)
!!$        call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,b_mc,&
!!$             & theta_e(k,ns),radcor_e(k,ns),bval0)
!!$     enddo
!!$     v_e(:,ns)=sqrt(vpar_e(:,ns)**2+mu_e(:,ns)*two*bval0)
     call gyro_ring(ns, touch_bdry_e(:,ns), theta_e(:,ns), radcor_e(:,ns),&
          &  alpha_e(:,ns), mu_e(:,ns),  x_ring(:,:,ns), y_ring(:,:,ns)) !output are used by gyro_averaging_field and deposit_electrons
  enddo

  call prepare_poisson_matrix()
  call prepare_ampere_matrix()

  do kt=kstart,kend !time-advancing loop
     apara_s_old = apara_s
     t_omega_i_axis = dtao_omega_i_axis*kt
     if(myid==0) write(*,*) 'time-step No.', kt 
     if(fk_ions_switch==1) then
        call fk_push_first_step()
        call deposit_ions(nmarker_i, active_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid,phi_i_mid,w_i_mid)
     endif
     my_den_e_left=zero;  my_den_e_right=zero !before the deposition, set the array to zero
     my_jpar_e_left=zero; my_jpar_e_right=zero
     do ns=1,nsm
        call gyro_average(ns,nm_gk(ns), radcor_e(:,ns), theta_e(:,ns), alpha_e(:,ns), &
             & x_ring(:,:,ns), y_ring(:,:,ns), touch_bdry_e(:,ns), &
             &  phix_e, phiy_e, phiz_e, ax_e,ay_e,az_e, ahx_e,ahy_e,ahz_e, ah_e)
        call compute_drift(ns,radcor_e(:,ns), theta_e(:,ns), mu_e(:,ns), vpar_e(:,ns), &
             & touch_bdry_e(:,ns),   phix_e(:), phiy_e(:), phiz_e(:), ax_e(:), ay_e(:), az_e(:),&
             & equ_radial_drift, equ_theta_drift, equ_alpha_drift, mirror_force, equ_theta_drift0,&
             & perturbed_radial_drift, perturbed_theta_drift, perturbed_alpha_drift)
        !if((ns==1) .and. (mod(kt-1,10)==0)) call compute_heat_flux(t_omega_i_axis, ns, nm_gk(ns), touch_bdry_e(:,ns), &
        !          & mu_e(:,ns), vpar_e(:,ns), w_e(:,ns), radcor_e(:,ns), bval_e, grad_psi_e, perturbed_radial_drift)
        call push_gk_weight(ns, one_half*dtao_e(ns), touch_bdry_e(:,ns), radcor_e(:,ns),&
             & vpar_e(:,ns), v_e(:,ns), equ_radial_drift, equ_alpha_drift, equ_theta_drift, equ_theta_drift0, &
             & mirror_force,perturbed_radial_drift,&
             & phix_e,phiy_e,phiz_e, ahx_e, ahy_e, ahz_e, ah_e, w_e(:,ns), w_e_mid(:,ns))
        call push_gc_half_step(ns,equ_radial_drift, equ_theta_drift, equ_alpha_drift, mirror_force,&
             & perturbed_radial_drift, perturbed_theta_drift, perturbed_alpha_drift)
        call sort_gk_markers(ns, theta_e_mid(:,ns), step=1)
        call gyro_ring(ns,touch_bdry_e(:,ns), theta_e_mid(:,ns), radcor_e_mid(:,ns), alpha_e_mid(:,ns), &
             & mu_e(:,ns),  x_ring(:,:,ns), y_ring(:,:,ns))
        call deposit_gk(ns, touch_bdry_e(:,ns), theta_e_mid(:,ns), vpar_e_mid(:,ns),&
             & w_e_mid(:,ns), x_ring(:,:,ns), y_ring(:,:,ns))
     enddo
     call apara_s_evolution(apara_s_old(:,:,1), apara_s(:,:,1), 0.5_p_*dtao_main)
     call solve_poisson()
     call solve_ampere()
     !---------------second step of 2nd order Runge-Kutta------------------------
     if(fk_ions_switch==1) then
        call fk_push_second_step()
        call deposit_ions(nmarker_i, active_i,radcor_i,theta_i,alpha_i,phi_i,w_i_star)
     endif
     my_den_e_left=zero; my_den_e_right=zero
     my_jpar_e_left=zero; my_jpar_e_right=zero

     do ns=1,nsm
        call gyro_average(ns,nm_gk(ns), radcor_e_mid(:,ns),theta_e_mid(:,ns), alpha_e_mid(:,ns),&
             & x_ring(:,:,ns), y_ring(:,:,ns), touch_bdry_e(:,ns),  &
             & phix_e, phiy_e, phiz_e, ax_e,ay_e,az_e, ahx_e,ahy_e,ahz_e, ah_e)
        call compute_drift(ns,radcor_e_mid(:,ns),theta_e_mid(:,ns),mu_e(:,ns),vpar_e_mid(:,ns),&
             & touch_bdry_e(:,ns), phix_e, phiy_e, phiz_e, ax_e, ay_e, az_e,&
             & equ_radial_drift, equ_theta_drift, equ_alpha_drift, mirror_force, equ_theta_drift0,&
             & perturbed_radial_drift, perturbed_theta_drift, perturbed_alpha_drift)
        call push_gk_weight(ns, dtao_e(ns), touch_bdry_e(:,ns), radcor_e_mid(:,ns),&
             & vpar_e_mid(:,ns), v_e(:,ns), equ_radial_drift, equ_alpha_drift, equ_theta_drift, equ_theta_drift0, &
             & mirror_force,perturbed_radial_drift,&
             & phix_e,phiy_e,phiz_e, ahx_e, ahy_e, ahz_e, ah_e, w_e(:,ns), w_e(:,ns))
        call push_gc_full_step(ns, equ_radial_drift, equ_theta_drift, equ_alpha_drift, mirror_force, &
             & perturbed_radial_drift, perturbed_theta_drift, perturbed_alpha_drift)
        call sort_gk_markers(ns,theta_e(:,ns),step=2)
        call gyro_ring(ns,touch_bdry_e(:,ns), theta_e(:,ns),radcor_e(:,ns), alpha_e(:,ns), mu_e(:,ns),  &
             & x_ring(:,:,ns), y_ring(:,:,ns))
        call deposit_gk(ns,touch_bdry_e(:,ns), theta_e(:,ns), vpar_e(:,ns), w_e(:,ns), x_ring(:,:,ns),y_ring(:,:,ns))
        !if(mod(kt,500)==0) call count_lost_markers_gk(ns)
     enddo
     call apara_s_evolution(apara_s_old(:,:,1), apara_s(:,:,1), dtao_main)
     call solve_poisson()
     call solve_ampere()
     
     call apara_resplit_and_weight_pullback() !at the end of each time-step

     !if(myid.eq.1) call report(t_omega_i_axis)
     if(myid.eq.id_writing_evolution) call mode_evolution_analysis5(t_omega_i_axis, potential(:,:,1), mtor,nflux2, u_evolution)
     if(TCLR.eq.0 .and. mod((kt-1),iplot_mode_structure).eq.0) then
        call mode_structure_on_xy_plane(kt,GCLR,potential(:,:,1),'epa')
        call mode_structure_on_xz_plane(kt,potential(:,:,1),'epa')
        call mode_structure_on_yz_plane(kt,potential(:,:,1),'epa')
        call mode_structure_on_poloidal_plane(kt,potential(:,:,1))
     endif
  enddo !------------------------------------!main time advancing loop
  if(myid.eq.0) close(u_evolution)

  !call write_data_for_restarting(kend)

!!$  write(*,*) 'myid=',myid,'ntouch_bdry_i,e=',ntouch_bdry_i,ntouch_bdry_e
!!$  call MPI_Allreduce( ntouch_bdry_i,  total_ntouch_bdry_i, 1, MPI_integer, MPI_sum, MPI_Comm_World, ierr)
!!$  call MPI_Allreduce( ntouch_bdry_e,  total_ntouch_bdry_e, 1, MPI_integer, MPI_sum, MPI_Comm_World, ierr)
!!$  if(myid.eq.0) write(*,*) 'total_ntouch_bdry_i,e in all processes=',total_ntouch_bdry_i,total_ntouch_bdry_e, ', fraction=',&
!!$       & real(total_ntouch_bdry_i)/total_nmarker_i,real(total_ntouch_bdry_e)/total_nm_gk


!!$  write(filename,'(a1,i4.4)') 'e',myid
!!$  file_unit=myid+100
!!$  open(file_unit,file=filename)
!!$  do k=1,nm_gk
!!$     write(file_unit,*)  radcor_e(k) ,theta_e(k),rg_e(k),zg_e(k)
!!$  enddo
!!$  close(file_unit)


!!$  write(filename,'(a1,i4.4)') 'k',myid
!!$  file_unit=myid+100
!!$  open(file_unit,file=filename)
!!$  do k=1,nmarker_i
!!$    if(touch_bdry_i(k).eqv..true.)  write(file_unit,*)  radcor_i(k) ,theta_i(k),r_i(k),z_i(k)
!!$!     write(file_unit,*)  vr_i(k),vphi_i(k),vz_i(k)
!!$  enddo
!!$  close(file_unit)

1234 call cpu_time(tarray(2))  !total CPU time used by all threads of each process
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  t2 = mpi_wtime()
  if (myid.eq.0) write(*,*) 'Total CPU time (seconds) of all threads of a process', tarray(2)-tarray(1), &
       & 'Wall time used (seconds) ',  (count2-count1)/count_rate, &
       & 'MPI_wall_time (seconds)=',  t2-t1

  call MPI_FINALIZE(ierr)
end program main



subroutine read_parameters()
  use constants,only:p_
  use control_parameters,only:kstart,kend,dtao_omega_i_axis,niter,ion_spatial_loading_scheme,ion_velocity_loading_scheme,&
       &  gk_spatial_loading_scheme,gk_velocity_loading_scheme,poloidal_angle_type,iplot_mode_structure,&
       & filter_toroidal,filter_radial,radial_harmonics_included, &
       & fk_ions_switch, gk_species_switch, space_charge_switch, fk_nonlinear, gk_nonlinear,  remove_n0,&
       &   diagnosis, ismooth, nh
  use magnetic_coordinates,only:nflux2,mpol,mtor,pfn_inner, pfn_bdry,nsegment
  use domain_decomposition,only: ntube,myid
  use gk_module, only : nsm
  use mpi
  implicit none
  integer:: u,ierr
  namelist/control_nmlt/kstart,kend,dtao_omega_i_axis,niter,ion_spatial_loading_scheme,ion_velocity_loading_scheme,&
       & gk_spatial_loading_scheme,gk_velocity_loading_scheme,iplot_mode_structure,&
       & filter_toroidal,filter_radial,radial_harmonics_included,&
       & poloidal_angle_type,nsegment,nflux2,mpol,mtor,pfn_inner, pfn_bdry,ntube, &
       & fk_ions_switch, gk_species_switch, space_charge_switch,&
       &  fk_nonlinear, gk_nonlinear,nh, remove_n0,diagnosis, ismooth, nsm

  open(newunit=u,file='input.nmlt')
  read(u,control_nmlt)
  close(u)
  if(myid==0)  write(*,control_nmlt)

  if((fk_ions_switch==0) .and. (gk_species_switch==0)) stop "****error: can not use adiabatic responce for both species"
  !write(*,*)    'charge_sign_e=',   charge_sign_e
end subroutine read_parameters
