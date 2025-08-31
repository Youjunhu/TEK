program main
  use constants, only: p_, twopi, pi, kev, two, mu0, one_half, zero
  use control_parameters, only: kstart, kend, dt_omega_i_axis, iplot_mode_structure, &
       &  fk_switch, diagnosis
  use normalizing, only: bn,ln, dtao_main, qu, tu, vu, nu
  use magnetic_coordinates, only: mpol,nrad, xlow, xupp, nsegment, dtheta, vol, grad_psi, &
       & mpol2, mtor, zgrid, tor_shift_mc, toroidal_range, GSpsi_prime
  use magnetic_field, only: qfunc, pfn_func
  use func_in_mc, only: minor_r_radcor
  use table_in_mc, only: prepare_table_in_mc
  use radial_module, only: baxis, r_axis, minor_a, radcor_fixed, j_fixed
  use fk_module, only: mass_i,charge_i,dtao_fk, nmarker_i, ni0,ti0,kappa_ti, &
       & nmarker_i_per_cell, r_i,z_i,phi_i, radcor_i,theta_i,alpha_i,tor_shift_i,ps_vol_i,&
       & touch_bdry_i,active_i, touch_bdry_i_mid,active_i_mid, &
       & radcor_i_mid,theta_i_mid,alpha_i_mid, &
       & vr_i,vz_i,vphi_i, ntouch_bdry_i, total_ntouch_bdry_i, &
       & vr_i_integer,vz_i_integer,vphi_i_integer,&
       & vr_i_mid,vphi_i_mid,vz_i_mid, &
       & r_i_mid,z_i_mid,phi_i_mid,vr_i_integer_mid,vz_i_integer_mid,vphi_i_integer_mid, &
       & grad_psi_i_mid,grad_alpha_i_mid,grad_psi_dot_grad_alpha_i_mid,bval_i_mid, &
       & v_i_mid,vpar_i_mid,vx_i_mid,vy_i_mid, &
       & w_i,w_i_star,w_i_mid, initialize_fk, tn_fk, vn_fk, omegan_fk
  use gk_module, only : nsm, nmmax, nm_gk, mass_gk, charge_gk, tgk0, ngk0, &
       & dtao_gk, vn_gk, w_unit, touch_bdry_gc, &
       & xgc, zgc, ygc,vpar_gk,mu_gk, v_gk, ps_vol_gk, ptcl_num0_gk, w_gk, w_gk_mid, &
       & xgc_mid, zgc_mid, ygc_mid, vpar_gk_mid, x_ring, y_ring, z_ring, &
       & initialize_gk, sort_gk_markers
  use load_gk_mod, only: load_gk
  use gyro_ring_mod, only : set_gyro_phase, gyro_ring
  use gyro_average_mod, only : gyro_average
  use drift, only : compute_drift
  use gk_trajectory_pusher, only: push_gc
  use gk_weight_pusher, only: push_gk_weight
  use perturbation_field, only: allocate_field_matrix, potential, phix, phiy, phiz, &
       &  apara, apara_h, apara_s, apara_s_old, ax, ay, az, ahx, ahy, ahz
  use domain_decomposition, only: myid,numprocs,tube_comm,grid_comm,ntube,gclr,tclr,GCLR_cut, GCLR_cut_left,&
       & dtheta2,theta_start,my_right,my_left, my_right2, my_left2, multi_eq_cells, ipol_eq, dvol
  use misc, only: calculate_dvol
  use pputil, only: ppinit
  use fk_particle_coordinates_transform_module, only: compute_particle_magnetic_coordinates, &
       & clean_up_lost_markers_fk, count_lost_markers_fk
  use mode_structure
  use diagnosis_mod, only: report, mode_evolution_analysis2, mode_evolution_analysis6, compute_heat_flux
  use deposit_fk_module, only: deposit_fk
  use deposit_gk_module, only: deposit_gk
  use FFTW3
  use spectrum_diagnostic
  use poisson, only: prepare_poisson_matrix, solve_poisson
  use ampere, only: prepare_ampere_matrix, solve_ampere, apara_s_evolution, &
       & apara_resplit_and_weight_pullback
  use restart_mod
  use push_ion_weight_module
  use initial_half_step_for_boris
  use initialization_mod
  use interpolate_module, only: linear_2d_interpolate
  use mpi
  implicit none

  integer :: ierr, id_writing_evolution, k, kt, ns
  integer :: count1,count2, count3, count_rate, count_max
  integer :: file_unit_i, file_unit_e, file_unit, u_evolution
  real(p_) :: stime, omega_i_axis, rho_i, k_binormal1, k_binormal2, beta_ni
  real(p_) :: minor_r_min, minor_r_max, minor_r_width, va0
  real(p_) :: t1, t2, tarray(2) !store the cpu clock time
  real(p_), allocatable :: xdrift0(:), ydrift0(:), zdrift0(:), zdrift00(:), mirror_force(:)
  real(p_), allocatable :: xdrift1(:), ydrift1(:), zdrift1(:)
  real(p_), allocatable :: phix_ga(:), phiy_ga(:), phiz_ga(:), ax_ga(:), ay_ga(:), az_ga(:)
  real(p_), allocatable :: ahx_ga(:), ahy_ga(:), ahz_ga(:), ah_ga(:) !gyro-averaged perturbations
  real(p_), dimension(:,:), allocatable :: density_left, density_right !gk density
  real(p_), dimension(:,:), allocatable :: jpar_left, jpar_right !gk parallel current

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  if(myid==0) write(*,*) 'numprocs=', numprocs, 'myid=', myid
  CALL SYSTEM_CLOCK(count1, count_rate, count_max)
  call cpu_time(tarray(1))  !f95 intrinsic subroutine

  call read_parameters()
  call ppinit(ntube, tube_comm, grid_comm)

  t1 = mpi_wtime() !measure the wall time
  GCLR = INT(myid/ntube)
  TCLR = MOD(myid,ntube)
  mpol2 = numprocs/ntube !poloidal grids for perturbed field
  if(mod(numprocs,ntube) .ne. 0) then
     write(*,*) "mod(numprocs,ntube) must be zero, please adjust numprocs or ntube"
     goto 1234 !end the job
  endif
  if(mod(mpol-1,mpol2) .ne. 0) then
     write(*,*) '***, mod(mpol-1,mpol2) must be zero, please adjust poloidal gridpoint number'
     goto 1234 !end the job
  endif
  dtheta2=twopi/mpol2 !the poloidal angle spacing of grids for perturbations
  my_right = GCLR+1
  my_right2 = GCLR+2
  if(my_right==mpol2) my_right=0
  if(my_right2==mpol2) my_right2=0
  if(my_right2==mpol2+1) my_right2=1

  my_left=GCLR-1
  my_left2=GCLR-2
  if(my_left==-1) my_left=mpol2-1
  if(my_left2==-1) my_left2=mpol2-1
  if(my_left2==-2) my_left2=mpol2-2
  !Domain decomposition.
  !A mpi process is responsible for the poloidal range [theta_start : theta_start+dtheta2]
  theta_start = -pi+GCLR*dtheta2 
  GCLR_cut=0
  GCLR_cut_left = mpol2-1

  call read_and_process_equilibrium() !in cylindrical coordinates
  call construct_magnetic_coordinates()

  !dtheta is the poloidal angle spacing of the equilibrium grids,
  multi_eq_cells = NINT(dtheta2/dtheta)
  !equilibrium poloidal index of the present MPI processes
  ipol_eq = 1+nint((theta_start-zgrid(1))/dtheta)
  if(myid==0) write(*,*) 'multi_eq_cells=', multi_eq_cells
  call calculate_dvol(multi_eq_cells, dvol)
  call mapping_cylindrical_to_magnetic_coordinates() !for fk species
  !if(myid.eq.0) call field_lines_analyse()
  call prepare_table_in_mc() 

  if(myid==0 .and. diagnosis .eqv. .true.) call visualize_grid()
  minor_r_min = minor_r_radcor(xlow)
  minor_r_max = minor_r_radcor(xupp)
  minor_r_width = minor_r_max-minor_r_min
  if(myid==0) write(*,*) 'minor_r_min, minor_r_max, minor_r_width, center (m)=', &
       & minor_r_min, minor_r_max, minor_r_width, (minor_r_min+minor_r_max)/two

  if(fk_switch==1) call initialize_fk()
  call initialize_gk(baxis, minor_a, dt_omega_i_axis) !set radial profiles etc.
  qu = charge_gk(1)
  tu = tgk0(1)
  vu = vn_gk(1)
  nu = ngk0(1)
  dtao_main = dtao_gk(1)
  if(myid==0) write(*,*) 'vu/vthermal=', vu/sqrt(2*tgk0(1)*kev/mass_gk(1))

  omegan_fk = bn*charge_i/mass_i !cyclotron angular frequency in Hz
  tn_fk = twopi/omegan_fk !time unit used in this program
  vn_fk = Ln/tn_fk !the value of the normalizing velocity in SI unit m/s
  beta_ni = ni0*mass_i*vn_fk**2/two/(bn**2/(two*mu0)) !following the notation used in my notes

  if(fk_switch==1) then
     call load_fk()
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,radcor_i,active_i,touch_bdry_i,theta_i,alpha_i)
  endif

  if(kstart==0) then
     do ns = 1, nsm
        call load_gk(ns, nm_gk(ns), xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns), &
             & mu_gk(:,ns),v_gk(:,ns), ps_vol_gk(:,ns), ptcl_num0_gk(:,ns), touch_bdry_gc(:,ns))
     enddo

     call initialize_gk_weight()     
     w_unit = ngk0(1)*vol/nm_gk(1)
     if(myid==0) write(*,*) 'w_unit=', w_unit
     w_gk = w_gk / w_unit
     ptcl_num0_gk = ptcl_num0_gk / w_unit
  endif

  call set_gyro_phase()
  call allocate_field_matrix()

  omega_i_axis = abs(baxis)*charge_gk(1)/mass_gk(1)
  rho_i = sqrt(tgk0(1)*kev/mass_gk(1))/omega_i_axis
  dtao_fk = dt_omega_i_axis/omega_i_axis/tn_fk !time step in unit of tn_fk

  if(myid.eq.0) write(*,*) 'baxis*charge_gk/mass_gk (MHz): ', baxis*charge_gk/mass_gk/10**6
  !if(myid.eq.0) write(*,*) 'dt (seconds)=',dtao_fk*tn_fk, 'dtao_fk=dt/tn_fk',dtao_fk
  !if(myid.eq.0) write(*,*) 'Total_evolution_time (seconds)=',dtao_fk*tn_fk*(kend-kstart)

  if(myid.eq.0) write(*,*) 'first radial sine harmonic: kr*rho_i=',pi/minor_r_width*rho_i
  if(myid.eq.0) write(*,*) 'R0/rho_i=',r_axis/rho_i
  !if(myid==0) write(*,*) radcor_minor_r(0.5d0)
  k_binormal1=nsegment*abs(qfunc(radcor_fixed))/minor_r_radcor(radcor_fixed)*rho_i
  k_binormal2=nsegment*abs(baxis)/(GSpsi_prime*&
       & (grad_psi(1,j_fixed)+grad_psi(mpol/2,j_fixed))/two)*rho_i
  if(myid.eq.0) write(*,*) 'k_binorm*rhoi1=',k_binormal1,'k_binorm*rhoi2=',k_binormal2
  if(myid.eq.0) write(*,*) 'number of radial harmonics that should be included (pi*shear0*k_binormal/kr1)=',&
       & pi*0.84*k_binormal1/(pi/minor_r_width*rho_i)

  !if(myid.eq.0) write(*,*) 'omega_star1/twopi (kHz)=', k_binormal1*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_
  !if(myid.eq.0) write(*,*) 'omega_star2/twopi (kHz)=', k_binormal2*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_

  if(myid==0) write(*,*) 'species beta=', tgk0*kev*ngk0/(baxis**2/(two*mu0))
  if(myid==0) write(*,*) 'tgk0 (keV)=', tgk0(:)
  if(myid==0) write(*,*) 'ngk0 (m^-3)=', ngk0(:)
  va0=abs(baxis)/sqrt(mu0*mass_gk(1)*ngk0(1))
  block
    real(p_), allocatable :: cs(:)
    allocate(cs(nsm))
    cs = sqrt(tgk0*kev/mass_gk)
    if(myid==0) write(*,*) 'VA0 (10^6m/s)=',va0/(1.d6)
    if(myid==0) write(*,*) 'sound speed (10^6m/s)=',cs/(1.d6)
    if(myid==0) write(*,*) 'VA0/Cs=',va0/cs
    if(myid==0) write(*,*) 'R0/Cs (seconds) =', r_axis/cs(:)
  end block
  if(myid==0) write(*,*) 'compressional Alfven wave freq/omega_i=k_binorma*va0/omega_axis=',k_binormal1/rho_i*va0/omega_i_axis

  if(kstart==0 .and. fk_switch==1) then
     vr_i_integer = vr_i !vr_i_integer is the projection of velocity at t_{n} to the basis vectors at t_{n}
     vz_i_integer = vz_i 
     vphi_i_integer = vphi_i 

     do k=1,nmarker_i !vr_i is initially the the projection of velocity at t_{0} to basis vectors at t_{0}, after the backward pushing, vr_i is the projection of velocity at t_{-1/2} to basis vector at t_{0}
        call backward_half_step_for_boris(sign(1._p_,charge_i),dtao_fk,r_i(k),z_i(k),phi_i(k),vr_i(k),vz_i(k),vphi_i(k)) !push only velocity, to set initial condition for the first step of boris algorithm, using mutiple steps (2nd RK)
!!$        call forward_half_step_for_boris(sign(1._p_,charge_i),dtao_fk,r_i(k),z_i(k),phi_i(k),&
!!$             & vr_i_integer(k),vz_i_integer(k),vphi_i_integer(k),&
!!$             & r_i_mid(k),z_i_mid(k),phi_i_mid(k),vr_i_integer_mid(k),vz_i_integer_mid(k),vphi_i_integer_mid(k)) !push location half-step and then the velocity is projected onto the new local vector basis
     enddo
     call initialize_weight_i(w_i)
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,&
          & radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !so that markers can be sorted and deposited in magnetic coordinates
     call push_ion_orbit_first_step(dtao_fk)
  endif

  call initialize_fft()
  if(kstart > 0)  then
     call read_data_for_restarting(kstart)
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,&
          & radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !theta of lost markers are not calculated
     call clean_up_lost_markers_fk() !for fully kinetic markers, drop lost markers so that we we have smaller particle arrays
  endif

  id_writing_evolution = numprocs/2 !corresponding to theta=0, roughly the low-field-side midplane
  if(myid.eq.id_writing_evolution) open(newunit=u_evolution, file='evolution.txt')

  t2=mpi_wtime(); if(myid.eq.0) write(*,*) 'wall_time (seconds) used before calculating polarization matrix=', t2-t1
  call prepare_poisson_matrix()
  call prepare_ampere_matrix()
  t2=mpi_wtime(); if(myid.eq.0) write(*,*) 'wall_time (seconds) used before main time loop=', t2-t1

  allocate(phix_ga(nmmax))
  allocate(phiy_ga(nmmax))
  allocate(phiz_ga(nmmax))
  allocate(ax_ga(nmmax))
  allocate(ay_ga(nmmax))
  allocate(az_ga(nmmax))
  allocate(ahx_ga(nmmax))
  allocate(ahy_ga(nmmax))
  allocate(ahz_ga(nmmax))
  allocate(ah_ga(nmmax)) 
  
  allocate(xdrift0(nmmax), ydrift0(nmmax), zdrift0(nmmax))
  allocate(xdrift1(nmmax), ydrift1(nmmax), zdrift1(nmmax))
  allocate(zdrift00(nmmax))
  allocate(mirror_force(nmmax))

  allocate(density_left(mtor,nrad), density_right(mtor,nrad)) 
  allocate(jpar_left(mtor,nrad), jpar_right(mtor,nrad))

  do ns = 1, nsm
     call gyro_ring(ns, touch_bdry_gc(:,ns), xgc(:,ns), ygc(:,ns), zgc(:,ns), &
          & mu_gk(:,ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns)) 
  enddo

  apara_s = 0
  do kt = kstart+1, kend !time-advancing loop
     apara_s_old = apara_s
     stime = dt_omega_i_axis*kt/omega_i_axis
     if(myid==0 .and. mod(kt,10)==0) write(*,*) 'time-step No.', kt 
     if(fk_switch==1) then
        call fk_push_first_step()
        call deposit_fk(nmarker_i, active_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid,phi_i_mid,w_i_mid)
     endif
     density_left = zero;  density_right = zero !before the deposition, set the array to zero
     jpar_left = zero; jpar_right = zero
     do ns = 1, nsm
        call gyro_average(ns, nm_gk(ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns), &
             & touch_bdry_gc(:,ns), phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga)
        call compute_drift(ns, xgc(:,ns), zgc(:,ns), mu_gk(:,ns), vpar_gk(:,ns), &
             & touch_bdry_gc(:,ns), phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, &
             & xdrift0, zdrift0, ydrift0, mirror_force, zdrift00, xdrift1, zdrift1, ydrift1)
        if(mod(kt-1,10)==0) call compute_heat_flux(stime, ns, nm_gk(ns), touch_bdry_gc(:,ns), &
             & mu_gk(:,ns), vpar_gk(:,ns), w_gk(:,ns), xgc(:,ns), zgc(:,ns), xdrift1(:))

        call push_gk_weight(ns, nm_gk(ns), one_half*dtao_gk(ns), touch_bdry_gc(:,ns), xgc(:,ns),&
             & vpar_gk(:,ns), v_gk(:,ns), xdrift0, ydrift0, zdrift0, zdrift00, &
             & mirror_force, xdrift1, ydrift1, zdrift1, &
             & phix_ga, phiy_ga, phiz_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga, w_gk(:,ns), w_gk_mid(:,ns))
        call push_gc(ns, one_half*dtao_gk(ns), xdrift0, zdrift0, ydrift0, mirror_force, &
             & xdrift1, zdrift1, ydrift1, xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns),&
             & xgc_mid(:,ns), zgc_mid(:,ns), ygc_mid(:,ns), vpar_gk_mid(:,ns), touch_bdry_gc(:,ns))

        call sort_gk_markers(ns, zgc_mid(:,ns), step=1)
        call gyro_ring(ns, touch_bdry_gc(:,ns), xgc_mid(:,ns), ygc_mid(:,ns),  zgc_mid(:,ns), &
             & mu_gk(:,ns),  x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns))
        call deposit_gk(ns, touch_bdry_gc(:,ns), vpar_gk_mid(:,ns), w_gk_mid(:,ns), &
             & x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns), &
             & density_left, density_right, jpar_left, jpar_right)
     enddo

     call solve_poisson(density_left, density_right, potential, phix, phiy, phiz)
     call apara_s_evolution(apara_s_old(:,:,1), apara_s(:,:,1), phiz, 0.5_p_*dtao_main)
     call solve_ampere(1, jpar_left, jpar_right, apara_s, apara_h, apara, ax, ay, az, ahx, ahy, ahz)

     !---------------second step of 2nd order Runge-Kutta------------------------
     if(fk_switch==1) then
        call fk_push_second_step()
        call deposit_fk(nmarker_i, active_i, radcor_i, theta_i, alpha_i, phi_i, w_i_star)
     endif
     density_left = zero; density_right = zero
     jpar_left = zero; jpar_right = zero
     do ns = 1, nsm
        call gyro_average(ns, nm_gk(ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns), &
             & touch_bdry_gc(:,ns), phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga)
        call compute_drift(ns, xgc_mid(:,ns), zgc_mid(:,ns), mu_gk(:,ns), vpar_gk_mid(:,ns), &
             & touch_bdry_gc(:,ns), phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, &
             & xdrift0, zdrift0, ydrift0, mirror_force, zdrift00, xdrift1, zdrift1, ydrift1)
        call push_gk_weight(ns, nm_gk(ns), dtao_gk(ns), touch_bdry_gc(:,ns), xgc_mid(:,ns), &
             & vpar_gk_mid(:,ns), v_gk(:,ns), xdrift0, ydrift0, zdrift0, zdrift00, &
             & mirror_force, xdrift1, ydrift1, zdrift1, &
             & phix_ga, phiy_ga, phiz_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga, w_gk(:,ns), w_gk(:,ns))
        call push_gc(ns, dtao_gk(ns), xdrift0, zdrift0, ydrift0, mirror_force, &
             & xdrift1, zdrift1, ydrift1, xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns), &
             & xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns), touch_bdry_gc(:,ns))

        call sort_gk_markers(ns, zgc(:,ns), step=2)
        call gyro_ring(ns, touch_bdry_gc(:,ns), xgc(:,ns), ygc(:,ns), zgc(:,ns), &
             & mu_gk(:,ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns))
        call deposit_gk(ns, touch_bdry_gc(:,ns), vpar_gk(:,ns), w_gk(:,ns), &
             & x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns), &
             density_left, density_right,  jpar_left, jpar_right)
        !if(mod(kt,500)==0) call count_lost_markers_gk(ns)
     enddo

     call solve_poisson(density_left, density_right, potential, phix, phiy, phiz)
     call apara_s_evolution(apara_s_old(:,:,1), apara_s(:,:,1), phiz, dtao_main)
     call solve_ampere(2, jpar_left, jpar_right, apara_s, apara_h, apara, ax, ay, az, ahx, ahy, ahz)
     call apara_resplit_and_weight_pullback(w_gk, apara_s, apara_h, ahx, ahy, ahz) !at the end of each time-step

     if(myid==id_writing_evolution .and. mod(kt,10)==0) then
        call mode_evolution_analysis6(stime, potential(:,:,1), mtor, nrad, u_evolution)
     endif
     if(TCLR==0 .and. mod((kt-1),iplot_mode_structure)==0) then
        call mode_structure_on_xy_plane(kt,GCLR,potential(:,:,1),'epa')
        call mode_structure_on_xz_plane(kt,potential(:,:,1),'epa')
        call mode_structure_on_yz_plane(kt,potential(:,:,1),'epa')
        call mode_structure_on_poloidal_plane(kt, potential(:,:, :), 'Phi')
        call mode_structure_on_poloidal_plane(kt, apara(:,:, :), 'Apara')
     endif
  enddo !-----------time advancing loop---------
  if(myid.eq.0) close(u_evolution)

  !call write_data_for_restarting(kend)

1234 call cpu_time(tarray(2))  !total CPU time used by all threads of each process
  CALL SYSTEM_CLOCK(count2, count_rate, count_max)
  t2 = mpi_wtime()
  if (myid.eq.0) write(*,*) 'Total CPU time (seconds) of all threads of a process', tarray(2)-tarray(1), &
       & 'Wall time used (seconds) ',  (count2-count1)/count_rate, &
       & 'MPI_wall_time (seconds)=',  t2-t1

  call MPI_FINALIZE(ierr)
end program main


subroutine read_parameters()
  use constants, only: p_
  use control_parameters, only: kstart,kend,dt_omega_i_axis, &
       &  poloidal_angle_type,iplot_mode_structure,&
       & filter_toroidal,filter_radial, &
       & fk_switch, space_charge_switch, adiabatic_electrons, &
       &   diagnosis, ismooth, nh_min, nh_max
  use magnetic_coordinates, only: nrad,mpol,mtor,pfn_inner, pfn_bdry,nsegment
  use domain_decomposition, only: ntube,myid
  use gk_module, only : nsm
  use mpi
  implicit none
  integer:: u,ierr
  namelist/control_nmlt/kstart,kend,dt_omega_i_axis, &
       & iplot_mode_structure,&
       & filter_toroidal,filter_radial, &
       & poloidal_angle_type,nsegment,nrad,mpol,mtor,pfn_inner, pfn_bdry,ntube, &
       & fk_switch, space_charge_switch, adiabatic_electrons, &
       &  nh_min, nh_max, diagnosis, ismooth, nsm

  open(newunit=u,file='input.nmlt')
  read(u,control_nmlt)
  close(u)
  if(myid==0)  write(*,control_nmlt)
end subroutine read_parameters
