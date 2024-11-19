module fk_module !fully-kinetic for ions
  use constants,only:p_
  implicit none
  save
  real(p_) :: mass_i,charge_i,dtao_i
  real(p_) :: ni0,ti0,kappa_ni,kappa_ti !ti0 is a typical value of ion temperature in the compuational region, in flux tube model, ti0 is the temperature at the reference magnetic surface
  real(p_) :: vt_i,vmin_i,vmax_i
  real(p_) ::normalizing_factor
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

  real(p_),allocatable :: vr_i(:),vz_i(:),vphi_i(:) ! projection of velocity at t_{n-1/2} to the basis cylindrical vector at integer-time-step t_{n}, unit vn_i=ln/tn_i, 
  real(p_),allocatable :: vr_i_old(:),vz_i_old(:),vphi_i_old(:) !temporary working arrays for computing averaing
  real(p_),allocatable :: vr_i_integer_mid(:),vz_i_integer_mid(:),vphi_i_integer_mid(:) !projection of velocity at t_{n} to the basis cylindrical vector at t_{n+1/2}, initial condition for the second boris pusher, unit vn_i=ln/tn_i, 

  real(p_),allocatable :: vr_i_mid(:),vz_i_mid(:),vphi_i_mid(:) !projection of velocity at t_{n+1/2} to the local basis vectors at t_{n+1/2}
  real(p_),allocatable :: vr_i_integer(:),vz_i_integer(:),vphi_i_integer(:) !projection of velocity at t_{n} to the local basis vectors at t_{n}

  real(p_),allocatable :: v_i(:),vpar_i(:),vx_i(:),vy_i(:) !vx is defined by vx=v_dot_grad_x, vy is defined by vy=v_dot_grad_y,note that grad_x and grad_y are not perpendicular to each other
  real(p_),allocatable :: grad_psi_i(:),grad_alpha_i(:), grad_psi_dot_grad_alpha_i(:),bval_i(:)
  real(p_),allocatable :: v_i_mid(:),vpar_i_mid(:),vx_i_mid(:),vy_i_mid(:)
  real(p_),allocatable :: grad_psi_i_mid(:),grad_alpha_i_mid(:), grad_psi_dot_grad_alpha_i_mid(:),bval_i_mid(:)

  logical,allocatable :: touch_bdry_i(:),active_i(:) !indicates whether the orbit of a marker touches the boundary
  logical,allocatable :: touch_bdry_i_mid(:),active_i_mid(:)
contains
  subroutine allocate_fk_particle_array()
    use domain_decomposition,only: numprocs, myid
    use magnetic_coordinates, only : nflux2, mpol2,mtor
    namelist/fk_nmlt/mass_i,charge_i,ni0,ti0,kappa_ni,kappa_ti,nmarker_i_per_cell
    integer:: fixed_large_size, u

     open(newunit=u,file='input.nmlt')
     read(u,fk_nmlt)
     close(u)
     if(myid==0)  write(*,fk_nmlt)

    
    total_nmarker_i=nmarker_i_per_cell*nflux2*mpol2*mtor
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

  end subroutine allocate_fk_particle_array

end module fk_module
