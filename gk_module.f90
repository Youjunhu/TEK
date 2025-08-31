module gk_module
  use constants, only: p_, twopi, elementary_charge, kev
  use normalizing, only : ln, bn
  implicit none
  save
  integer :: nsm !total number of species, get its value from the input namelist
  integer :: nmmax !maximal number of markers in a processor
  integer, allocatable ::  nm_gk(:) !actual number of markers in a processor, fluctuates with time
  integer, allocatable  :: total_nm_gk(:)  !total number of markers (in all the processors) for a species
  integer  :: gk_spatial_loading_scheme, gk_velocity_loading_scheme 
  integer :: gk_nonlinear
  real(p_), allocatable :: mass_gk(:), charge_gk(:), charge_sign_gk(:)
  real(p_), allocatable :: tgk0(:) !in kev, a typical value of temperature in the compuational region
  real(p_), allocatable :: ngk0(:) !in m^(-3) a typical value of density in the compuational region
  real(p_), allocatable :: vn_gk(:), dtao_gk(:)
  logical, allocatable :: gk_flr(:)
  integer, parameter ::  gyro_npt = 4 !number of points on a gyro-ring
  real(p_) :: w_unit
  real(p_), allocatable :: xgc(:,:), zgc(:,:), ygc(:,:) !magnetic coordinates of guiding centers
  real(p_), allocatable :: vpar_gk(:,:), w_gk(:,:) !weight of gk markers

  real(p_), allocatable :: xgc_mid(:,:), zgc_mid(:,:), ygc_mid(:,:)
  real(p_), allocatable :: vpar_gk_mid(:,:), w_gk_mid(:,:) !vlaues at half-time-step
  
  real(p_), allocatable :: mu_gk(:,:), v_gk(:,:)
  real(p_), allocatable :: ps_vol_gk(:,:), ptcl_num0_gk(:,:)
  logical, allocatable :: touch_bdry_gc(:,:) !indicates whether a marker touches the boundary
  !computed on fly, i.e., temporary working arrays. Hence, no time subfix (_mid), and no need to sort
  real(p_), allocatable :: x_ring(:,:,:), y_ring(:,:,:), z_ring(:,:,:)   
  !---for testing
  !integer :: ntouch_bdry_gc=0, total_ntouch_bdry_gc !numbe of markers that touch the boundary in each process and all processes

contains
  subroutine initialize_gk(baxis, minor_a, dt_omega_i_axis)
    use domain_decomposition,only : numprocs, myid
    use magnetic_coordinates, only : nrad, mpol2, mtor, pfn, pfn_inner, pfn_bdry
    use gk_radial_profiles, only : initialize_gk_radial_profiles
    use gk_profile_funcs, only: gkt_func, gkn_func, gkdtdx_func, gkdndx_func
    real(p_), intent(in) :: baxis, minor_a, dt_omega_i_axis
    real(p_), allocatable ::  omega_gk_axis(:), rho_gk(:)
    character(100) :: density_file(nsm), density_radcor(nsm), temperature_file(nsm), temperature_radcor(nsm)
    real(p_) :: density_unit(nsm), temperature_unit(nsm), pfn0

    integer,allocatable  :: nm_gk_per_cell(:)
    namelist/gk_nmlt/ mass_gk, charge_gk, gk_flr, nm_gk_per_cell, density_file, density_unit, &
         & density_radcor, temperature_file, temperature_unit, temperature_radcor, &
         & gk_spatial_loading_scheme, gk_velocity_loading_scheme, gk_nonlinear
    integer :: u, i, k

    allocate(mass_gk(nsm), charge_gk(nsm), charge_sign_gk(nsm))
    allocate(ngk0(nsm), tgk0(nsm))
    allocate(omega_gk_axis(nsm), dtao_gk(nsm), vn_gk(nsm))
    allocate(rho_gk(nsm))
    allocate(nm_gk_per_cell(nsm))
    allocate(total_nm_gk(nsm))
    allocate(nm_gk(nsm))
    allocate(gk_flr(nsm))

    open(newunit=u,file='input.nmlt')
    read(u,gk_nmlt)
    close(u)
    if(myid==0)  write(*,gk_nmlt)
    charge_gk = charge_gk*elementary_charge
    charge_sign_gk(:) = sign(1._p_, charge_gk(:))
    vn_gk(:) = ln/(twopi/(bn*abs(charge_gk(:))/mass_gk(:)))
    omega_gk_axis=abs(baxis*charge_gk)/mass_gk

    total_nm_gk(:) = nm_gk_per_cell(:) * nrad * mpol2 * mtor
    nm_gk(:) = total_nm_gk(:)/numprocs !the number of markers initially loaded per processor
    !Later, nm_gk will be set to the actual number of markers per proc 
    !the value of which will be differnt for differnt processors and at differnt time

    !marker number in one process fluctuates with time, so choose a large number (fixed) to be safe
    nmmax = int((maxval(total_nm_gk)/numprocs)*1.3) 

    allocate(xgc(nmmax, nsm)) 
    allocate(zgc(nmmax,nsm))
    allocate(ygc(nmmax,nsm))

    allocate(v_gk(nmmax, nsm))
    allocate(vpar_gk(nmmax, nsm))
    allocate(mu_gk(nmmax, nsm))
    allocate(w_gk(nmmax, nsm)) 
    allocate(w_gk_mid(nmmax, nsm)) 
    allocate(ptcl_num0_gk(nmmax, nsm))
    allocate(ps_vol_gk(nmmax, nsm))

    allocate(xgc_mid(nmmax,nsm)) 
    allocate(zgc_mid(nmmax,nsm))
    allocate(ygc_mid(nmmax,nsm))
    allocate(vpar_gk_mid(nmmax,nsm))
    allocate(x_ring(gyro_npt, nmmax, nsm))
    allocate(y_ring(gyro_npt, nmmax, nsm))
    allocate(z_ring(gyro_npt, nmmax, nsm))
    allocate(touch_bdry_gc(nmmax,nsm))
    if(myid.eq.0) write(*,*) 'marker number of each gk species=',total_nm_gk

    call initialize_gk_radial_profiles(nsm, density_file, density_unit, density_radcor, &
         &   temperature_file, temperature_unit, temperature_radcor)
    

    block !diagnostic information
      use magnetic_field, only: tfn_func_pfn
      use func_in_mc, only: minor_r_radcor, minor_r_prime
      real(p_) :: tfn, c
      if(myid==0) then 
         open(newunit=u,file='profiles.txt') 
         do i=1, size(pfn)
            tfn= tfn_func_pfn(pfn(i))
            c = minor_r_prime(pfn(i))
            write(u,'(80ES18.6E4)') pfn(i), tfn, minor_r_radcor(pfn(i)), &
                 & (gkt_func(pfn(i),k), k=1,nsm), &
                 & (gkn_func(pfn(i),k), k=1,nsm), &
                 & (gkdndx_func(pfn(i),k)/c, k=1,nsm), &
                 & (gkdtdx_func(pfn(i),k)/c, k=1,nsm)
         enddo
         close(u)
      endif
    endblock
    
    do i = 1, nsm
       !pfn0 =  3.915686E-1
       pfn0 =  0.1
       tgk0(i) = gkt_func(pfn0, i) !typical value of temperature a species
       ngk0(i) = gkn_func(pfn0, i) !typical value of number density for a apecies
    enddo
    rho_gk = sqrt(tgk0*kev/mass_gk)/omega_gk_axis
    if(myid==0) write(*,*) 'a/rho_gk(:)=',minor_a/rho_gk
    !time step in unit of the gyro-period twopi/abs(omegan_gk(:)):
    dtao_gk(:) = dt_omega_i_axis*(bn*abs(charge_gk(:))/mass_gk(:))/omega_gk_axis(1)/twopi
    if(myid==0) write(*,*) 'dt/tn_gk=',dtao_gk
    if(myid==0) write(*,*) 'tgk0=', tgk0
  end subroutine initialize_gk

  subroutine sort_gk_markers(ns, theta, step)
    !assign particles to different processors according to their poloidal coordinates theta, using the subroutines provided in pputil
    use constants,only: twopi
    use pputil, only : init_pmove, pmove, ppexit, pmove2
    implicit none
    integer, intent(in) :: ns
    real(p_),intent(in) :: theta(:)
    integer, intent(in) :: step
    integer :: np_old, np_new, ierr

    np_old=nm_gk(ns)
    call init_pmove(theta(:), np_old, twopi, ierr)

    call pmove(xgc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(zgc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ygc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(vpar_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(mu_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(v_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove2(touch_bdry_gc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ps_vol_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit

    call pmove(ptcl_num0_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(w_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit

!    if (step==1) then !need sorting only at the half time-step
       call pmove(xgc_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(zgc_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(ygc_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(vpar_gk_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(w_gk_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
 !   endif

    nm_gk(ns)=np_new
  end subroutine sort_gk_markers

end module gk_module
