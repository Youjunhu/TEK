module gk_module
  use constants,only:p_
  implicit none
  save
  integer :: nsm !total number of species, get its value from the input namelist
  integer :: nmmax !set in the subroutine below
  real(p_),allocatable :: mass_e(:), charge_e(:), charge_sign_e(:)
  real(p_),allocatable :: ne0(:), te0(:) !te0 is a typical value of temperature in the compuational region, in flux tube model, it is the temperature at the reference magnetic surface, used in setting the velocity range when loading markers
  real(p_),allocatable :: kappa_ne(:), kappa_te(:)
  real(p_),allocatable :: omega_e_axis(:), omegan_e(:), tn_e(:), dtao_e(:), vn_e(:)
  real(p_),allocatable :: beta_ne(:), rho_e(:)
  integer,allocatable  :: total_nm_gk(:)  ! total number of markers in all the processors.
  integer,allocatable ::  nm_gk(:) !number of markers in a processor, fluctuates with time
  logical,allocatable  :: gk_flr(:)
  integer, parameter ::  gyro_npt=4 !number of points on a gyro-ring
  real(p_) :: w_unit
  real(p_),allocatable :: radcor_e(:,:),theta_e(:,:),alpha_e(:,:) !magnetic coordinates, alpha_e is the generalized toroidal angle
  real(p_),allocatable :: vpar_e(:,:), mu_e(:,:), v_e(:,:)
  real(p_),allocatable :: ptcl_num0_e(:,:), w_e(:,:) !weight of electron markers
  logical,allocatable :: touch_bdry_e(:,:) !indicates whether the orbit of a marker touches the boundary
  real(p_),allocatable :: radcor_e_mid(:,:), theta_e_mid(:,:), alpha_e_mid(:,:), &
       & vpar_e_mid(:,:), w_e_mid(:,:) !vlaues at half-time-step

  !computed on fly, i.e., temporary working arrays and hence, no time subfix (_mid), and no need to sorting
  real(p_),allocatable :: x_ring(:,:,:), y_ring(:,:,:)   
  !reused for species and hence no species index
  real(p_),allocatable :: equ_radial_drift(:), equ_theta_drift(:),equ_alpha_drift(:),mirror_force(:), equ_theta_drift0(:)
  real(p_),allocatable :: perturbed_radial_drift(:), perturbed_alpha_drift(:), perturbed_theta_drift(:)
  real(p_),allocatable :: phix_e(:), phiy_e(:), phiz_e(:), ax_e(:), ay_e(:), az_e(:), &
  & ahx_e(:), ahy_e(:), ahz_e(:), ah_e(:) !gyro-averaged perturbed field, used in the mixed-variable method
  !---for testing
  !integer :: ntouch_bdry_e=0, total_ntouch_bdry_e !numbe of markers that touch the boundary in each process and all processes
  !-----
contains
  subroutine allocate_gk_particle_array()
    use domain_decomposition,only : numprocs, myid
    use magnetic_coordinates, only : nflux2, mpol2, mtor
    integer,allocatable  :: nm_gk_per_cell(:)
    namelist/gk_nmlt/mass_e, charge_e, gk_flr, te0,ne0,kappa_ne,kappa_te,nm_gk_per_cell
    integer :: u

    allocate(mass_e(nsm), charge_e(nsm), charge_sign_e(nsm))
    allocate(ne0(nsm), te0(nsm))
    allocate(kappa_ne(nsm), kappa_te(nsm))
    allocate(omega_e_axis(nsm), omegan_e(nsm), tn_e(nsm), dtao_e(nsm), vn_e(nsm))
    allocate(beta_ne(nsm), rho_e(nsm))
    allocate(nm_gk_per_cell(nsm))
    allocate(total_nm_gk(nsm))
    allocate(nm_gk(nsm))
    allocate(gk_flr(nsm))

!!$    mass_e(:)=[3.3452d-27 , 9.1094d-31]
!!$    charge_e(:)= [1.6022d-19 , -1.6022d-19] !in SI unit coulumb, can be positive or negative
!!$    gk_flr(:) = [.true. , .false.]
!!$    te0(:)= [2.14d0 , 2.14d0] !in kev
!!$    ne0(:)= [4.66d19 , 4.66d19] !unit m^-3
!!$    kappa_ne(:)=[1.3174d0 , 1.3174d0] !in unit Ln^-1
!!$    kappa_te(:)=[4.1317d0 , 4.1317d0] !in unit Ln^-1
!!$    nm_gk_per_cell(:)=[20 ,20]  

    open(newunit=u,file='input.nmlt')
    read(u,gk_nmlt)
    close(u)
    if(myid==0)  write(*,gk_nmlt)

    total_nm_gk(:)=nm_gk_per_cell(:)*nflux2*mpol2*mtor
    nm_gk(:) = total_nm_gk(:)/numprocs !the number of markers initially loaded per processor, Later, actual number of markers per proc will be assigned to nm_gk, the value of which will be differnt for differnt processors and at differnt time
    nmmax=int((maxval(total_nm_gk)/numprocs)*1.6) !marker number in one process fluctuates with time, so choose a large number (fixed) to be safe

    allocate(radcor_e(nmmax, nsm)) 
    allocate(theta_e(nmmax,nsm))
    allocate(alpha_e(nmmax,nsm))

    allocate(v_e(nmmax, nsm))
    allocate(vpar_e(nmmax,nsm))
    allocate(mu_e(nmmax,nsm))
    allocate(w_e(nmmax,nsm)) 
    allocate(w_e_mid(nmmax,nsm)) 
    allocate(ptcl_num0_e(nmmax,nsm))

    allocate(radcor_e_mid(nmmax,nsm)) 
    allocate(theta_e_mid(nmmax,nsm))
    allocate(alpha_e_mid(nmmax,nsm))
    allocate(vpar_e_mid(nmmax,nsm))

    allocate(touch_bdry_e(nmmax,nsm))


    !--------
    allocate(equ_radial_drift(nmmax))
    allocate(equ_theta_drift(nmmax))
    allocate(equ_theta_drift0(nmmax))
    allocate(equ_alpha_drift(nmmax))
    allocate(mirror_force(nmmax))
    allocate(perturbed_radial_drift(nmmax))
    allocate(perturbed_alpha_drift(nmmax))
    allocate(perturbed_theta_drift(nmmax))

    allocate(phix_e(nmmax))
    allocate(phiy_e(nmmax))
    allocate(phiz_e(nmmax))
    allocate(ax_e(nmmax))
    allocate(ay_e(nmmax))
    allocate(az_e(nmmax))
    allocate(ahx_e(nmmax))
    allocate(ahy_e(nmmax))
    allocate(ahz_e(nmmax))
    allocate(ah_e(nmmax))  

    allocate(x_ring(gyro_npt, nmmax, nsm))
    allocate(y_ring(gyro_npt, nmmax, nsm))
    if(myid.eq.0) write(*,*) 'marker number of each gk species=',total_nm_gk
  end subroutine allocate_gk_particle_array

  subroutine sort_gk_markers(ns,theta, step)
    !assign particles to different processors according to their poloidal coordinates theta, using the subroutines provided in pputil
    use constants,only: twopi
    use pputil, only : init_pmove, pmove, ppexit, pmove2, end_pmove
    implicit none
    integer, intent(in) :: ns
    real(p_),intent(in) :: theta(:)
    integer, intent(in) :: step
    integer :: np_old, np_new, ierr

    np_old=nm_gk(ns)
    call init_pmove(theta(:), np_old, twopi, ierr)

    call pmove(radcor_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(theta_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(alpha_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(vpar_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(mu_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(v_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove2(touch_bdry_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ptcl_num0_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(w_e(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit

    !    if (step==1) then !need sorting only at the half time-step
    call pmove(radcor_e_mid(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(theta_e_mid(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(alpha_e_mid(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(vpar_e_mid(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(w_e_mid(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    !   endif
    call end_pmove(ierr)
    nm_gk(ns)=np_new
  end subroutine sort_gk_markers

end module gk_module
