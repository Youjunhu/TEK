module constants
  implicit none
  save
  integer,parameter:: p_=kind(1.0d0) !precision
  real(p_),parameter :: coulomb_log=15._p_ !assumed to be a constant
  real(p_),parameter :: kev=1.6022d-16   !unit J
  real(p_),parameter :: elementary_charge=1.6022d-19   !unit C
  real(p_),parameter :: electron_mass=9.1094d-31 !in unit of kg
  real(p_),parameter :: epsilon0=8.8542d-12 !Permittivity of free space 
  real(p_),parameter :: pi=3.1415926_p_
  real(p_),parameter :: twopi=pi*2.0_p_
  real(p_),parameter :: fourpi=pi*4.0_p_
  real(p_),parameter :: half_pi=pi*0.5_p_
  real(p_),parameter :: mu0=fourpi*1.0d-7 !permeability in SI unit
  real(p_),parameter :: zero=0.0_p_
  real(p_),parameter :: one=1.0_p_
  real(p_),parameter :: two=2.0_p_
  real(p_),parameter :: three=3.0_p_
  real(p_),parameter :: four=4.0_p_
  real(p_),parameter :: five=5.0_p_
  real(p_),parameter :: six=6.0_p_
  real(p_),parameter :: seven=7.0_p_
  real(p_),parameter :: eight=8.0_p_
  real(p_),parameter :: nine=9.0_p_
  real(p_),parameter :: ten=10.0_p_
  real(p_),parameter :: eleven=11.0_p_
  real(p_),parameter :: twelve=12.0_p_
  real(p_),parameter :: thirteen=13.0_p_
  real(p_),parameter :: fourteen=14.0_p_
  real(p_),parameter :: fifteen=15.0_p_
  real(p_),parameter :: one_half=0.5_p_
  real(p_),parameter :: one_third=one/three
  real(p_),parameter :: one_fifth=0.2_p_
  real(p_),parameter :: three_halfs=1.5_p_
  real(p_),parameter :: c=299792458d0 !speed of light in vaccuum
end module constants


module normalizing
  use constants, only : p_
  implicit none
  save
  real(p_), parameter :: tu=2.14d0  !temperature unit in keV
  real(p_), parameter :: qu=1.6022d-19   !elementary_charge  in SI unit (Coulomb)
  real(p_), parameter :: nu=4.66d19 !SI unit (m^-3)
  real(p_), parameter :: Ln=1.0d0 !length unit (in unit of meter), please do not change the value because I have assumed this value (1.0) in some parts of the code
  real(p_), parameter :: bn=1.0d0 !nmagnetic field strength (in unit of Tesla), do not change the value
  real(p_) :: omegan_i, tn_i, vn_i, beta_ni !omegan_i is the ions gyro-angular-frequency in magnetic field bn
  real(p_) :: tn_main, dtao_main, vu
end module normalizing


module  poloidal_flux_2d !poloidal flux and its partial derivatives in (R,Z) plane
  use constants,only:p_
  implicit none
  save
  integer:: nx, nz !number of gridpoints in R and Z directions (specified in G-file)
  real(p_), dimension(:), allocatable :: xarray, zarray ! R and Z array
  real(p_), dimension(:,:), allocatable :: psi !psi is related to poloidal magnetic field by Bp=\nabla{psi}\times\nabla{phi}
  real(p_),dimension(:,:),allocatable :: psi_x, psi_z  !partial derivatives
  real(p_),dimension(:,:),allocatable :: psi_xx, psi_zz, psi_xz, psi_zx !partial derivatives
  real(p_),dimension(:,:),allocatable :: psi_gradient, y2a_gradient
  !strength of the gradient of the poloidal flux, y2a_gradient is the second derivative array used in 2D spline interpolation
  real(p_), dimension(:,:), allocatable :: y2a_psi, y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx
  ! y2a_* is the second derivative array used in 2D spline interpolation
end module poloidal_flux_2d


module radial_module
  use constants,only:p_
  implicit none
  save
  integer :: npsi
  real(p_) :: r_axis, z_axis, baxis, psi_axis, psi_lcfs
  real(p_) :: sign_bphi !sign of the toroidal component of the magnetic field
  real(p_),dimension(:),allocatable :: fpsi, ffprime, fprime, qpsi
  real(p_),dimension(:),allocatable :: pressure, pprime
  real(p_),dimension(:),allocatable :: psi_1d, pfn_npsi, tfn_npsi
  real(p_),dimension(:),allocatable :: q_with_sign
end module radial_module


module magnetic_coordinates
  use constants,only:p_
  implicit none
  save
  integer :: mpol,nflux, mpol2, nflux2 !mpol2 is number of poloidal gridpoints for perturbed field
  real(p_) :: pfn_inner, pfn_bdry, radial_width
  real(p_) :: dtheta, dradcor
  real(p_) :: radcor_low0, radcor_upp0
  real(p_) :: radcor_low2, radcor_upp2
  integer :: j_low2, j_upp2
  integer :: i_theta_zero !poloidal index corresponding to theta=0

  real(p_) :: GSpsi_prime !dGSpsi/dx, x is the normalized poloidal magnetic flux
  real(p_) :: abs_jacobian_min, abs_jacobian_max, vol
  real(p_) :: sign_of_jacobian, sign_of_GSpsi_prime

  real(p_),dimension(:,:), allocatable :: r_mc, z_mc !SI units
  real(p_),dimension(:,:), allocatable :: tor_shift_mc, jacobian, qhat
  real(p_),dimension(:), allocatable :: tor_shift_mc_left_bdry_minus_one 
  real(p_),dimension(:), allocatable :: av_jacobian
  !GSpsi_array is Grad-Shafranov poloidal flux in SI units, i.e. poiloidal_magnetic_flux/twopi
  real(p_),dimension(:), allocatable :: GSpsi_array, pfn, minor_r_array, minor_r_prime_array, &
       & radcor_1d_array, radcor_1d_array2, & !shranked radial array
       & theta_1d_array

  real(p_),dimension(:,:),allocatable :: dl_mc, &
       & rth, zth, rpsi, zpsi, &
       & grad_psi_r, grad_psi_z, &
       & grad_theta_r, grad_theta_z, &
       & grad_alpha_r, grad_alpha_z, grad_alpha_phi, &
       & grad_psi, grad_alpha, grad_theta,&
       & grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, grad_alpha_dot_grad_theta
  real(p_),dimension(:),allocatable :: grad_alpha_r_left_bdry_minus_one,grad_alpha_z_left_bdry_minus_one

  integer :: nsegment !computational region is 1/nsegment of the full torus
  integer ::  mtor !number of gridpoints along the toroidal direction
  real(p_) :: dtor !interval of toroidal gridpoints
  real(p_) :: toroidal_range !toroidal_range=twopi/nsegment
  real(p_),dimension(:), allocatable :: tor_1d_array !toroidal grids for perturbations
end  module magnetic_coordinates


module domain_decomposition
  use constants,only: p_
  implicit none
  integer :: numprocs, myid
  integer :: GRID_COMM, TUBE_COMM
  integer :: GCLR, TCLR, ntube, my_left, my_right, my_left2,my_right2
  integer :: GCLR_cut !the value of GCLR at theta cut
  integer :: GCLR_cut_left !the GCLR_id of the left neighbour of GCLR_cut
  integer :: multi_eq_cells
  !1D domain decomposion along theta coordinates, a processor treats [theta_start, theta_start+dtheta2]
  real(p_) :: dtheta2, theta_start 
  integer :: ipol_eq !index of theta_start in the original magnetic coordinate grids
end module domain_decomposition

module control_parameters
  use constants,only:p_
  use constants,only: elementary_charge,kev
  implicit none
  save
  integer  :: kstart,kend,niter
  real(p_) :: dtao_omega_i_axis
  character(100):: poloidal_angle_type
  integer  :: ion_spatial_loading_scheme  !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
  integer  :: ion_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
  integer  :: gk_spatial_loading_scheme !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
  integer  :: gk_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
  integer  :: iplot_mode_structure
  logical  :: remove_n0, filter_toroidal, filter_radial, diagnosis !if ture, output additional testing information
  integer  :: radial_harmonics_included, order_faraday_scheme
  integer :: space_charge_switch, fk_ions_switch, gk_species_switch
  integer :: fk_nonlinear, gk_nonlinear, ismooth
  integer :: nh !toroidal harmonics included are exp(i*nsegment*y) with i=(-nh, -nh+1,..., -1,0,1,2,...,nh-1,nh) , if remove_n0=.true. the 0 harmonic will be removed
end module control_parameters


module perturbation_field_matrix
  use constants,only:p_
  implicit none
  save

  real(p_),dimension(:,:,:),allocatable :: potential,phix, phiy, phiz !electrostatic potential and its derivatives
  real(p_),dimension(:,:,:),allocatable :: apara, ax, ay, az !parallel component of the vector potnetial and its derivative
  real(p_),dimension(:,:,:),allocatable :: apara_h, apara_s, apara_s_old !used in the mixed-variable pullback method
  real(p_),dimension(:,:,:),allocatable :: ahx, ahy, ahz

  real(p_),dimension(:,:),allocatable :: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
  real(p_),dimension(:,:),allocatable :: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right
  real(p_),dimension(:,:),allocatable :: ef_cyl_r_left_old,ef_cyl_z_left_old,ef_cyl_phi_left_old
  real(p_),dimension(:,:),allocatable :: ef_cyl_r_right_old,ef_cyl_z_right_old,ef_cyl_phi_right_old

  real(p_),dimension(:,:),allocatable :: my_den_i_left, my_den_i_right
  real(p_),dimension(:,:),allocatable :: my_den_e_left, my_den_e_right !density
  real(p_),dimension(:,:),allocatable :: my_jpar_e_left, my_jpar_e_right !parallel  current


contains

  subroutine allocate_field_matrix()
    use constants, only: p_, zero
    use magnetic_coordinates, only: m=>mtor, n=>nflux2
    use domain_decomposition, only: myid
    implicit none

    if(myid.eq.0) write(*,*) 'mtor, nflux2=', m, n

    allocate(my_den_i_left(m,n), my_den_i_right(m,n))
    allocate(my_den_e_left(m,n), my_den_e_right(m,n)) 
    allocate(my_jpar_e_left(m,n),my_jpar_e_right(m,n))

    allocate(potential(m,n, 2), source=zero)
    allocate(phix(m,n,2), source=zero) !dphi/dx
    allocate(phiy(m,n,2), source=zero) !dphi/dy
    allocate(phiz(m,n,2), source=zero) !dphi/dz
    allocate(apara  (m,n,2), source=zero)
    allocate(apara_h(m,n,2), source=zero)
    allocate(apara_s(m,n,2), source=zero)
    allocate(apara_s_old(m,n,2), source=zero)
    allocate(ax(m,n,2), source=zero) !dApar/dx
    allocate(ay(m,n,2), source=zero) !dApar/dy
    allocate(az(m,n,2), source=zero) !dApar/dz
    allocate(ahx(m,n,2), source=zero) !dApar_h/dx
    allocate(ahy(m,n,2), source=zero) !dApar_h/dy
    allocate(ahz(m,n,2), source=zero) !dApar_h/dz


    allocate(ef_cyl_r_left(m+1,n), source=zero)
    allocate(ef_cyl_z_left(m+1,n), source=zero)
    allocate(ef_cyl_phi_left(m+1,n), source=zero)
    allocate(ef_cyl_r_right(m+1,n), source=zero)
    allocate(ef_cyl_z_right(m+1,n), source=zero)
    allocate(ef_cyl_phi_right(m+1,n), source=zero)


    allocate(ef_cyl_r_left_old(m+1,n))
    allocate(ef_cyl_z_left_old(m+1,n))
    allocate(ef_cyl_phi_left_old(m+1,n))
    allocate(ef_cyl_r_right_old(m+1,n))
    allocate(ef_cyl_z_right_old(m+1,n))
    allocate(ef_cyl_phi_right_old(m+1,n))

  end subroutine allocate_field_matrix

end module perturbation_field_matrix

module flux_tube_model
  use constants,only:p_
  implicit none
  save
  integer:: j_fixed
  real(p_):: radcor_fixed !the radcor of the center of computational region
end module flux_tube_model

module boundary
  use constants,only:p_
  save
  integer:: nlim, np_lcfs
  real(p_),dimension(:),allocatable :: rlim, zlim, x_lcfs, z_lcfs
end module boundary


module mapping_module !from cylindrical coordinates to magnetic coordinates
  use constants,only:p_
  implicit none
  save
  integer,parameter:: nx_mapping=100,nz_mapping=100
  real(p_):: r_cyl(nx_mapping),z_cyl(nz_mapping)
  real(p_):: radcor(nx_mapping,nz_mapping)
  real(p_):: theta_a(nx_mapping,nz_mapping),theta_b(nx_mapping,nz_mapping)
  real(p_):: tor_shift_a(nx_mapping,nz_mapping),tor_shift_b(nx_mapping,nz_mapping)
  real(p_):: dr,dz
  integer:: i0,j0 !index of the point at magnetic axis
  real(p_):: dtheta_dr(nx_mapping,nz_mapping),dtheta_dz(nx_mapping,nz_mapping)
  real(p_):: ddelta_dr_a(nx_mapping,nz_mapping),ddelta_dz_a(nx_mapping,nz_mapping)
  real(p_):: ddelta_dr_b(nx_mapping,nz_mapping),ddelta_dz_b(nx_mapping,nz_mapping)
  real(p_):: dradial_dr(nx_mapping,nz_mapping),dradial_dz(nx_mapping,nz_mapping)
end module mapping_module
