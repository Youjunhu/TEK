module poisson
  use constants, only : p_
  save
  complex(p_), allocatable :: poisson_matrix(:,:,:) 
  complex(p_), allocatable :: poisson_matrix0(:,:,:) !pure gk contribution, used in solving n=0 component when using adiabatic electrons, 
  integer, allocatable :: ipiv(:,:), ipiv0(:,:)

contains
  subroutine prepare_poisson_matrix()
    use constants, only: zero, elementary_charge, atom_mass_unit, Mev
    use magnetic_coordinates, only: nrad, xgrid, mtor
    use gk_radial_profiles, only: density_object, temperature_object
    use control_parameters, only: dt_omega_i_axis, nh_max, adiabatic_electrons
    use gk_module, only: nsm
    use gk_polarization
    use normalizing, only: tu,qu, nu
    use domain_decomposition, only: myid, numprocs, tclr, gclr
    use adiabatic_e_profiles, only: initialize_adiabatic_electron, ne_object, te_object
    implicit none
    complex(p_), allocatable :: polarization(:,:,:)
    complex(p_), allocatable :: mmm(:,:,:)
    real(p_) ::  x, tmp
    integer :: nx, i, j, jeq, info, kn, ns, u

    nx = nrad - 2
    allocate(polarization(nx, nx, 0:nh_max), source=(zero,zero))
    allocate(mmm(nx, nx, 0:nh_max))

    do ns = 1, nsm
       !if(ns<30) then
          call prepare_polarization_matrix(ns, mmm)
          !call prepare_polarization_matrix2(ns, mmm, nx)
          !call prepare_polarization_matrix20(ns, mmm, nx)
          !call prepare_polarization_matrix3(ns, mmm)
       !else
       !   call prepare_slowing_down_polarization_matrix(mmm, 4*atom_mass_unit, 2*elementary_charge, 3.5*Mev)
       !endif
       polarization(:,:,:) = polarization + mmm
    enddo
if(tclr==0) write(*,*) gclr, real(polarization(nx/2,nx/2, nh_max)), imag(polarization(nx/2,nx/2, nh_max))

    allocate(poisson_matrix(nx, nx, 0:nh_max))
    allocate(ipiv(nx, 0:nh_max))

    poisson_matrix = polarization

    if(adiabatic_electrons .eqv. .true.) then
       allocate(poisson_matrix0(nx, nx, 0:0))
       allocate(ipiv0(nx,0:0))
       poisson_matrix0 = polarization(:,:, 0:0)
       call initialize_adiabatic_electron()
       do j = 1, nx 
          jeq = j+1
          x = xgrid(jeq)
          tmp = (ne_object%func(x)/nu)*(elementary_charge/qu)**2/(te_object%func(x)/tu)
          poisson_matrix(j, j, :) = poisson_matrix(j, j, :) + tmp !add to the diagonal elemets of the matrix
       enddo
       kn = 0
       call ZGETRF(nx, nx, poisson_matrix0(:,:,kn), nx, ipiv0(:,kn), info) 
    endif

    do kn = 0, nh_max !for each toroidal Fourier component, LU factorize the radial coeff matrix
       ! Computes LU factorization of a general matrix using partial pivoting with row interchanges.
       call ZGETRF(nx, nx, poisson_matrix(:,:,kn), nx, ipiv(:,kn), info) 
    enddo

  end subroutine prepare_poisson_matrix


  subroutine solve_poisson(density_left, density_right, potential, phix, phiy, phiz, phi_dft)
    use constants, only: p_, zero, elementary_charge
    use magnetic_coordinates, only: m=>mtor, n=>nrad, jacobian_av, abs_jacobian, &
         & mpol2, xgrid, dradcor, dtor
    use control_parameters, only: fk_switch, filter_radial, dt_omega_i_axis, &
         & ismooth, nh_min, nh_max, adiabatic_electrons
    use gk_radial_profiles, only : density_object, temperature_object
    use normalizing,only: tu,qu, nu
    !use fk_module,only: mass_i,charge_i
    use gk_module, only : w_unit
    use domain_decomposition,only: myid, ipol_eq, tube_comm, multi_eq_cells, dtheta2, dvol
    use math, only : ZGETRS_wrapper, solver_toroidal_mode_number_parallel
    use filter_module
    use smoothing_module
    use derivatives_in_xyz
    use communication_connection
    use transform_module
    use mpi
    implicit none
    real(p_), intent(in) :: density_left(:,:), density_right(:,:)
    real(p_), intent(out), dimension(:,:,:) :: potential, phix, phiy, phiz
    complex(p_), intent(out) :: phi_dft(0:m-1, n-2)
    real(p_) :: density(0:m-1,n), source(0:m-1,n-2), source_dst(0:m-1,n-2), phi(0:m-1,n-2)
    complex(p_) :: source_dft(0:m-1, n-2), phix_dft(0:m-1, n-2)
    real(p_) :: signal_dst(1:1, n-2), signal(1:1, n-2)
    real(p_) :: phi_n0(n-2), sum_phi(n-2), av_phi(n-2)
    integer :: kn, ierr, j,i, it, jeq
    real(p_) :: coeff, x

    call merge_source(density_left, density_right, density)
    do j = 1, n
       density(:,j) = density(:,j)*w_unit/dvol(j)/nu
    enddo
    source = density(:,2:n-1)
    
    call oned_DFT_parallel_version(source, source_dft, m, n-2)

    if(nh_min == 0) then !remove harmonics of low kx to suppress numerical instabilities
       call oned_sine_transform2(real(source_dft(0:0,:)), signal_dst, 1, n-2)
       call radial_sine_high_pass(signal_dst, 1, n-2)
       call oned_inverse_sine_transform2(signal_dst, signal, 1, n-2)
       source_dft(0,:) = signal(1,:)
    endif
    
    phi_dft = 0 !only some toroidal harmonics will be solved, others are assumed to be zero

    if((adiabatic_electrons .eqv. .true.) .and. (nh_min==0)) then
       ! Solve for n=0 component while neglecting the electron contribution in the poisson_matrix.
       !The solution is incoorrect, but the surface averaging of it equals to the surface average of the correct solution
       kn=0; call ZGETRS_wrapper(kn, poisson_matrix0, IPIV0, source_dft(kn,:), phi_dft(kn,:))
       ! Convert to real type and include the Jacobian (to perform magnetic surface average):
       phi_n0(:) = real(phi_dft(0,:)) * abs_jacobian(ipol_eq, 2:n-1)

       call MPI_Allreduce(phi_n0, sum_phi, n-2, MPI_Double, MPI_sum, tube_comm, ierr)
       av_phi(:) = (sum_phi(:)/mpol2)/jacobian_av(2:n-1) !magnetic surface averaging
       do j = 1, n-2 !add <phi> to the right-hand side
          x = xgrid(j+1)
          coeff = (density_object(1)%func(x)/nu)*(elementary_charge/qu)**2/(temperature_object(1)%func(x)/tu)
          source_dft(0,j) = source_dft(0,j) + av_phi(j)*coeff
       enddo
    endif

    call solver_toroidal_mode_number_parallel(poisson_matrix, IPIV, source_dft, phi_dft)
    do kn = nh_min, nh_max
       !call mfilter_for_each_n(phi_dft(kn,:), n-2, kn)
    enddo
    do kn = 1, nh_max ! For negative toroidal mode number
       phi_dft(m-kn,:) = conjg(phi_dft(kn,:)) !use the Complex conjugate relation
    enddo
    !if(nh_min == 0) call mfilter_for_n0(phi_dft(0,:), n-2)
    if(nh_min == 0) call surface_average_of_n0(phi_dft(0,:), n-2)

    call oned_backward_DFT_parallel_version(phi_dft, phi, m, n-2) !radial task decomposion
    
    call x_derivative0(phi_dft, phix_dft) !more accurate than taking x derivative in real space
    call oned_backward_DFT_parallel_version(phix_dft, phix(:, 2:n-1, 1), m, n-2) 
    
    potential(:, 2:n-1, 1) = phi(:,:)
    potential(:, 1, 1) = 0._p_ !zero boundary condition
    potential(:, n, 1) = 0._p_  !zero boundary condition
    call update_scalar_at_right_boundary(potential)

    do j = 1, ismooth
       call smoothing_along_field_line_core(potential(:,:,1))
    enddo

    !call x_derivative(potential(:,:,1), phix(:,:,1))
    phix(:, 1, 1) = 0; phix(:, n, 1) = 0
    call y_derivative(potential(:,:,1), phiy(:,:,1))
    call z_derivative(potential(:,:,:), phiz(:,:,1)) 
    call update_derivatives_at_right_boundary(phix, phiy, phiz)

    !call electric_potential_to_field() !this version uses epar, ex, and ey rather than the cylindrical components.
    if(fk_switch==1) call potential_to_cylindrical_field() !this version use cylindrical components Er,Ephi,Ez

  end subroutine solve_poisson

end module poisson

subroutine potential_to_cylindrical_field()
  use mpi
  use constants,only: p_, one,two,one_half
  use magnetic_coordinates,only: mtor,nrad,dradcor,zgrid,dtheta,r_mc
  use magnetic_coordinates,only: grad_psi_r,grad_psi_z,grad_theta_r,grad_theta_z, &
       & grad_alpha_r,grad_alpha_z,grad_alpha_phi

  use perturbation_field,only: potential
  use domain_decomposition,only: theta_start,ipol_eq
  use perturbation_field,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left !as output
  use control_parameters,only: filter_radial
  use constants,only:kev
  use normalizing,only:bn,ln, tu,qu
  use fk_module, only : vn_fk
  use domain_decomposition,only:myid
  use transform_module,only:oned_fourier_transform1, oned_backward_fourier_transform1,&
       &oned_sine_transform2,oned_inverse_sine_transform2
  use filter_module,only: toroidal_filter,radial_sine_filter_core
  use derivatives_in_xyz,only: x_derivative,y_derivative,z_derivative, z_derivative0
  use communication_connection, only: communicate_field_value_between_neighbour_cells2
  implicit none

  real(p_):: potential_psi(mtor,nrad),potential_alpha(mtor,nrad),potential_theta(mtor,nrad)
  real(p_):: ef_cyl_r(mtor,nrad),ef_cyl_z(mtor,nrad),ef_cyl_phi(mtor,nrad)
  complex(p_):: potential_dft(mtor,nrad),out(mtor,nrad)
  real(p_):: potential_dst(mtor,nrad) !discrete sine transform
  integer:: i,j,i_left,i_right,j_left,j_right,ierr
  integer:: jeq
  real(p_):: normal
  real(p_):: tmp1(mtor,nrad),tmp2(mtor,nrad),tmp3(mtor,nrad)


  !smoothing is moved outside of this subroutine
!!$  if(filter_toroidal.eqv..true.) then !filter over the toroidal mode number, keeping the perturbation with desired toroidal mode number
!!$     call oned_fourier_transform1(potential,potential_dft,mtor,nrad) !calculating 1d DFT of s(:,:) along the first dimension
!!$     call toroidal_filter(potential_dft,mtor,nrad)
!!$     call oned_backward_fourier_transform1(potential_dft,potential,mtor,nrad)
!!$  endif

!!$  if(filter_radial.eqv..true.) then !filter over the radial mode number, keeping only low-radial-harmonics of the perturbation
!!$     call oned_sine_transform2(potential,potential_dst,mtor,nrad) !calculating 1d DST of s(:,:) along the second dimension
!!$     call radial_sine_filter_core(potential_dst,mtor,nrad)
!!$     call oned_inverse_sine_transform2(potential_dst,potential,mtor,nrad)
!!$  endif

  call x_derivative  (potential(:,:,1),potential_psi)
  call y_derivative(potential(:,:,1),potential_alpha)
  call z_derivative(potential,potential_theta) !derivative along the magnetic field line
  !write(*,*) potential_theta(10,40),tmp1(10,40),'myid=',myid
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do i=1,mtor
     do j=1,nrad
        jeq=j
        ef_cyl_r(i,j)=-potential_psi(i,j)*grad_psi_r(ipol_eq,jeq)&
             & -potential_theta(i,j)*grad_theta_r(ipol_eq,jeq)-potential_alpha(i,j)*grad_alpha_r(ipol_eq,jeq)

        ef_cyl_z(i,j)=-potential_psi(i,j)*grad_psi_z(ipol_eq,jeq) &
             & -potential_theta(i,j)*grad_theta_z(ipol_eq,jeq)-potential_alpha(i,j)*grad_alpha_z(ipol_eq,jeq)
        !ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mc(ipol_eq,jeq)*grad_alpha_phi(ipol_eq,jeq) !wrong, an additional 1/R factor is wrongly included, a bug found on 2018-Jan.-4 evening
        ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mc(ipol_eq,jeq) !corrected
     enddo
  enddo


  normal=tu*kev/(ln*bn*vn_fk*qu) !;write(*,*) 'normal=',normal
  !write(*,*) 'normal=',normal
  ef_cyl_r=ef_cyl_r*normal !normalized by the unit used in the code (Bn*vn_fk)
  ef_cyl_z=ef_cyl_z*normal
  ef_cyl_phi=ef_cyl_phi*normal

  do j=1,nrad
     do i=1,mtor
        ef_cyl_r_left(i,j)=ef_cyl_r(i,j) !store the data in the proper arrays
        ef_cyl_z_left(i,j)=ef_cyl_z(i,j)
        ef_cyl_phi_left(i,j)=ef_cyl_phi(i,j)
     enddo
     ef_cyl_r_left(mtor+1,j)=ef_cyl_r_left(1,j)  !peroidic toroidal boundary condition
     ef_cyl_z_left(mtor+1,j)=ef_cyl_z_left(1,j) 
     ef_cyl_phi_left(mtor+1,j)=ef_cyl_phi_left(1,j)
  enddo

  ef_cyl_r_left(:,1)=0._p_
  ef_cyl_z_left(:,1)=0._p_
  ef_cyl_phi_left(:,1)=0._p_

  ef_cyl_r_left(:,nrad)=  0._p_
  ef_cyl_z_left(:,nrad)=  0._p_
  ef_cyl_phi_left(:,nrad)= 0._p_

  call communicate_field_value_between_neighbour_cells2() !use ER,Ephi, EZ
end subroutine potential_to_cylindrical_field



