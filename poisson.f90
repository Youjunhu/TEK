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
    use gk_polarization, only: prepare_polarization_matrix, prepare_polarization_matrix2, &
         & prepare_slowing_down_polarization_matrix
    use normalizing, only: tu,qu, nu
    use adiabatic_e_profiles, only: initialize_adiabatic_electron, ne_object, te_object
    implicit none
    complex(p_), allocatable :: polarization(:,:,:)
    complex(p_), allocatable :: mmm(:,:,:)
    real(p_) ::  x, tmp
    integer :: nx, j, jeq, info, kn, ns

    nx = nrad - 2
    allocate(polarization(nx,nx,0:nh_max),source=(zero,zero))
    allocate(mmm(nx,nx,0:nh_max))

    do ns = 1, nsm
       if(ns<30) then
          !call prepare_polarization_matrix(ns, mmm)
          call prepare_polarization_matrix2(ns, mmm)
       else
          call prepare_slowing_down_polarization_matrix(mmm, 4*atom_mass_unit, 2*elementary_charge, 3.5*Mev)
       endif
       polarization(:,:,:) = polarization + mmm
    enddo

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

 subroutine merge_source(my_jpar_left, my_jpar_right, jpar)
    use constants, only: p_
    use magnetic_coordinates, only: m=>mtor, n=>nrad
    use domain_decomposition, only: GRID_COMM, TUBE_COMM, GCLR, my_left, my_right
    use connection_condition
    use mpi
    implicit none
    real(p_), intent(in) :: my_jpar_left(m,n), my_jpar_right(m,n)
    real(p_), intent(out) :: jpar(m,n)
    real(p_) :: jpar_left(m,n), jpar_right(m,n), jpar_left0(m,n)
    integer :: status(MPI_STATUS_SIZE), ierr

    !summing over all those procs in a z-cell
    call MPI_ALLREDUCE(my_jpar_left,  jpar_left,  m*n, MPI_REAL8, MPI_SUM, GRID_COMM,ierr)
    call MPI_ALLREDUCE(my_jpar_right, jpar_right, m*n, MPI_REAL8, MPI_SUM, GRID_COMM,ierr)

    call MPI_Sendrecv(jpar_right, m*n, MPI_real8, my_right, 3, &
         &            jpar_left0, m*n, MPI_real8, my_left,  3, Tube_comm, status, ierr)

    if(GCLR==0) call connection_condition_at_theta_cut_for_deposition(jpar_left0) 
    jpar(:,:) = jpar_left(:,:) + jpar_left0(:,:) !add the contribution from the neighbour cell

  end subroutine merge_source

  subroutine solve_poisson(density_left, density_right, potential, phix, phiy, phiz)
    use constants, only: p_, zero, elementary_charge
    use magnetic_coordinates, only: m=>mtor, n=>nrad, av_jacobian, abs_jacobian, &
         & mpol2, xgrid, dradcor, dtor
    use control_parameters, only: fk_switch, filter_radial, dt_omega_i_axis, &
         & ismooth, nh_min, nh_max, adiabatic_electrons
    use gk_radial_profiles, only : density_object, temperature_object
    use normalizing,only: tu,qu, nu
    !use fk_module,only: mass_i,charge_i
    use gk_module, only : w_unit
    use filter_module, only: radial_sine_filter_core
    use domain_decomposition,only: myid, ipol_eq, tube_comm, multi_eq_cells, dtheta2, dvol
    use math, only : ZGETRS_wrapper
    use smoothing_module
    use derivatives_in_xyz
    use communication_connection
    use transform_module
    use mpi
    implicit none
    real(p_), intent(in) :: density_left(:,:), density_right(:,:)
    real(p_), intent(out), dimension(:,:,:) :: potential, phix, phiy, phiz
    real(p_) :: density(0:m-1,n), source(0:m-1,n-2), source_dst(0:m-1,n-2), phi(0:m-1,n-2)
    complex(p_) :: source_dft(0:m-1,n-2), phi_dft(0:m-1,n-2)
    real(p_) :: phi_n0(n-2), sum_phi(n-2), av_phi(n-2)
    integer :: kn, ierr, j,i, it, jeq
    real(p_) :: coeff, x

    call merge_source(density_left, density_right, density)
    do j=1,n
       density(:,j) = density(:,j)*w_unit/dvol(j)/nu
    enddo
    source = density(:,2:n-1)

    if(filter_radial .eqv. .true.) then !keeping only low-kr-harmonics
       call oned_sine_transform2(source,source_dst,m,n-2) !DST of s(:,:) along the 2nd dimension
       call radial_sine_filter_core(source_dst,m,n-2)
       call oned_inverse_sine_transform2(source_dst,source,m,n-2)
    endif

    !call oned_fourier_transform1(source, source_dft, m, n-2) !calculating DFT of source1(:,:) along the 1st dimension
    call oned_fourier_transform1_parallel_version(source, source_dft, m, n-2)

    if((adiabatic_electrons .eqv. .true.) .and. (nh_min==0)) then
       ! Solve for n=0 component while neglecting the electron contribution in the poisson_matrix. The solution is incoorrect,
       ! but the surface averaging of it equals to the surface average of the correct solution
       kn=0; call ZGETRS_wrapper(kn, poisson_matrix0, IPIV0, source_dft(kn,:), phi_dft(kn,:))
       ! Convert to real type and include the Jacobian (to perform magnetic surface average):
       phi_n0(:) = real(phi_dft(0,:)) * abs_jacobian(ipol_eq, 2:n-1)

       call MPI_Allreduce(phi_n0, sum_phi, n-2, MPI_Double, MPI_sum, tube_comm, ierr)
       av_phi(:) = (sum_phi(:)/mpol2)/av_jacobian(2:n-1) !magnetic surface averaging
       do j = 1, n-2 !add <phi> to the right-hand side
          x = xgrid(j+1)
          coeff = (density_object(1)%func(x)/nu)*(elementary_charge/qu)**2/(temperature_object(1)%func(x)/tu)
          source_dft(0,j) = source_dft(0,j) + av_phi(j)*coeff
       enddo
    endif

    phi_dft = 0 !only some toroidal harmonics will be solved, others are assumed to be zero
    ! do kn = nh_min, nh_max !solve the field equation for a single toroidal harmonic
    !    call ZGETRS_wrapper(kn, poisson_matrix, IPIV, source_dft(kn,:), phi_dft(kn,:)) 
    ! enddo
    ! do kn = nh_min, nh_max !for negative toroidal mode number, no need to solve, use the relation of complex conjogate 
    !    phi_dft(m-kn,:) = conjg(phi_dft(kn,:)) 
    ! enddo
    call solver_toroidal_mode_number_parallel(poisson_matrix, IPIV, source_dft, phi_dft)

    !call oned_backward_fourier_transform1(phi_dft, phi, m, n-2)
    call oned_backward_fourier_transform1_parallel_version(phi_dft, phi, m, n-2) !radial task decomposion

    potential(:,2:n-1,1) = phi(:,:)
    potential(:,1, 1) = 0._p_ !fixed zero boundary condition
    potential(:,n,1) = 0._p_  !fixed zero boundary condition

    call communicate_between_neighbour_cells(potential)
    do i=1, ismooth
       call smoothing_along_field_line_core5(potential(:,:,1))
    enddo
    call radial_derivative  (potential(:,:,1), phix(:,:,1))
    call toroidal_derivative(potential(:,:,1), phiy(:,:,1))
    call theta_derivative   (potential(:,:,:), phiz(:,:,1)) 
    call update_field_at_right_boundary_of_present_cell(phix, phiy, phiz)

    !call electric_potential_to_field() !this version uses epar, ex, and ey rather than the cylindrical components.
    if(fk_switch==1) call potential_to_cylindrical_field() !this version use cylindrical components Er,Ephi,Ez

  end subroutine solve_poisson
  


  subroutine solver_toroidal_mode_number_parallel(matrix, IPIV, source_dft, potential_dft)
    use constants, only : p_
    use magnetic_coordinates, only : mtor
    use control_parameters, only : nh_min, nh_max
    use domain_decomposition, only : TCLR, ntube, grid_comm
    use math, only : ZGETRS_wrapper
    use mpi
    implicit none
    complex(p_), intent(in) :: matrix(:,:,0:) !LU factorization of the radial coefficient matrix
    integer, intent(in) :: ipiv(:,0:) 
    complex(p_), intent(in) :: source_dft(0:,:)
    complex(p_), intent(inout) :: potential_dft(0:,:)
    logical, save :: is_first=.true.
    integer, save :: my_kn_start, my_kn_end, sendcount, n, nht
    integer, save, allocatable :: recvcounts(:), displace(:)
    complex(p_), allocatable  :: my_field(:,:), field(:,:)
    integer :: i, ierr, kn, ntask

    if (is_first.eqv..true.) then !needs to be genreated only for the first time and used multiple times
       is_first=.false.
       allocate(recvcounts(0:ntube-1))
       allocate(displace(0:ntube-1))
       n=size(source_dft,2) !radial grid number
       nht = nh_max - nh_min + 1 !total number of toroidal harmonics
       if(nht <= ntube) then
          my_kn_start = TCLR + nh_min
          my_kn_end = my_kn_start
          if(TCLR>nht-1) then
             my_kn_start = 0
             my_kn_end = -1
          endif
          recvcounts(0:nht-1)=1*n
          recvcounts(nht:)= 0
          do i=0,nht-1
             displace(i)=i*1*n
          enddo
       else
          ntask = nht/ntube
          my_kn_start = TCLR*ntask + nh_min
          my_kn_end = my_kn_start + ntask-1
          if(TCLR .eq. (ntube-1)) my_kn_end = nht-1 !the last process handles all the remainder part
          recvcounts(:)=ntask*n
          recvcounts(ntube-1)= (nht-ntask*(ntube-1))*n  !last process contains additional elements
          do i=0,ntube-1
             displace(i)=i*ntask*n
          enddo
       endif
       sendcount= (my_kn_end - my_kn_start +1)*n
    endif

    allocate(my_field(n, my_kn_start:my_kn_end))
    allocate(field(n, nh_min:nh_max))

    do kn = my_kn_start, my_kn_end !for non-negative toroidal mode number
       call ZGETRS_wrapper(kn, matrix, IPIV, source_dft(kn,:),my_field(:,kn)) 
    enddo

    call MPI_Allgatherv(my_field, sendcount, MPI_Double_Complex, &
         &      field, recvcounts, displace, MPI_Double_Complex, Grid_comm, ierr)

    do kn = nh_min, nh_max !re-arange the gathered data
       potential_dft(kn,:) = field(:,kn)
    enddo

    do kn = 1, nh_max !for negative toroidal mode number, no need to solve Poisson equation.
       potential_dft(mtor-kn,:) = conjg(potential_dft(kn,:)) !use the Complex conjugate relation
       !result is put in the usual order adopted by FFTW (and most FFT libraries) so that it can be backward transformed by FFTW,
    enddo
  end subroutine solver_toroidal_mode_number_parallel

end module poisson


subroutine potential_to_cylindrical_field()
  use mpi
  use constants,only: p_, one,two,one_half
  use magnetic_coordinates,only: mtor,nrad,dradcor,zgrid,dtheta,r_mc
  use magnetic_coordinates,only: xgrid,ygrid,&
       & grad_psi_r,grad_psi_z,grad_theta_r,grad_theta_z, &
       & grad_alpha_r,grad_alpha_z,grad_alpha_phi

  use perturbation_field,only: potential
  use domain_decomposition,only: theta_start,ipol_eq
  use perturbation_field,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left !as output
  use control_parameters,only: filter_toroidal,filter_radial
  use constants,only:kev
  use normalizing,only:bn,ln, tu,qu
  use fk_module, only : vn_fk
  use domain_decomposition,only:myid
  !use delta_ne_mod,only: delta_ne, delta_ne_theta,delta_ne_psi,delta_ne_alpha !for testing
  use transform_module,only:oned_fourier_transform1, oned_backward_fourier_transform1,&
       &oned_sine_transform2,oned_inverse_sine_transform2
  use filter_module,only: toroidal_filter,radial_sine_filter_core
  use derivatives_in_xyz,only: radial_derivative,toroidal_derivative,theta_derivative
  use communication_connection
  implicit none

  real(p_):: potential_psi(mtor,nrad),potential_alpha(mtor,nrad),potential_theta(mtor,nrad)
  real(p_):: ef_cyl_r(mtor,nrad),ef_cyl_z(mtor,nrad),ef_cyl_phi(mtor,nrad)
  !complex(p_):: ef_cyl_r_dft(mtor,nrad),ef_cyl_z_dft(mtor,nrad),ef_cyl_phi_dft(mtor,nrad)
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

  call radial_derivative  (potential(:,:,1),potential_psi)
  call toroidal_derivative(potential(:,:,1),potential_alpha)
  call theta_derivative(potential,potential_theta) !derivative along the magnetic field line
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


subroutine communicate_field_value_between_neighbour_cells2() !use cylindrical components of the electric field
  use mpi
  use constants,only:p_
  use perturbation_field,only: ef_cyl_r_left, ef_cyl_z_left, ef_cyl_phi_left !already known before entering this subroutine
  use perturbation_field,only: ef_cyl_r_right, ef_cyl_z_right, ef_cyl_phi_right !as output
  use magnetic_coordinates,only: m=>mtor,n=>nrad
  use domain_decomposition,only:myid,GCLR,Tube_comm,ntube,GCLR_cut,GCLR_cut_left,my_left,my_right
  use connection_condition,only:  connection_condition_at_theta_cut
  implicit none

  integer:: status(MPI_STATUS_SIZE),ierr
  !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc, the field on the right-boundary is received from the neighbour proc. Note that the definition of cell here is different from the the definition of cell in PIC: the grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.

  call MPI_Sendrecv(ef_cyl_r_left,  (m+1)*n,  MPI_real8, my_left,  4,&
       &            ef_cyl_r_right, (m+1)*n,  MPI_real8, my_right, 4,Tube_COMM,status,ierr)
  call MPI_Sendrecv(ef_cyl_z_left,  (m+1)*n,  MPI_real8, my_left,  5,&
       &            ef_cyl_z_right, (m+1)*n,  MPI_real8, my_right, 5,Tube_COMM,status,ierr)
  call MPI_Sendrecv(ef_cyl_phi_left, (m+1)*n,  MPI_real8, my_left,  6,&
       &            ef_cyl_phi_right,(m+1)*n,  MPI_real8, my_right, 6,Tube_COMM,status,ierr)

  !if(GCLR.eq.GCLR_cut) then !special treatment at theta cut
  if(GCLR.eq.GCLR_cut_left) then !special treatment at theta cut
     call connection_condition_at_theta_cut(ef_cyl_r_right) 
     call connection_condition_at_theta_cut(ef_cyl_z_right) 
     call connection_condition_at_theta_cut(ef_cyl_phi_right) 
  endif

end subroutine communicate_field_value_between_neighbour_cells2
