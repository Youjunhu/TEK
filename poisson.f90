module poisson
  use constants, only : p_
  save
  complex(p_),allocatable:: poisson_matrix(:,:,:)
  integer, allocatable:: ipiv(:,:)
  real(p_),allocatable :: adiabatic(:)
  real(p_),parameter :: damping_coeff= 0.0d0 !must be zero or positive, does not work 
  real(p_),parameter :: fraction_adiabatic=0.0_p_
contains
  subroutine prepare_poisson_matrix()
    use constants, only : elementary_charge
    use magnetic_coordinates,only: nflux2
    use magnetic_coordinates,only: j_low2, radcor_1d_array,mtor
    use density_temperature_profile_mod, only : ti_func, ni_func
    use fk_module,only: mass_i,charge_i 
    use control_parameters,only: dtao_omega_i_axis, nh
    use gk_polarization,only: polarization, prepare_polarization_matrix, prepare_polarization_matrix2
    use normalizing, only : tu,qu, nu
    implicit none
    real(p_) :: s, x
    integer :: nx, j, jeq, info, kn

    !call prepare_polarization_matrix()
    call prepare_polarization_matrix2()
    nx=nflux2-2
    allocate(poisson_matrix(nx,nx,0:nh))
    allocate(ipiv(nx,0:nh))
    allocate(adiabatic(nx))
    poisson_matrix = polarization

    do j=1,nx  !add artifical damping term for high-frequency modes
       poisson_matrix(j,j,:)=poisson_matrix(j,j,:) + damping_coeff/dtao_omega_i_axis
    enddo
    
    goto 100
    do j=1,nx !for testing correctness of polarization matrix by using the simple adiabatic electron model
       jeq=j+j_low2
       x=radcor_1d_array(jeq)
       adiabatic(j)=(ni_func(x)/nu)*(elementary_charge/qu)**2/(ti_func(x)/tu)
       poisson_matrix(j,j, :) = poisson_matrix(j,j, :) + adiabatic(j) !add the contribution from adiabatic species to the diagonal elemets of the matrix
    enddo

100 do kn=0,nh !for each toroidal Fourier component, LU factorize the radial coeff matrix
       call ZGETRF(nx,nx,poisson_matrix(:,:,kn),nx,ipiv(:,kn),info) !ZGETRF computes LU factorization of a general matrix using partial pivoting with row interchanges.
    enddo

  end subroutine prepare_poisson_matrix

  subroutine solve_poisson()
    use constants,only:p_, zero
    use magnetic_coordinates,only: m=>mtor, n=>nflux2, jacobian, av_jacobian, j_low2, j_upp2, mpol2, radcor_1d_array
    use control_parameters,only: fk_ions_switch, remove_n0, filter_radial, dtao_omega_i_axis, ismooth, nh
    use perturbation_field_matrix,only: potential, phix, phiy, phiz  !output
    use density_temperature_profile_mod, only : ti_func, ni_func
    use normalizing,only: tu,qu, nu
    use fk_module,only: mass_i,charge_i 
    use transform_module
    use filter_module, only: radial_sine_filter_core
    use domain_decomposition,only: myid, ipol_eq, tube_comm
    use smoothing_module
    use derivatives_in_xyz
    use communication_connection
    use math, only : ZGETRS_wrapper
    use mpi
    implicit none
    real(p_):: source(0:m-1,n-2), source_dst(0:m-1,n-2), potential0(0:m-1,n-2)
    complex(p_):: source_dft(0:m-1,n-2), potential_dft(0:m-1,n-2)
    real(p_) :: potential_n0(n-2), sum_potential(n-2), av_potential(n-2)
    integer:: kn, ierr, j,i, it
    real(p_):: coeff, radcor

    call merge_source_term(source)
!!$    do j=1,n-2
!!$       source(:,j) = source(:,j) + potential_old(:,j+1,1)*damping_coeff/dtao_omega_i_axis !adding the artificial damping term                           
!!$    enddo

    if(filter_radial.eqv..true.) then !filter over the radial mode number, keeping only low-radial-harmonics of the perturbation
       call oned_sine_transform2(source,source_dst,m,n-2) !calculating 1d DST of s(:,:) along the second dimension
       call radial_sine_filter_core(source_dst,m,n-2)
       call oned_inverse_sine_transform2(source_dst,source,m,n-2)
    endif
    !call oned_fourier_transform1_parallel_version(source,source_dft,m,n-2) !calculating DFT of source1(:,:) along the first dimension
    call oned_fourier_transform1(source,source_dft,m,n-2) !calculating DFT of source1(:,:) along the first dimension
    potential_dft=0._p_ !initialized to zero
    !solve for n=0 component (which is incoorrect) but the surface averaging of the solution equals to the surface average of the correct solution
!!$    kn=0; call ZGETRS_wrapper(kn, poisson_matrix0, IPIV0, source_dft(kn,:),potential_dft(kn,:)) !solve the field equation for n=0 toroidal harmonic
!!$    potential_n0=potential_dft(0,:)*abs(jacobian(ipol_eq, j_low2+1:j_upp2-1)) !convert to real and include the Jacobian factor
!!$    call MPI_Allreduce(potential_n0,  sum_potential, n-2, MPI_Double,MPI_sum,  tube_comm, ierr)
!!$    av_potential(:)=(sum_potential(:)/mpol2)/av_jacobian(j_low2+1:j_upp2-1) !magnetic surface averaging
!!$    do j=1,n-2 !add magnetic surface average of potential to the right-hand side
!!$       radcor=radcor_1d_array(j+j_low2)
!!$       coeff=(ni_func(radcor)/nu)*(charge_i/qu)**2/(ti_func(radcor)/tu)
!!$       source_dft(0,j) = source_dft(0,j) + av_potential(j)*coeff
!!$    enddo

!!$    call solver_toroidal_mode_number_parallel(poisson_matrix, IPIV, source_dft, potential_dft)
!!$
    do kn=0, nh !for non-negative toroidal mode number
       call ZGETRS_wrapper(kn, poisson_matrix, IPIV, source_dft(kn,:),potential_dft(kn,:)) !solve the field equation for a single toroidal harmonic
    enddo
    do kn=1,nh !for negative toroidal mode number, no need to solve Poisson equation because of the relation of complex conjogate to the positive toroidal mode number
       potential_dft(m-kn,:)=conjg(potential_dft(kn,:)) 
    enddo

    if(remove_n0 .eqv. .true.) potential_dft(0,:)=zero
    !call oned_backward_fourier_transform1_parallel_version(potential_dft, potential0, m, n-2) !radial task decomposion
    call oned_backward_fourier_transform1(potential_dft, potential0, m, n-2) 
    potential(:,2:n-1,1)= potential0(:,:)
    potential(:,1, 1)=0._p_ !fixed zero boundary condition
    potential(:,n,1)=0._p_ !fixed zero boundary condition

    call communicate_field_value_between_neighbour_cells(potential)
    do i=1, ismooth
       call smoothing_along_field_line_core5(potential(:,:,1))
    enddo
    call radial_derivative  (potential(:,:,1),phix(:,:,1))
    call toroidal_derivative(potential(:,:,1),phiy(:,:,1))
    call theta_derivative   (potential(:,:,:),phiz(:,:,1)) 
    call update_field_at_right_boundary_of_present_cell(phix,phiy,phiz)
    !call electric_potential_to_field() !this version uses epar, ex, and ey rather than the cylindrical components. ex, ey, epar are needed in gk electron pusher
    if(fk_ions_switch==1) call potential_to_cylindrical_field() !this version use cylindrical components Er,Ephi,Ez

  end subroutine solve_poisson

  subroutine solver_toroidal_mode_number_parallel(matrix, IPIV, source_dft, potential_dft)
    use constants, only : p_
    use magnetic_coordinates, only : mtor
      use control_parameters, only : nh
    use domain_decomposition, only : TCLR, ntube, grid_comm
    use math, only : ZGETRS_wrapper
    use mpi
    implicit none
    complex(p_), intent(in) :: matrix(:,:,0:) !LU factorization of the radial coefficient matrix
    integer, intent(in) :: ipiv(:,0:) 
    complex(p_),intent(in):: source_dft(0:,:)
    complex(p_),intent(inout):: potential_dft(0:,:)
    logical, save :: is_first=.true.
    integer, save :: my_kn_start, my_kn_end, sendcount, n
    integer, save, allocatable :: recvcounts(:), displace(:)
    complex(p_), allocatable  :: my_field(:,:), field(:,:)
    integer :: i, ierr, kn, ntask

    if (is_first.eqv..true.) then !needs to be genreated only for the first time and used multiple times
       is_first=.false.
       allocate(recvcounts(0:ntube-1))
       allocate(displace(0:ntube-1))
       n=size(source_dft,2) !radial grid number
       if((nh+1)<ntube) then
          my_kn_start=TCLR
          my_kn_end=TCLR
          if(TCLR>nh) then
             my_kn_start=0
             my_kn_end=-1
          endif
          recvcounts(0:nh)=1*n
          recvcounts(nh+1:)= 0
          do i=0,nh
             displace(i)=i*1*n
          enddo
       else
          ntask=(1+nh)/ntube
          my_kn_start=TCLR*ntask
          my_kn_end=my_kn_start+ntask-1
          if(TCLR.eq.(ntube-1)) my_kn_end=nh !the last process handles all the remainder part, in the case that nh+1 is not a perfect multiple of ntube
          recvcounts(:)=ntask*n
          recvcounts(ntube-1)= (nh+1-ntask*(ntube-1))*n  !last process contains additional elements
          do i=0,ntube-1
             displace(i)=i*ntask*n
          enddo
       endif
       sendcount= (my_kn_end - my_kn_start +1)*n
    endif
    allocate(my_field(n, my_kn_start:my_kn_end))
    allocate(field(n, 0:nh))

    do kn=my_kn_start, my_kn_end !for non-negative toroidal mode number
       call ZGETRS_wrapper(kn, matrix, IPIV, source_dft(kn,:),my_field(:,kn)) !solve the field equation for a single toroidal harmonic
    enddo
    call MPI_Allgatherv(my_field, sendcount, MPI_Double_Complex, &
         &    field, recvcounts, displace, MPI_Double_Complex, Grid_comm, ierr)
    do kn=0,nh !re-arange the gathered data
       potential_dft(kn,:)=field(:,kn)
    enddo
    do kn=1,nh !for negative toroidal mode number, no need to solve Poisson equation because of the relation of complex conjogate to the positive toroidal mode number
       potential_dft(mtor-kn,:)=conjg(potential_dft(kn,:)) !corresponding to the negative toroidal mode number, result is put in the usual order adopted by FFTW (and most FFT libraries) so that it can be backward transformed by FFTW,
    enddo
  end subroutine solver_toroidal_mode_number_parallel

  subroutine merge_source_term(source)
    use constants,only: p_
    use magnetic_coordinates,only:m=>mtor,n=>nflux2, j_low2, dtor,dradcor, jacobian, radcor_1d_array
    use control_parameters,only: fk_ions_switch, gk_species_switch
    use gk_module, only : w_unit
    use normalizing, only : tu,qu, nu
    use fk_module,only: charge_i, ni0,ti0
    use density_temperature_profile_mod, only : ti_func, ni_func
    use perturbation_field_matrix,only : my_den_i_left, my_den_i_right,&
         &   my_den_e_left, my_den_e_right
    use domain_decomposition,only: GRID_COMM,TUBE_COMM, ntube, GCLR, GCLR_cut,my_left,my_right,&
         &  dtheta2,ipol_eq, multi_eq_cells
    use connection_condition
    use mpi
    implicit none
    real(p_), intent(out) :: source(:,:)
    real(p_), dimension(m,n) :: my_source_left, my_source_right
    real(p_), dimension(m,n) :: source_left, source_right, source_left_tmp
    integer:: status(MPI_STATUS_SIZE),ierr, j, jeq
    real(p_) :: dv1, dv2

    do j=1,n  !divided by the spatial volume of a cell, in which a grid is the center.
       jeq=j_low2+(j-1)
       dv1=abs(jacobian(ipol_eq,jeq))*dradcor*dtheta2*dtor !volume of the cell (the center of the cell is the grid)
       dv2=abs(jacobian(ipol_eq+multi_eq_cells,jeq))*dradcor*dtheta2*dtor !volume of the cell 
       my_den_e_left (:,j)  = my_den_e_left (:,j)*(w_unit/dv1/nu)
       my_den_e_right(:,j)  = my_den_e_right(:,j)*(w_unit/dv2/nu)
    enddo

    if(fk_ions_switch==1) then  !include contribution from fk species
       my_source_left (:,:) = my_den_e_left (:,:)+ (charge_i/qu)*my_den_i_left (:,:) 
       my_source_right(:,:) = my_den_e_right(:,:)+ (charge_i/qu)*my_den_i_right(:,:)
    else
       my_source_left (:,:) = my_den_e_left (:,:)
       my_source_right(:,:) = my_den_e_right(:,:)
    endif

    call MPI_ALLREDUCE(my_source_left,   source_left,   m*n,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
    call MPI_ALLREDUCE(my_source_right,  source_right,  m*n,MPI_REAL8,MPI_SUM,GRID_COMM,ierr)
    call MPI_Sendrecv(source_right,   m*n, MPI_real8, my_right, 3,&
         &            source_left_tmp,m*n, MPI_real8, my_left,  3, Tube_comm, status, ierr)
    if(GCLR.eq.GCLR_cut) then !special treatment at the theta cut
       call connection_condition_at_theta_cut_for_deposition(source_left_tmp) 
    endif
    source_left(:,:) = source_left(:,:) + source_left_tmp(:,:) !add the contribution from the neighbour cell
    source(:,:) = source_left(:,2:n-1) !2:n-1 rather than 1:n, is to remove two boundary points, sine we do not solve the field equation at these two points
  end subroutine merge_source_term
  
end module poisson


subroutine potential_to_cylindrical_field()
  use mpi
  use constants,only: one,two,one_half
  use constants,only:p_
  use magnetic_coordinates,only: mtor,nflux2,dradcor,j_low2,theta_1d_array,dtheta,r_mc
  use magnetic_coordinates,only: radcor_1d_array2,tor_1d_array,&
       & grad_psi_r,grad_psi_z,grad_theta_r,grad_theta_z, &
       & grad_alpha_r,grad_alpha_z,grad_alpha_phi

  use perturbation_field_matrix,only: potential
  use domain_decomposition,only: theta_start,ipol_eq
  use perturbation_field_matrix,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left !as output
  use control_parameters,only: filter_toroidal,filter_radial
  use constants,only:kev
  use normalizing,only:bn,ln,vn_i, tu,qu
  use domain_decomposition,only:myid
  use delta_ne_mod,only: delta_ne, delta_ne_theta,delta_ne_psi,delta_ne_alpha !for testing
  use transform_module,only:oned_fourier_transform1, oned_backward_fourier_transform1,&
       &oned_sine_transform2,oned_inverse_sine_transform2
  use filter_module,only: toroidal_filter,radial_sine_filter_core
  use derivatives_in_xyz,only: radial_derivative,toroidal_derivative,theta_derivative
  use communication_connection
  implicit none

  real(p_):: potential_psi(mtor,nflux2),potential_alpha(mtor,nflux2),potential_theta(mtor,nflux2)
  real(p_):: ef_cyl_r(mtor,nflux2),ef_cyl_z(mtor,nflux2),ef_cyl_phi(mtor,nflux2)
  !complex(p_):: ef_cyl_r_dft(mtor,nflux2),ef_cyl_z_dft(mtor,nflux2),ef_cyl_phi_dft(mtor,nflux2)
  complex(p_):: potential_dft(mtor,nflux2),out(mtor,nflux2)
  real(p_):: potential_dst(mtor,nflux2) !discrete sine transform
  integer:: i,j,i_left,i_right,j_left,j_right,ierr
  integer:: jeq
  real(p_):: normal
  real(p_):: tmp1(mtor,nflux2),tmp2(mtor,nflux2),tmp3(mtor,nflux2)

  !----testing
!!$  do i=1,mtor
!!$     do j=1,nflux2
!!$        den_left(i,j)=delta_ne(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$        tmp1(i,j)=delta_ne_theta(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$        tmp2(i,j)=delta_ne_psi(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$        tmp3(i,j)=delta_ne_alpha(radcor_1d_array2(j),theta_start,tor_1d_array(i))/ne0
!!$     enddo
!!$  enddo !---testing end

  !  potential=den_i_left !assuming adiabatic electrons and quasineutrality, potential normalized by Te/e

  !smoothing is moved outside of this subroutine
!!$  if(filter_toroidal.eqv..true.) then !filter over the toroidal mode number, keeping the perturbation with desired toroidal mode number
!!$     call oned_fourier_transform1(potential,potential_dft,mtor,nflux2) !calculating 1d DFT of s(:,:) along the first dimension
!!$     call toroidal_filter(potential_dft,mtor,nflux2)
!!$     call oned_backward_fourier_transform1(potential_dft,potential,mtor,nflux2)
!!$  endif

!!$  if(filter_radial.eqv..true.) then !filter over the radial mode number, keeping only low-radial-harmonics of the perturbation
!!$     call oned_sine_transform2(potential,potential_dst,mtor,nflux2) !calculating 1d DST of s(:,:) along the second dimension
!!$     call radial_sine_filter_core(potential_dst,mtor,nflux2)
!!$     call oned_inverse_sine_transform2(potential_dst,potential,mtor,nflux2)
!!$  endif

  call radial_derivative  (potential(:,:,1),potential_psi)
  call toroidal_derivative(potential(:,:,1),potential_alpha)
  call theta_derivative(potential,potential_theta) !derivative along the magnetic field line
  !write(*,*) potential_theta(10,40),tmp1(10,40),'myid=',myid
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do i=1,mtor
     do j=1,nflux2
        jeq=j-1+j_low2
        ef_cyl_r(i,j)=-potential_psi(i,j)*grad_psi_r(ipol_eq,jeq)&
             & -potential_theta(i,j)*grad_theta_r(ipol_eq,jeq)-potential_alpha(i,j)*grad_alpha_r(ipol_eq,jeq)

        ef_cyl_z(i,j)=-potential_psi(i,j)*grad_psi_z(ipol_eq,jeq) &
             & -potential_theta(i,j)*grad_theta_z(ipol_eq,jeq)-potential_alpha(i,j)*grad_alpha_z(ipol_eq,jeq)
        !ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mc(ipol_eq,jeq)*grad_alpha_phi(ipol_eq,jeq) !wrong, an additional 1/R factor is wrongly included, a bug found on 2018-Jan.-4 evening
        ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mc(ipol_eq,jeq) !corrected
     enddo
  enddo


  normal=tu*kev/(ln*bn*vn_i*qu) !;write(*,*) 'normal=',normal
  !write(*,*) 'normal=',normal
  ef_cyl_r=ef_cyl_r*normal !normalized by the unit used in the code (Bn*vn_i)
  ef_cyl_z=ef_cyl_z*normal
  ef_cyl_phi=ef_cyl_phi*normal

  do j=1,nflux2
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

  ef_cyl_r_left(:,nflux2)=  0._p_
  ef_cyl_z_left(:,nflux2)=  0._p_
  ef_cyl_phi_left(:,nflux2)= 0._p_

  call communicate_field_value_between_neighbour_cells2() !use ER,Ephi, EZ
end subroutine potential_to_cylindrical_field

