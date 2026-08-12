module filter_module
  implicit none
  integer, parameter ::  radial_harmonics_included = 25
contains

  subroutine surface_average_of_n0(s, nx)
    ! Retain only the magnetic surface average of n=0 harmonic
    use constants, only: p_, twopi
    use magnetic_coordinates, only: nrad, nz => mpol2, jacobian_av, abs_jacobian
    use domain_decomposition, only: tube_comm, ipol_eq, dtheta2
    use mpi
    implicit none
    integer, intent(in) :: nx !radial
    complex(p_), intent(inout) :: s(nx) ! Amplitude of n=0 Fourier harmonic
    real(p_) :: my_sz(nx), sz(nx, 0:nz-1), sx(nx)
    integer :: j, ierr

    my_sz(:) = real(s(:))*abs_jacobian(ipol_eq, 2:nrad-1)
    call mpi_allgather(my_sz(:), nx, MPI_real8, sz, nx, MPI_real8, tube_comm, ierr)

    do j = 1, nx
       sx(j) = sum(sz(j,:))*dtheta2/(jacobian_av(j+1)*twopi)
    enddo

    s(:) = sx(:) !store the result in the complex array.
  end subroutine 
  
  subroutine mfilter_for_n0(s, n)
    ! The n=0 mode is often numerically unstable if we keep all m harmonics
    ! m filtering (usually retain only m=0) is needed to make the simulation stable
    use constants, only: p_, ii, twopi, pi
    use magnetic_coordinates, only: nz => mpol2
    use domain_decomposition, only: tclr, gclr, ntube, grid_comm, tube_comm
    use mpi
    implicit none
    integer, intent(in) :: n !radial
    complex(p_), intent(inout) :: s(n)
    integer, parameter :: mupp = 0 !poloidal harmonics m in [-mupp:mupp], basis funtions: exp(ii*m*theta)
    logical, save :: is_first = .true.
    integer, save :: my_start, my_end, spacing, my_range
    integer, allocatable, save :: recvcounts(:), displacement(:)
    complex(p_), allocatable, save :: wexp1(:,:), wexp2(:)
    complex(p_), allocatable :: coef(:,:)
    real(p_), allocatable :: my_sz(:), sz(:,:), my_sx(:), sx(:)
    integer :: i, j, k, ierr

    if (is_first .eqv. .true.) then !needs to be genreated only for the first time
       is_first = .false.
       spacing = n/ntube
       my_start = TCLR*spacing + 1
       my_end = my_start + (spacing-1)
       if(TCLR == ntube-1) my_end = n !the last process handles all the remainder if n is not a perfect multiple of ntube
       my_range = my_end - my_start + 1
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:) = spacing
       recvcounts(ntube-1) = recvcounts(ntube-1)+(n-spacing*ntube) !last process contains additional elements
       do i = 0, ntube-1
          displacement(i) = i * spacing
       enddo
       allocate(wexp1(0:nz-1, -mupp:mupp))
       allocate(wexp2(-mupp:mupp))
       do k = -mupp, mupp
          do i = 0, nz-1
             wexp1(i,k) = exp(-ii*k*(-pi+i*twopi/nz))
          enddo
          wexp2(k) = exp(+ii*k*(-pi+gclr*twopi/nz))
       enddo
    endif

    allocate(my_sz(my_start:my_end))
    allocate(sz(my_start:my_end, 0:nz-1))
    allocate(coef(-mupp:mupp, my_start:my_end))
    allocate(my_sx(my_start:my_end))
    allocate(sx(n))

    my_sz(:) = real(s(my_start:my_end))
    call mpi_allgather(my_sz(:), my_range, MPI_real8, sz, my_range, MPI_real8, tube_comm, ierr)

    do j = my_start, my_end
       do k = -mupp, mupp
          coef(k,j) = sum(sz(j,:)*wexp1(:,k)) ! Fourier expansion coefficient
       enddo
       my_sx(j) = real(sum(coef(:,j)*wexp2(:))/nz)  ! Inverse Fourier Transform, and convert the type to real
    enddo

    call mpi_allgatherv(my_sx, my_range, MPI_real8, &
         &     sx, recvcounts, displacement, MPI_real8, grid_comm, ierr)

    s(:) = sx(:) !store the result in the complex array.
  end subroutine mfilter_for_n0

  subroutine mfilter_for_each_n(s, nx, kn)
    use constants, only: p_, ii, twopi, pi
    use magnetic_coordinates, only: nsegment, nrad, nz => mpol2
    use radial_module, only: qrad
    use domain_decomposition, only: tclr, gclr, ntube, grid_comm, tube_comm, theta_start
    use mpi
    implicit none
    integer, intent(in) :: nx, kn !radial, toroidal_harmonic
    complex(p_), intent(inout) :: s(nx)
    integer, parameter :: mupp = 5 
    complex(p_) :: coef(-mupp:mupp)
    logical, save :: is_first = .true.
    integer, save :: my_start, my_end, span, my_range
    integer, allocatable, save :: recvcounts(:), displacement(:)
    complex(p_), allocatable, save :: wexp1(:,:,:), wexp2(:,:)
    complex(p_), allocatable :: my_sz(:), sz(:,:), my_filtered(:), filtered(:)
    integer :: i, j, km, ierr, nq0
    real(p_) :: th

    if (is_first .eqv. .true.) then !do this only the first time of entering this subroutine
       is_first = .false.
       span = nx/ntube   ! radial decomposition
       my_start = TCLR*span + 1
       my_end = my_start + (span-1)
       if(TCLR == ntube-1) my_end = nx !last process handles all the remainder if nx is not a perfect multiple of ntube
       my_range = my_end - my_start + 1
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:) = span
       recvcounts(ntube-1) = recvcounts(ntube-1) + (nx-span*ntube) !last process contains additional elements
       do i = 0, ntube-1
          displacement(i) = i * span
       enddo
       allocate(wexp1(0:nz-1, -mupp:mupp, my_start:my_end))
       allocate(wexp2(-mupp:mupp, my_start:my_end))

       do j = my_start, my_end
          nq0 = nint(nsegment*kn*qrad(j+1))
          do km = -mupp, mupp

             do i = 0, nz-1
                th = -pi + i*twopi/nz
                wexp1(i, km, j) = exp(-ii*(km - nq0)*th)
             enddo

             th = -pi + gclr*twopi/nz
             wexp2(km, j) = exp(+ii*(km - nq0)*th)

          enddo
       enddo
    endif

    allocate(my_sz(my_start:my_end))
    allocate(sz(my_start:my_end, 0:nz-1))
    my_sz(:) = s(my_start:my_end) / exp(ii*nsegment*kn*qrad(my_start+1:my_end+1)*(theta_start+pi))
    call mpi_allgather(my_sz(:), my_range, MPI_complex16, sz, my_range, MPI_complex16, tube_comm, ierr)

    allocate(my_filtered(my_start:my_end))
    allocate(filtered(nx))
    do j = my_start, my_end
       do km = -mupp, mupp
          coef(km) = sum(sz(j,:)*wexp1(:,km, j))/nz ! Fourier expansion coefficient
       enddo
       my_filtered(j) = sum(coef(:)*wexp2(:,j))  ! Reconstruction (i.e., Inverse Fourier Transform)
    enddo
    call mpi_allgatherv(my_filtered, my_range, MPI_complex16, &
         &     filtered, recvcounts, displacement, MPI_complex16, grid_comm, ierr)

    s(:) = filtered(:) * exp(ii*nsegment*kn*qrad(2:nrad-1)*(theta_start+pi))
  end subroutine mfilter_for_each_n


  subroutine toroidal_filter(field,m,n)
    use constants,only:p_
    use magnetic_coordinates,only:nsegment !computational toroidal region is 1/nsegment of the full torus
    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(inout)::   field(0:m-1,0:n-1)
    integer:: nharmonic
    integer:: i,i_negative_wn

    nharmonic=1 !harmonic that needs to be kept in the computational toroidal segment
    !write(*,*) 'nharmonic=',nharmonic
    i_negative_wn=m-nharmonic
    if(nharmonic.eq.0) i_negative_wn=0

    do i=0,m-1
       if((i.ne. nharmonic) .and. (i.ne.i_negative_wn)) field(i,:)=(0._p_,0._p_)
    enddo

!!$    do j=0,n-1
!!$     do i=0,m/2 !scanning over the positive wavenumber
!!$        if(i.eq.nharmonic) then 
!!$           out(i,j)=in(i,j) !positive wavenumber
!!$           i_negative_wn=m-i
!!$           if(i.eq.0) i_negative_wn=0
!!$           out(i_negative_wn,j)=in(i_negative_wn,j) !negative or zero wavenumber
!!$        endif
!!$     enddo
!!$    enddo
  end subroutine toroidal_filter


  subroutine toroidal_reconstruct_v0(s_dft,s,m,n) !manually reconstruct (instead of using inverse DFT) the real space function from the DFT data, the result was verified, which agrees with that given by inverse DFT
    use constants, only: p_, twopi
    use magnetic_coordinates, only: toroidal_range, ygrid, nsegment !computational toroidal region is 1/nsegment of the full torus
    implicit none
    integer, intent(in) :: m,n
    complex(p_), intent(in) :: s_dft(0:m-1,n)
    real(p_), intent(out) :: s(m,n)
    real(p_) :: tor, a, b
    integer :: i, j, i_positive, i_negative
    !complex(p_):: ii=(0._p_,1._p_)
    i_positive=1 !harmonic that needs to be kept in the computational toroidal segment
    i_negative=m-i_positive
    ! write(*,*) 'i_positive=',i_positive,'i_negative=',i_negative
    !s=0._p_
    do j=1,n
       a=2*real(s_dft(i_positive,j))/m
       b=-2*imag(s_dft(i_positive,j))/m
       !a=s_dft(i_positive,j)+s_dft(i_negative,j)
       !b=(s_dft(i_positive,j)-s_dft(i_negative,j))*ii
       do i=1,m
          tor=ygrid(i)
          s(i,j)=a*cos(i_positive*twopi*tor/toroidal_range)+b*sin(i_positive*twopi*tor/toroidal_range)
       enddo
    enddo
  end subroutine toroidal_reconstruct_v0


  subroutine toroidal_reconstruct(s_dft,s,m,n) !manually reconstruct (instead of using inverse DFT) the real space function from the DFT data, the result was verified, which agrees with that given by inverse DFT
    use constants,only:p_
    use constants,only:twopi
    use magnetic_coordinates,only:toroidal_range,ygrid,nsegment !computational toroidal region is 1/nsegment of the full torus
    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(in):: s_dft(0:m-1,n)
    real(p_),intent(out):: s(m,n)
    real(p_):: a(n),b(n),tor(m),cos_tor(m),sin_tor(m)
    integer:: i,j,i_positive,i_negative
    !complex(p_):: ii=(0._p_,1._p_)
    i_positive=1 !harmonic that needs to be kept in the computational toroidal segment
    ! i_negative=m-i_positive
    ! write(*,*) 'i_positive=',i_positive,'i_negative=',i_negative
    !s=0._p_

    a(:)=2*real(s_dft(i_positive,:))/m
    b(:)=-2*imag(s_dft(i_positive,:))/m

    tor(:)=i_positive*twopi*ygrid(1:m)/toroidal_range
    cos_tor(:)=cos(tor(:))
    sin_tor(:)=sin(tor(:))

    !$omp parallel do
    do j=1,n
       do i=1,m
          s(i,j)=a(j)*cos_tor(i)+b(j)*sin_tor(i)
       enddo
    enddo
    !$omp end parallel do
  end subroutine toroidal_reconstruct


  subroutine radial_fourier_filter(field,m,n)
    use constants,only:p_

    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(inout)::   field(0:m-1,0:n-1)
    complex(p_):: out(0:m-1,0:n-1)
    integer:: i,j,j_negative_wn

    out=(0._p_,0._p_) !initialized to zero
    do i=0,m-1
       do j=0,n/2  !scanning over the positive wavenumber
          if(j.le.radial_harmonics_included) then 
             out(i,j)=field(i,j) !positive wavenumber
             j_negative_wn=n-j
             if(j.eq.0) j_negative_wn=0
             out(i,j_negative_wn)=field(i,j_negative_wn) !negative or zero wavenumber
          endif
       enddo
    enddo

    field=out
  end subroutine radial_fourier_filter

  subroutine radial_sine_filter_core(s,m,n)
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(inout):: s(0:m-1,0:n-1)
    integer:: i,j
    do j=0,n-1
       if(j.gt.radial_harmonics_included) s(:,j)=0._p_
    enddo
  end subroutine radial_sine_filter_core

  subroutine radial_sine_high_pass(s, m, n)
    use constants, only: p_
    implicit none
    integer, intent(in) :: m, n
    real(p_), intent(inout) :: s(0:m-1, 0:n-1)
    integer, parameter :: jcut = 2
    integer :: i, j
    do j = 0, n-1
       if(j < jcut) s(:, j) = 0._p_
    enddo
  end subroutine radial_sine_high_pass

  
!!$subroutine radial_sine_filter_em_field()
!!$  use constants,only:p_
!!$  use perturbation_field,only: ex=>ex_left,ey=>ey_left,epar=>epar_left !as input and output
!!$  use perturbation_field,only: mx=>mf_x_left,my=>mf_y_left,mpar=>mf_par_left !as input and output
!!$  use magnetic_coordinates,only: m=>mtor,n=>nrad
!!$  use transform_module
!!$  implicit none
!!$  real(p_)::    epar_dst(m+1,n), ex_dst(m+1,n),  ey_dst(m+1,n)
!!$  real(p_)::    mpar_dst(m+1,n), mx_dst(m+1,n),  my_dst(m+1,n)
!!$
!!$  call oned_sine_transform2(ex,ex_dst,m+1,n) 
!!$  call oned_sine_transform2(ey,ey_dst,m+1,n) 
!!$  call oned_sine_transform2(epar,epar_dst,m+1,n) 
!!$  call radial_sine_filter_core(ex_dst,m+1,n)
!!$  call radial_sine_filter_core(ey_dst,m+1,n)
!!$  call radial_sine_filter_core(epar_dst,m+1,n)
!!$  call oned_inverse_sine_transform2(ex_dst,ex,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(ey_dst,ey,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(epar_dst,epar,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$
!!$  call oned_sine_transform2(mx,mx_dst,m+1,n) 
!!$  call oned_sine_transform2(my,my_dst,m+1,n) 
!!$  call oned_sine_transform2(mpar,mpar_dst,m+1,n) 
!!$  call radial_sine_filter_core(mx_dst,m+1,n)
!!$  call radial_sine_filter_core(my_dst,m+1,n)
!!$  call radial_sine_filter_core(mpar_dst,m+1,n)
!!$  call oned_inverse_sine_transform2(mx_dst,mx,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(my_dst,my,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(mpar_dst,mpar,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$
!!$
!!$end subroutine radial_sine_filter_em_field


  subroutine radial_sine_filter(s)
    use constants,only:p_
    use transform_module
    implicit none
    real(p_),intent(inout):: s(:,:)
    real(p_):: s_dst(size(s,1),size(s,2))
    integer:: m,n

    m=size(s,1)
    n=size(s,2)
    call oned_sine_transform2(s,s_dst,m,n) 
    call radial_sine_filter_core(s_dst,m,n)
    call oned_inverse_sine_transform2(s_dst,s,m,n) !computing 1d inverse DST of s(:,:) along the second dimension
  end subroutine radial_sine_filter


  subroutine radial_sine_reconstruct(s_dst,s,m,n)
    use constants,only:p_
    use constants,only: pi
    use magnetic_coordinates,only: xgrid
    implicit none
    integer,intent(in)::m,n
    real(p_),intent(in):: s_dst(m,0:n-1)
    real(p_),intent(out)::s(m,0:n-1)
    integer:: j,k
    real(p_):: radial_range,x,dx

    !  dx=xgrid(2)-xgrid(1)
    ! radial_range=dx*(n+1)
    s=0._p_
    do j=0,n-1
       !   x=dx*(j+1)
       do k=0,radial_harmonics_included
          s(:,j)=s(:,j)+s_dst(:,k)*sin((k+1)*pi*(j+1)/(n+1)) !see the formula in my notes on Fourier analysis
          !s(:,j)=s(:,j)+s_dst(:,k)*sin((k+1)*pi*x/radial_range) !see the formula in my notes on Fourier analysis
       enddo
    enddo
    s=s/(n+1)
  end subroutine radial_sine_reconstruct

  subroutine toroidal_fourier_radial_sine_filter(s,m,n)
    use transform_module,only: dst_dft
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(inout):: s(0:m-1,0:n-1)
    complex(p_):: s_spectrum(0:m-1,0:n-1)
    real(p_):: s_tmp(0:m-1,0:n-1)

    call dst_dft(s,s_spectrum,m,n)
    call toroidal_reconstruct(s_spectrum,s_tmp,m,n)
    call radial_sine_reconstruct(s_tmp,s,m,n)

  end subroutine toroidal_fourier_radial_sine_filter

  subroutine filter_basic(field)
    use constants,only:p_
    use transform_module
    implicit none
    real(p_),intent(inout):: field(:,:)
    complex(p_):: field_dft(size(field,1),size(field,2))
    real(p_):: field_dst(size(field,1),size(field,2))
    integer:: m,n
    m=size(field,1)
    n=size(field,2)

    call oned_fourier_transform1(field,field_dft,m,n) 
    call toroidal_reconstruct(field_dft,field,m,n)

    call oned_sine_transform2(field,field_dst,m,n) 
    call radial_sine_filter_core(field_dst,m,n)
    call oned_inverse_sine_transform2(field_dst,field,m,n) !computing 1d inverse DST of s(:,:) along the second dimension

  end subroutine filter_basic


  subroutine toroidal_fourier_radial_fourier_filter(s,m,n)
    !use fourn_module,only: twod_fourier_transform_nr,twod_inverse_fourier_transform_nr
    use transform_module
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(inout):: s(0:m-1,0:n-1)
    complex(p_):: s_spectrum(0:m-1,0:n-1)

    !  call twod_fourier_transform_nr(s,s_spectrum,m,n)
    call twod_fourier_transform(s,s_spectrum,m,n)
    call toroidal_filter(s_spectrum,m,n)
    call radial_fourier_filter(s_spectrum,m,n)
    !  call twod_inverse_fourier_transform_nr(s_spectrum,s,m,n) 
    call twod_inverse_fourier_transform(s_spectrum,s,m,n) 

  end subroutine toroidal_fourier_radial_fourier_filter

end module filter_module


module smoothing_module
  use constants, only: p_, one, two, six
  use communication_connection, only: get_nearby_field_along_field_line
  implicit none
  real(p_), parameter :: weight1 = 0.5_p_, weight2 = -one/six !got to know these values from the GEM code

contains

  subroutine smoothing_along_field_line_core(a)
    real(p_), intent(inout) :: a(:,:)
    real(p_) :: a_left(size(a,1), size(a,2)), a_right(size(a,1),size(a,2))

    call get_nearby_field_along_field_line(a, a_left, a_right) !get value of field on the two grids that are to the left/rightof the present grid

    a=(weight1*a_right + a + weight1*a_left)/(one+two*weight1) !smoothing using weight1

    !the smoothing can be appled multiple times (possibley with different weights):
    call get_nearby_field_along_field_line(a, a_left, a_right)

    a=(weight2*a_right + a +weight2*a_left)/(one+two*weight2) !smoothing using weight2

  end subroutine smoothing_along_field_line_core


  subroutine smoothing_along_field_line_core5(a)
    real(p_), intent(inout) :: a(:,:)
    real(p_) :: a_left(size(a,1),size(a,2)), a_right(size(a,1),size(a,2))
    real(p_) :: a_left2(size(a,1),size(a,2)), a_right2(size(a,1),size(a,2))

    call get_nearby_field_along_field_line(a, a_left, a_right, a_left2, a_right2) 

    a=(-a_left2 + 4*a_left + 10*a + 4*a_right -a_right2)/16.

  end subroutine smoothing_along_field_line_core5

end module smoothing_module
