module FFTW3
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  save
  type(C_PTR) :: plan_toroidal,plan_radial, plan_toroidal_backward,plan_radial_backward
  type(C_PTR) :: plan_sin,plan_sin_backward
  !integer*8:: plan_sin
  complex(p_),allocatable :: in1(:), out1(:)
  complex(p_),allocatable :: in2(:), out2(:)           
  real(p_),allocatable:: in3(:), out3(:)

contains
  subroutine  initialize_fft
    use magnetic_coordinates,only: m=>mtor,n=>nrad
    allocate(in1(0:m-1))
    allocate(out1(0:m-1))           
    allocate(in2(0:n-1))
    allocate(out2(0:n-1))

    !plan_toroidal = fftw_plan_dft_1d(m, in1,out1, FFTW_FORWARD,FFTW_ESTIMATE)
    !plan_toroidal = fftw_plan_dft_1d(m, in1,out1, FFTW_FORWARD,FFTW_measure)

    call dfftw_plan_dft_1d(plan_toroidal,m,in1,out1,FFTW_FORWARD,FFTW_measure)

    !plan_radial = fftw_plan_dft_1d(n, in2,out2, FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(plan_radial, n, in2,out2, FFTW_FORWARD,FFTW_MEASURE)

!call dfftw_plan_dft_1d(plan,N,in,out,FFTW_FORWARD,FFTW_ESTIMATE)

    call dfftw_plan_dft_1d(plan_toroidal_backward,m, in1,out1, FFTW_backward,FFTW_measure)
    plan_radial_backward = fftw_plan_dft_1d(n, in2,out2, FFTW_backward,FFTW_ESTIMATE)

    allocate(in3(n-2)) !
    allocate(out3(n-2))
    plan_sin =fftw_plan_r2r_1d(n-2, in3, out3, FFTW_RODFT00, FFTW_MEASURE)

!    call test_fft()
  end subroutine initialize_fft

!!$  subroutine test_fft()
!!$    use transform_module,only: oned_fourier_transform1, oned_backward_fourier_transform1
!!$    use magnetic_coordinates,only: m=>mtor,n=>nrad
!!$    real(p_):: source(m,n-2), source_tmp(m,n-2)
!!$    complex(p_):: source_dft(m,n-2)
!!$    integer:: i,j
!!$integer,parameter :: seed = 86456
!!$    call srand(seed)
!!$    do i=1,m
!!$       do j=1,n-2
!!$          source(i,j)= rand()
!!$       enddo
!!$    enddo
!!$    call oned_fourier_transform1(source,source_dft,m,n-2) !calculating DFT of source1(:,:) along the first dimension
!!$    call oned_backward_fourier_transform1(source_dft,source_tmp,m,n-2)
!!$write(*,*) maxval(source), maxval(source_tmp), maxval(source_tmp),maxval(abs(source-source_tmp))
!!$  end subroutine test_fft
end module FFTW3

module transform_module
implicit none

contains
  subroutine oned_fourier_transform1(s,s_fft,m,n) !calculating 1d DFT of s(:,:) along the first dimension
    use, intrinsic :: iso_c_binding
    use FFTW3,only: plan_toroidal, in1,out1
    use constants,only:p_
    implicit none
    include 'fftw3.f03'
    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1, n)
    complex(p_),intent(out):: s_fft(0:m-1, n)
    integer:: j

    do j=1,n  ! using openmp to parallize this loop gave wrong results (on Edison machine), strange!, possibly indicating fftw_execute_dft is not thread-safe.
       in1(:)=s(:,j) !copy in, meanwhile convert real array to complex array
       call fftw_execute_dft(plan_toroidal, in1(:), out1(:))     !Fourier transformation along the first dimension
       s_fft(:,j)=out1(:) !copy out
    enddo

!tests indicate the following is slower than the above, very interesting!
!!$    s_complex=s !convert real array to complex array
!!$    do j=1,n
!!$       call fftw_execute_dft(plan_toroidal, s_complex(:,j), s_fft(:,j))     !Fourier transformation along the first dimension
!!$    enddo

  end subroutine oned_fourier_transform1

  subroutine oned_fourier_transform1_parallel_version(s,s_fft,m,n) !calculating 1d DFT of s(:,:) along the first dimension
    use, intrinsic :: iso_c_binding
    use FFTW3,only: plan_toroidal, in1,out1
    use constants,only:p_
    use domain_decomposition,only: TCLR, ntube, grid_comm
    use mpi
    implicit none
    include 'fftw3.f03'
    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,n)
    complex(p_),intent(out):: s_fft(0:m-1,n)
    integer:: i,j,ierr
    logical,save:: is_first=.true.
    integer,save::my_part_start, my_part_end, spacing, my_range
    integer,allocatable,save:: recvcounts(:), displacement(:)
    complex(p_),allocatable,save:: my_s_fft(:,:)

    if (is_first.eqv..true.) then !needs to be genreated only for the first time and used multiple times
       spacing=n/ntube
       my_part_start=TCLR*spacing+1
       my_part_end=my_part_start+(spacing-1)
       if(TCLR.eq.(ntube-1)) my_part_end=n !the last process handles all the remainder part, in the case that n is not a perfect multiple of ntube
       my_range=my_part_end-my_part_start+1
       allocate(my_s_fft(0:m-1,my_range))
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:)=spacing*m
       recvcounts(ntube-1)= recvcounts(ntube-1)+(n-spacing*ntube)*m !last process contains additional elements
       do i=0,ntube-1
          displacement(i)=i*spacing*m
       enddo
       is_first=.false.
    endif

    do j=my_part_start,my_part_end  ! using openmp to parallize this loop gave wrong results (on Edison machine), strange!, possibly indicating fftw_execute_dft is not thread-safe.
       in1(:)=s(:,j) !copy in, meanwhile convert real array to complex array
       call fftw_execute_dft(plan_toroidal, in1(:), out1(:))     !Fourier transformation along the first dimension
       my_s_fft(:,j-my_part_start+1)=out1(:) !copy out
    enddo

    if(ntube.gt.1) then
       call mpi_allgatherv(my_s_fft, m*my_range, MPI_DOUBLE_COMPLEX, &
            &     s_fft, recvcounts, displacement, MPI_DOUBLE_COMPLEX, grid_comm, ierr)
    else
       s_fft=my_s_fft
    endif
  end subroutine oned_fourier_transform1_parallel_version

  subroutine oned_backward_fourier_transform1_parallel_version(field_fft,field,m,n) !radial task decomposion at a z gridpoint
    use, intrinsic :: iso_c_binding
    use FFTW3, only: plan_toroidal_backward
    use constants,only:p_
    use domain_decomposition,only: TCLR, ntube, grid_comm
    use mpi
    implicit none
    include 'fftw3.f03' 
    integer,intent(in) :: m,n
    complex(p_), intent(in) :: field_fft(0:m-1,n)
    real(p_), intent(out) :: field(0:m-1,n)
    complex(p_) :: in1(0:m-1), out1(0:m-1)  
    integer :: i, j, ierr
    logical, save :: is_first = .true.
    integer, save :: my_part_start, my_part_end, spacing, my_range
    integer, allocatable, save:: recvcounts(:), displacement(:)
    real(p_), allocatable, save :: my_field(:,:)

    if (is_first.eqv..true.) then !needs to be genreated only for the first time and used multiple times
       is_first = .false.
       spacing = n/ntube
       my_part_start = TCLR*spacing+1
       my_part_end = my_part_start+(spacing-1)
       !the last process handles all the remainder part, in the case that n is not a perfect multiple of ntube
       if(TCLR.eq.(ntube-1)) my_part_end=n 
       my_range=my_part_end-my_part_start+1
       allocate(my_field(0:m-1,my_range))
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:)=spacing*m
       recvcounts(ntube-1)= recvcounts(ntube-1)+(n-spacing*ntube)*m !last process contains additional elements
       do i=0,ntube-1
          displacement(i)=i*spacing*m
       enddo
    endif

    do j = my_part_start, my_part_end  !each processor handles its radial range
       in1(:)=field_fft(:,j) !copy in
       call fftw_execute_dft(plan_toroidal_backward, in1(:), out1(:)) !fourier transform along toroidal direction
       my_field(:,j-my_part_start+1) = real(out1(:))/m !copy out and inclued the 1/m factor
    enddo

    if(ntube.gt.1) then
       call mpi_allgatherv(my_field, m*my_range, MPI_DOUBLE, &
            &  field, recvcounts, displacement, MPI_DOUBLE, grid_comm, ierr)
    else
       field = my_field
    endif
  end subroutine oned_backward_fourier_transform1_parallel_version


  subroutine oned_backward_fourier_transform1(field_fft, field, m, n) 
    use, intrinsic :: iso_c_binding
    use FFTW3, only: plan_toroidal_backward
    use constants,only:p_
    implicit none
    include 'fftw3.f03' 
    integer, intent(in) :: m, n
    complex(p_), intent(in) :: field_fft(0:m-1, n)
    real(p_), intent(out):: field(0:m-1, n)
    complex(p_) :: in(0:m-1, n), out(0:m-1, n)  
    integer:: j

    in=field_fft
    !  plan = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_backward,FFTW_ESTIMATE)
    do j=1,n !fourier transform along the first direction
       call fftw_execute_dft(plan_toroidal_backward, in(:,j), out(:,j))
    enddo
    !  call fftw_destroy_plan(plan)  

    field=real(out)/m
  end subroutine oned_backward_fourier_transform1



  subroutine oned_fourier_transform2(s,s_fft,m,n) !calculating 1d DFT of s(:,:) along the second dimension
    use FFTW3
    use constants,only:p_
    implicit none

    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    type(C_PTR) :: plan
    integer:: i

    !Fourier transformation along the second dimension
    in=s
    plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan, in(i,:), out(i,:))
    enddo
    call fftw_destroy_plan(plan)  
    s_fft=out
  end subroutine oned_fourier_transform2


  subroutine twod_fourier_transform(s,s_fft,m,n) !calculating 2d DFT of source, tested
    use FFTW3
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    integer:: i,j

    !Fourier transformation along the first dimension
    in=s
    !plan_toroidal = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_FORWARD,FFTW_ESTIMATE)
    do j=0,n-1
       call fftw_execute_dft(plan_toroidal, in(:,j), out(:,j))
    enddo
    !call fftw_destroy_plan(plan_toroidal)  

    !  in=out

    !  call toroidal_filter(in,out,m,n)

!!$  write(*,*) 'given by fftw'
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j)),i=0,m-1)
!!$  enddo
!!$do j=0,n-1
!!$call my_fft(in(:,j),out(:,j),m)
!!$enddo
!!$ write(*,*) 'given by my_fft'
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j)),i=0,m-1)
!!$  enddo

    !Fourier transformation along the second dimension
    in=out
    !plan_radial = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan_radial, in(i,:), out(i,:))
    enddo
    !call fftw_destroy_plan(plan_radial)  
    s_fft=out

  end subroutine twod_fourier_transform


subroutine twod_fourier_transform_xz(s,s_fft,m,n) !calculating 2d DFT of source, tested
    use FFTW3
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    integer:: i,j
type(C_PTR) :: plan_z
    !Fourier transformation along the first dimension
    in=s
    do j=0,n-1
       call fftw_execute_dft(plan_radial, in(:,j), out(:,j))
    enddo

    !Fourier transformation along the second dimension
    in=out
    plan_z = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan_z, in(i,:), out(i,:))
    enddo
    call fftw_destroy_plan(plan_z)
    s_fft=out

  end subroutine twod_fourier_transform_xz


  subroutine twod_inverse_fourier_transform(field_fft,field,m,n) !tested
    use FFTW3
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    !  complex(C_DOUBLE_COMPLEX),intent(in):: field_fft(0:m-1,0:n-1)
    complex(p_),intent(in):: field_fft(0:m-1,0:n-1)
    real(p_),intent(out):: field(0:m-1,0:n-1)
    !  complex(C_DOUBLE_COMPLEX) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    integer:: i,j

    in=field_fft
    !plan_toroidal_backward = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_backward,FFTW_ESTIMATE)
    do j=0,n-1 !fourier transform along the first direction
       call fftw_execute_dft(plan_toroidal_backward, in(:,j), out(:,j))
    enddo
    !call fftw_destroy_plan(plan_toroidal_backward)  

    in=out
    !plan_radial_backward = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_backward,FFTW_ESTIMATE)
    do i=0,m-1 !fourier transform along the second direction
       call fftw_execute_dft(plan_radial_backward, in(i,:), out(i,:))
    enddo
    !call fftw_destroy_plan(plan_radial_backward)  

    field=real(out)/(m*n)
!!$  !check whether the results are the same
!!$  write(*,*) 'original data='
!!$  do j=0,n-1
!!$     write(*,*) (s(i,j),i=0,m-1)
!!$  enddo
!!$  write(*,*) 'DFT+backward-DFT data='
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j))/(m*n),i=0,m-1)
!!$  enddo
  end subroutine twod_inverse_fourier_transform


subroutine simple_fft(in,out,m) !for test
  use constants,only:p_
  use constants,only:twopi
  implicit none
  integer,intent(in):: m
  complex(p_),intent(in)::in(0:m-1)
  complex(p_),intent(out)::out(0:m-1)
  complex(p_),parameter::ii=(0.0_p_,1._p_)
  complex(p_):: sum
  integer:: i,ip

  do i=0,m-1
     sum=0._p_
     do ip=0,m-1
        sum=sum+in(ip)*exp(-twopi*ii/m*ip*i)
     enddo
     out(i)=sum
  enddo

end subroutine simple_fft

subroutine my_inverse_fft(in,out,m) !for test
  use constants,only:p_
  use constants,only:twopi
  implicit none
  integer,intent(in):: m
  complex(p_),intent(in)::in(0:m-1)
  complex(p_),intent(out)::out(0:m-1)
  complex(p_),parameter::ii=(0.0_p_,1._p_)
  complex(p_):: sum
  integer:: i,ip

  do i=0,m-1
     sum=0._p_
     do ip=0,m-1
        sum=sum+in(ip)*exp(twopi*ii/m*ip*i)
     enddo
     out(i)=sum
  enddo
out=out/m
end subroutine my_inverse_fft





subroutine oned_backward_fourier_transform2(field_fft,field,m,n) 
    use FFTW3
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  complex(p_),intent(in):: field_fft(0:m-1,0:n-1)
  real(p_),intent(out):: field(0:m-1,0:n-1)
  type(C_PTR) :: plan
  complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
  integer:: i

  in=field_fft
  plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_backward,FFTW_ESTIMATE)
  do i=0,m-1 !fourier transform along the second direction
     call fftw_execute_dft(plan, in(i,:), out(i,:))
  enddo
  call fftw_destroy_plan(plan)  
  field=real(out)/n
end subroutine oned_backward_fourier_transform2

subroutine oned_sine_transform2(s,s_dst,m,n) !calculating 1d DST of s(:,:) along the second dimension
  use FFTW3, only: plan_sin, in3, out3
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: s(0:m-1,0:n-1)
  real(p_),intent(out):: s_dst(0:m-1,0:n-1)
 ! type(C_PTR) :: plan !a pointer, which needs to be distoried, otherwise may cause memory leak
  integer:: i
  !plan =fftw_plan_r2r_1d(n, in(1,:), s_dst(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     in3(:)=s(i,:)
     !call dfftw_execute_r2r(plan_sin,s(i,:),s_dst(i,:)) !DST along the second dimension
     call dfftw_execute_r2r(plan_sin,in3,out3) !DST along the second dimension
     s_dst(i,:)=out3(:)
  enddo
!  call fftw_destroy_plan(plan)  

end subroutine oned_sine_transform2

subroutine oned_inverse_sine_transform2(s_dst,s,m,n) !computing 1d inverse DST of s(:,:) along the second dimension
  use FFTW3, only: plan_sin, in3, out3
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: s_dst(0:m-1,0:n-1)
  real(p_),intent(out):: s(0:m-1,0:n-1)
!  real(p_) :: in(0:m-1,0:n-1)
!  type(C_PTR) :: plan
  integer:: i

  !in=s_dst
  !plan =fftw_plan_r2r_1d(n, in(1,:), s(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     in3(:)=s_dst(i,:)
    ! call dfftw_execute_r2r(plan_sin,s_dst(i,:),s(i,:))   !DST along the second dimension
     call dfftw_execute_r2r(plan_sin,in3,out3)   !DST along the second dimension
     s(i,:)=out3(:)
  enddo
  !call fftw_destroy_plan(plan)  
  s=s/(2*(n+1))
end subroutine oned_inverse_sine_transform2


subroutine dst_dft(s,s_spectrum,m,n)
 ! use FFTW3,only:
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'
  integer,intent(in):: m,n
  real(p_),intent(in):: s(0:m-1,0:n-1)
  complex(p_),intent(out):: s_spectrum(0:m-1,0:n-1)
  real(p_):: in1(0:m-1,0:n-1),s_dst(0:m-1,0:n-1)
  complex(p_):: in2(0:m-1,0:n-1)
  type(C_PTR) :: plan
  integer:: i,j

  in1=s
  plan =fftw_plan_r2r_1d(n,   in1(1,:), s_dst(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     call dfftw_execute_r2r(plan,in1(i,:), s_dst(i,:)) !DST along the second dimension
  enddo
  call fftw_destroy_plan(plan)  

  in2=s_dst
  plan = fftw_plan_dft_1d(m,  in2(:,1), s_spectrum(:,1), FFTW_FORWARD,FFTW_ESTIMATE)
  do j=0,n-1
     call fftw_execute_dft(plan, in2(:,j), s_spectrum(:,j))   !DFT along the first dimension
  enddo
  call fftw_destroy_plan(plan)
!s_spectrum=s_spectrum/((n+1)*m) !the corresponding expansion coeficient is the dst_dft devided by (n+1)*m, see my notes on Fourier analysis
end subroutine dst_dft


end module transform_module


module fourn_module
contains
  subroutine twod_fourier_transform_nr(s,s_fft,m,n) !tested on Feb 27, 2018, note that the convension of FFT and inverse FFT definition is differnt between FFTW and Numerical-recipe book. I forgot about this when comparing the result given by this subroutine with that by FFTW subroutine, and it took about one hour for me to finally figure out why they are different.
    use constants,only:p_
!    use transform_module,only: my_fft
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(in):: s(m,n)
    complex(p_),intent(out):: s_fft(m,n)
    real(p_):: data(2*m*n)
    integer,parameter:: ndim=2
    complex(p_),parameter:: ii=(0._p_,1._p_)
    complex(p_):: in(m,n),out(m,n)
    integer:: i,j,k,nn(ndim),isign

    do j=1,n
       do i=1,m
          k=i+m*(j-1)
          data(2*k-1)=s(i,j)  !real part
          data(2*k)=0._p_ !imaginary part
       enddo
    enddo

!---for testing--------passed successfully
!!$    in=s
!!$    do j=1,n
!!$       call my_fft(in(:,j),out(:,j),m)
!!$    enddo
!!$    !Fourier transformation along the second dimension
!!$    in=out
!!$    do i=1,m
!!$       call my_fft(in(i,:),out(i,:),n)
!!$    enddo
!!$    s_fft=out
!----for testing----
    nn(1)=m
    nn(2)=n
    isign=-1 !isign=-1 corresponds to inverse FFT in numerical-recipe-book, but corresponds to the forward FFT in FFTW
    call fourn(data,nn,ndim,isign)

    do j=1,n
       do i=1,m
          k=i+m*(j-1)
          s_fft(i,j)=data(2*k-1)+ii*data(2*k)
       enddo
    enddo
 end subroutine twod_fourier_transform_nr


subroutine twod_inverse_fourier_transform_nr(s_fft,s,m,n) 
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  complex(p_),intent(in):: s_fft(m,n)
  real(p_),intent(out):: s(m,n)
  real(p_):: data(2*m*n)
  integer,parameter:: ndim=2
  integer:: i,j,k,nn(ndim),isign

  do j=1,n
     do i=1,m
        k=i+m*(j-1)
        data(2*k-1)=real(s_fft(i,j))  !real part
        data(2*k)=imag(s_fft(i,j)) !imaginary part
     enddo
  enddo

  nn(1)=m
  nn(2)=n
  isign=1 !isign=1 corresponds to FFT in numerical-recipe-book, but corresponds the inverse FFT in FFTW
  call fourn(data,nn,ndim,isign)
  do i=1,m
     do j=1,n
        k=i+m*(j-1)
        s(i,j)=data(2*k-1)
     enddo
  enddo
s=s/(m*n)
end subroutine twod_inverse_fourier_transform_nr

SUBROUTINE fourn(data,nn,ndim,isign)
  use constants,only:p_
  implicit none
  INTEGER isign,ndim,nn(ndim)
  REAL(p_):: data(*)
  !Replaces data by its ndim -dimensional discrete Fourier transform, if isign is input as 1.
  !nn(1:ndim) is an integer array containing the lengths of each dimension (number of complex values), which MUST all be powers of 2. 
  !data is a real array of length twice the product of these lengths, in which the data are stored as in a multidimensional complex FORTRAN array
  !If isign is input as -1, data is replaced by its inverse transform times the product of the lengths of all dimensions. 
  !code from the numerical recipe book
  INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,&
       & ip3,k1,k2,n,nprev,nrem,ntot
  REAL(p_):: tempi,tempr
  real(p_):: theta,wi,wpi,wpr,wr,wtemp !Double precision for trigonometric re-currences.

  ntot=1
  do idim=1,ndim !Compute total number of complex values.
     ntot=ntot*nn(idim)
  enddo
  nprev=1
  do idim=1,ndim !Main loop over the dimensions.
     n=nn(idim)
     nrem=ntot/(n*nprev)
     ip1=2*nprev
     ip2=ip1*n
     ip3=ip2*nrem
     i2rev=1
     do  i2=1,ip2,ip1 !This is the bit-reversal section of the routine.
        if(i2.lt.i2rev)then
           do  i1=i2,i2+ip1-2,2
              do  i3=i1,ip3,ip2
                 i3rev=i2rev+i3-i2
                 tempr=data(i3)
                 tempi=data(i3+1)
                 data(i3)=data(i3rev)
                 data(i3+1)=data(i3rev+1)
                 data(i3rev)=tempr
                 data(i3rev+1)=tempi
              enddo
           enddo
        endif
        ibit=ip2/2
1       if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
           i2rev=i2rev-ibit
           ibit=ibit/2
           goto 1
        endif
        i2rev=i2rev+ibit
     enddo
     ifp1=ip1 !Here begins the Danielson-Lanczos section of the routine.
2    if(ifp1.lt.ip2)then
        ifp2=2*ifp1
        theta=isign*6.28318530717959d0/(ifp2/ip1) !Initialize for the trig. recur-rence.
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do  i3=1,ifp1,ip1
           do  i1=i3,i3+ip1-2,2
              do  i2=i1,ip3,ifp2
                 k1=i2 !Danielson-Lanczos formula:
                 k2=k1+ifp1
                 tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                 tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                 data(k2)=data(k1)-tempr
                 data(k2+1)=data(k1+1)-tempi
                 data(k1)=data(k1)+tempr
                 data(k1+1)=data(k1+1)+tempi
              enddo
           enddo
           wtemp=wr !Trigonometric recurrence.
           wr=wr*wpr-wi*wpi+wr
           wi=wi*wpr+wtemp*wpi+wi
        enddo
        ifp1=ifp2
        goto 2
     endif
     nprev=n*nprev
  enddo
  return
END SUBROUTINE fourn
end module fourn_module
