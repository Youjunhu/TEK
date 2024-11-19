module filter_module
implicit none
contains


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
    use constants,only:p_
    use constants,only:twopi
    use magnetic_coordinates,only:toroidal_range,tor_1d_array,nsegment !computational toroidal region is 1/nsegment of the full torus
    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(in):: s_dft(0:m-1,n)
    real(p_),intent(out):: s(m,n)
    real(p_):: tor,a,b
    integer:: i,j,i_positive,i_negative
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
          tor=tor_1d_array(i)
          s(i,j)=a*cos(i_positive*twopi*tor/toroidal_range)+b*sin(i_positive*twopi*tor/toroidal_range)
       enddo
    enddo
  end subroutine toroidal_reconstruct_v0


  subroutine toroidal_reconstruct(s_dft,s,m,n) !manually reconstruct (instead of using inverse DFT) the real space function from the DFT data, the result was verified, which agrees with that given by inverse DFT
    use constants,only:p_
    use constants,only:twopi
    use magnetic_coordinates,only:toroidal_range,tor_1d_array,nsegment !computational toroidal region is 1/nsegment of the full torus
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

    tor(:)=i_positive*twopi*tor_1d_array(1:m)/toroidal_range
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
  use control_parameters,only:  radial_harmonics_included
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
  use control_parameters,only:  radial_harmonics_included
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(inout):: s(0:m-1,0:n-1)
  integer:: i,j
  do j=0,n-1
     if(j.gt.radial_harmonics_included) s(:,j)=0._p_
  enddo
end subroutine radial_sine_filter_core

!!$subroutine radial_sine_filter_em_field()
!!$  use constants,only:p_
!!$  use perturbation_field_matrix,only: ex=>ex_left,ey=>ey_left,epar=>epar_left !as input and output
!!$  use perturbation_field_matrix,only: mx=>mf_x_left,my=>mf_y_left,mpar=>mf_par_left !as input and output
!!$  use magnetic_coordinates,only: m=>mtor,n=>nflux2
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
  use magnetic_coordinates,only: radcor_1d_array2
  use control_parameters,only:  radial_harmonics_included
  implicit none
  integer,intent(in)::m,n
  real(p_),intent(in):: s_dst(m,0:n-1)
  real(p_),intent(out)::s(m,0:n-1)
  integer:: j,k
  real(p_):: radial_range,x,dx

!  dx=radcor_1d_array2(2)-radcor_1d_array2(1)
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
