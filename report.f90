module report_module
  implicit none
contains
subroutine report(t)
  use magnetic_coordinates,only: mtor,nflux2
!  use perturbation_field_matrix, only:ex=>ex_left,ey=>ey_left,epar=>epar_left

!  use perturbation_field_matrix,only:source_e1,source_e2 !,jper_x_i_left,jper_y_i_left
!  use perturbation_field_matrix,only:source1,source2,source3,source_faraday1, source_faraday2,source_faraday3
  use constants,only:p_

  real(p_),intent(in):: t
  integer:: i,j

i=mtor/2
j=nflux2/2

!write(*,'(7(1pe14.4))') t,ex(i,j),ey(i,j),epar(i,j),mf_x(i,j),mf_y(i,j),mf_par(i,j)
!write(*,*) mf_x(i,j),mf_y(i,j),mf_par(i,j)
!write(*,*) jper_x_i_left(mtor/2,nflux2/2),jper_y_i_left(mtor/2,nflux2/2)!,source_e1(mtor/2,nflux2/2),source_e2(mtor/2,nflux2/2)
!write(*,*) source1(mtor/3,nflux2/3),source2(mtor/3,nflux2/3),source3(mtor/3,nflux2/3)
!write(*,*) source_faraday1(mtor/3,nflux2/3),source_faraday2(mtor/3,nflux2/3),source_faraday3(mtor/3,nflux2/3)
end subroutine report


subroutine mode_evolution_analysis(t)
  use constants,only: one
  use constants,only:p_
  use magnetic_coordinates,only: dtor,dradcor,m=>mtor,n=>nflux2
!  use perturbation_field_matrix,only: ex_left
  use domain_decomposition,only:myid,numprocs
  use transform_module,only: twod_fourier_transform
  implicit none
  real(p_),intent(in):: t
  real(p_):: a(0:m-1,0:n-1)
  complex(p_):: a_fft(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative,jp,jn
  do i=0,m-1
     do j=0,n-1
 !       a(i,j)=ex_left(i+1,j+1)
     enddo
  enddo
  call twod_fourier_transform(a,a_fft,m,n)

ipositive=1
inegative=m-ipositive
!!$jp=9
!!$jn=n-jp
!!$if(myid.eq.1) write(*,*) t,real(a_fft(ipositive,jp)),imag(a_fft(ipositive,jp)),real(a_fft(inegative,jn)),imag(a_fft(inegative,jn))


! write(*,'(13(1pe14.4))') t,(real(a_fft(ipositive,j)),imag(a_fft(ipositive,j)),j=0,5)


end subroutine mode_evolution_analysis


subroutine mode_evolution_analysis2(t) !for adiabatic electron model, using DFT for both toroidal and raial direction
  use constants,only: one
  use constants,only:p_
  use magnetic_coordinates,only: m=>mtor,n=>nflux2
  use perturbation_field_matrix,only: ef_cyl_phi_left,ef_cyl_r_left,ef_cyl_z_left
  use perturbation_field_matrix,only: ef_cyl_phi_right,ef_cyl_r_right,ef_cyl_z_right
  use transform_module,only: twod_fourier_transform
  implicit none
  real(p_),intent(in):: t
  real(p_):: a(0:m-1,0:n-1)
  complex(p_):: a_fft(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative

  do i=0,m-1
     do j=0,n-1
        a(i,j)=ef_cyl_phi_left(i+1,j+1)
     enddo
  enddo
  call twod_fourier_transform(a,a_fft,m,n)

  ipositive=1
  inegative=m-ipositive

  write(*,'(20(1pe14.4))') t,(real(a_fft(ipositive,j)),imag(a_fft(ipositive,j)),j=0,4),&
       & ef_cyl_phi_left(m/2,n/2)

end subroutine mode_evolution_analysis2


subroutine mode_evolution_analysis3(t,a,m,n,file_unit) 
  use constants,only: one
  use constants,only:p_
  use transform_module,only: twod_fourier_transform,dst_dft
  use fourn_module,only: twod_fourier_transform_nr
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: a(0:m-1,0:n-1)
  real(p_),intent(in):: t
  integer,intent(in):: file_unit
  complex(p_):: a_spectrum(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative

!m=size(a,1)
!n=size(a,2)
!call dst_dft(a,a_spectrum,m,n) !using DST for radial direction and DFT for toroidal direction
!a_spectrum=a_spectrum/((n+1)*m) !the corresponding expansion coeficient is the dst_dft devided by (n+1)*m, see my notes on Fourier analysis
call twod_fourier_transform(a,a_spectrum,m,n) 
!call twod_fourier_transform_nr(a,a_spectrum,m,n) 
a_spectrum=a_spectrum/(n*m) 

  ipositive=1
  !inegative=m-ipositive

  write(file_unit,'(20(1pe20.8))') t,(real(a_spectrum(ipositive,j)),imag(a_spectrum(ipositive,j)),j=0,3)

end subroutine mode_evolution_analysis3


subroutine mode_evolution_analysis4(t,a,m,n,file_unit) 
  use constants,only: one
  use constants,only:p_
  use transform_module
!  use fourn_module,only: twod_fourier_transform_nr
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: a(0:m-1,0:n-1)
  real(p_),intent(in):: t
  integer,intent(in):: file_unit
  real(p_):: a_dst(0:m-1,0:n-1)
  complex(p_):: a_spectrum(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative

!m=size(a,1)
!n=size(a,2)
!call dst_dft(a,a_spectrum,m,n) !using DST for radial direction and DFT for toroidal direction
!a_spectrum=a_spectrum/((n+1)*m) !the corresponding expansion coeficient is the dst_dft devided by (n+1)*m, see my notes on Fourier analysis

!  call oned_sine_transform2(a,a_dst,m,n) 
!  call oned_fourier_transform1(a_dst,a_spectrum,m,n) 

call dst_dft(a,a_spectrum,m,n)
a_spectrum=a_spectrum/(2*(n+1)*m)

  ipositive=1
  !inegative=m-ipositive

  write(file_unit,'(20(1pe20.8))') t,(real(a_spectrum(ipositive,j)),imag(a_spectrum(ipositive,j)),j=0,3)

!  write(file_unit,'(2(1pe20.8))') t, a(m/2,n/2)

  
end subroutine mode_evolution_analysis4

subroutine mode_evolution_analysis5(t,s,m,n,file_unit) 
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  use FFTW3,only: plan_toroidal, in1,out1
  use constants,only:p_
    use control_parameters, only : nh
  implicit none
  include 'fftw3.f03'
  real(p_),intent(in):: t
  integer,intent(in):: m,n
  real(p_),intent(in):: s(0:m-1,0:n-1)
  integer,intent(in):: file_unit
  integer:: itor, jrad

!  call oned_fourier_transform1(a,a_spectrum,m,n)

  jrad=n/2

  in1(:)=s(:,jrad) !copy in, meanwhile convert real array to complex array
  call fftw_execute_dft(plan_toroidal, in1(:), out1(:))     !Fourier transformation along the first dimension
 
  !itor=1 !select the first harmonic
  write(file_unit,'(200(1pe20.8))') t, (real(out1(itor)),imag(out1(itor)), itor=0,nh)
end subroutine mode_evolution_analysis5


end module report_module
