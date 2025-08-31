module diagnosis_mod
contains

subroutine parallel_current_diagnostic(niter,iter)
use fk_module,only: ni0,charge_i, vn_fk
use gk_module,only: ngk0
!use perturbation_field,only: jpar_i_old, jpar_e_old

integer,intent(in):: niter,iter

!if(iter.eq.niter) write(*,*) iter, jpar_e_old(5,5)*vn_gk*ngk0*(-charge_i),jpar_i_old(5,5)*vn_fk*ni0*charge_i,&
!     & jpar_e_old(5,5)*vn_gk*ngk0*(-charge_i)/(jpar_i_old(5,5)*vn_fk*ni0*charge_i)

 !write(*,*) jpar_e_old(15,15)*vn_gk*ngk0*(-charge_i),jpar_i_old(15,15)*vn_fk*ni0*charge_i,&
  !   & jpar_e_old(15,15)*vn_gk*ngk0*(-charge_i)/(jpar_i_old(15,15)*vn_fk*ni0*charge_i)

end subroutine parallel_current_diagnostic

  
subroutine check_rms_ion_weight(w_i_star,nmarker_i)
  use constants,only: p_
  use domain_decomposition,only: myid,numprocs
  use mpi
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: w_i_star(nmarker_i)
  real(p_):: my_sum,all_sum
  integer:: ierr
  all_sum=0.
  my_sum=sum(w_i_star(1:nmarker_i)**2)/nmarker_i
  call MPI_Reduce(my_sum, all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr);
  if(myid.eq.0) write(*,*) 'rms_of_weight_averaged_over_all_particles=',sqrt(all_sum)/numprocs !assume that number of parcles in each process is equal to each other
end subroutine check_rms_ion_weight

subroutine compute_heat_flux(time,ns,nm, touch_bdry_gc, mu_gk, vpar_gk, w_gk,&
     &  xgc, zgc, xdrift1)
  use constants, only : p_, two, kev, three
  use mpi
  use magnetic_coordinates,only: radial_width, dradcor, dtheta, mpol, toroidal_range, xgrid
  use magnetic_coordinates,only: jacobian, vol, grad_psi
  use func_in_mc, only: b_mc_func, grad_psi_func
  use gk_module,only: mass_gk, charge_gk, vn_gk, w_unit, nsm
  use gk_profile_funcs, only : gkt_func, gkn_func, gkdtdx_func
  use radial_module, only : baxis, minor_a
  use domain_decomposition,only: myid
  implicit none
  real(p_), intent(in) :: time
  integer, intent(in) :: ns, nm
  logical, intent(in) :: touch_bdry_gc(:)
  real(p_), intent(in) :: mu_gk(:), vpar_gk(:), w_gk(:), xgc(:), zgc(:), xdrift1(:)
  integer, parameter :: nsub=20
  real(p_) :: tmp_xa(nsub), dv(nsub), myheat_flux(nsub), heat_flux(nsub)
  real(p_) :: x, dx, kinetic, gradient
  real(p_) :: diffusivity(nsub), gyro_bohm(nsub)
  integer :: i,j,k, jeq,ierr
  character(len=64) :: fn
  logical, save :: is_first = .true.
  integer, allocatable, save :: u(:)

  if ((is_first .eqv. .true.) .and. (myid.eq.0)) then 
     is_first=.false.
     allocate(u(nsm))
     do i = 1, nsm
        fn = 'heat_flux_nsx.txt'
        write(fn(13:13),'(i1.1)') i
        open(newunit=u(i), file=fn)
     enddo
  endif

  dx = radial_width/nsub
  do j=1,nsub
     tmp_xa(j) = xgrid(1)+dx*(j-1)
  enddo

  myheat_flux=0._p_
  do k=1,nm
     if(touch_bdry_gc(k) .eqv. .true.) cycle
     j = floor((xgc(k)- tmp_xa(1))/dx)+1
     kinetic = vpar_gk(k)**2/two + mu_gk(k)*b_mc_func(zgc(k), xgc(k))
     myheat_flux(j)=myheat_flux(j) + w_gk(k)*kinetic*xdrift1(k)/grad_psi_func(zgc(k), xgc(k)) !corrected
  enddo

  call MPI_Reduce(myheat_flux, heat_flux, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr)

  if(myid==0) then
     dv(:)=0._p_
     do j=1,nsub
        jeq=1+(j-1)*int(dx/dradcor)
        x =  tmp_xa(j)
        do i=1,mpol-1
           dv(j)=dv(j)+abs(jacobian(i,jeq))*dx*dtheta*toroidal_range
        enddo
        heat_flux(j) = heat_flux(j)*w_unit*mass_gk(ns)*vn_gk(ns)**3/dv(j) !in SI unit
        gyro_bohm(j)=sqrt(mass_gk(ns))*sqrt(gkt_func(x,ns)*kev)**3/(minor_a*Baxis**2*charge_gk(ns)**2) !in SI unit: m^2/s
         gradient = -gkn_func(x,ns)*gkdtdx_func(x,ns)*kev*sum(grad_psi(:,jeq))/mpol
         diffusivity(j) = 2./3.0*heat_flux(j)/gradient/gyro_bohm(j)

     enddo

     write(u(ns),'(50(1pe20.8))') time, sum(diffusivity)/nsub, &
          & sum(heat_flux)/nsub, (heat_flux(j), j=1,nsub)
  endif

end subroutine compute_heat_flux




subroutine draw_grids_on_theta_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !on theta=constant surface, the subroutine name is wrong, not necessarily top view
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nrad
  real(p_),intent(in)::tor_shift_mc(mpol,nrad),r_mc(mpol,nrad),z_mc(mpol,nrad)
  real(p_):: phi,alpha,dalpha
  integer:: i,j,itor,mtor

  i=10
  i=20
  i=60
  i=1
  !  i=mpol
  mtor=20
  dalpha=twopi/mtor

  open(113,file='grids_on_theta_isosurface.txt')
  do itor=1,mtor+1
     !  do itor=1,2
     alpha=dalpha*(itor-1)
     !     do j=1,nrad,10
     !    do j=1,nrad,1
     do j=1,1
        phi=alpha+tor_shift_mc(i,j) !phi is changing due to the radial dependene of tor_shift
        write(113,*) phi,r_mc(i,j),z_mc(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

  open(113,file='grids_on_theta_isosurface2.txt')
  do j=1,nrad,10
     do itor=1,mtor+1
        alpha=dalpha*(itor-1)
        phi=alpha+tor_shift_mc(i,j)
        write(113,*) phi,r_mc(i,j),z_mc(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)
end subroutine draw_grids_on_theta_isosurface


subroutine draw_alpha_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !on a nonzero radial range
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nrad
  real(p_),intent(in)::tor_shift_mc(mpol,nrad),r_mc(mpol,nrad),z_mc(mpol,nrad)
  real(p_):: phi,alpha0
  integer:: i,j

  alpha0=twopi/8._p_
  open(113,file='alpha_isosurface.txt')
  do j=1,101,5
     do i=1,mpol
        phi=alpha0+tor_shift_mc(i,j) 
        !phi=alpha0
        write(113,*) phi,r_mc(i,j),z_mc(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

end subroutine draw_alpha_isosurface


subroutine draw_alpha_contours_on_a_magnetic_surface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !field lines on a magnetic surface
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nrad
  real(p_),intent(in)::tor_shift_mc(mpol,nrad),r_mc(mpol,nrad),z_mc(mpol,nrad)
  real(p_):: phi,alpha0,tor_range
  integer:: i,j,ialpha,ishift,u,iphi
  integer,parameter::nalpha=10,nphi=20

  open(newunit=u,file='alpha_contours_on_magnetic_surface.txt')
  !do j=1,30,2
  j=nrad/2 !select a radial location, i.e., a magnetic surface 
  tor_range=twopi/4
  do ialpha=1,nalpha
     !alpha0=0._p_+tor_range/(nalpha-1)*(ialpha-1)
     alpha0=0._p_
     do i=1,mpol
        phi=alpha0+tor_shift_mc(i,j)
!!$        if(phi<0. .or. phi>tor_range) then !shift into the range [0:tor_range]
!!$           ishift=floor(phi/tor_range)
!!$           phi=phi-ishift*tor_range 
!!$           alpha0=alpha0-ishift*tor_range
!!$           write(u,*)
!!$           write(u,*)
!!$        endif
        !phi=alpha0
        write(u,*) phi,z_mc(i,j),r_mc(i,j)
     enddo
     write(u,*)
     write(u,*)
  enddo
  close(u)
  open(newunit=u,file='ref_surface')
  do iphi=1,nphi
     phi=0.+twopi/(nphi-1)*(iphi-1)
     do i=1,mpol
        write(u,*) phi,z_mc(i,j),r_mc(i,j)
     enddo
     write(u,*)
     write(u,*)
  enddo
  close(u)
phi=alpha0+tor_shift_mc(1,j)

call field_line_tracing_simplified(r_mc(1,j),z_mc(1,j),phi)
end subroutine draw_alpha_contours_on_a_magnetic_surface


subroutine calculate_possibility_density(v,total_number,sample_number,starting_value,ending_value)
  use constants,only:p_
  use  constants,only:one,two,four,five,twopi,eight
  implicit none

  integer,intent(in):: total_number,sample_number
  real(p_),intent(in):: v(total_number)
  real(p_),intent(in):: starting_value,ending_value
  real(p_):: xcenter(sample_number-1)
  real(p_):: possibility_density(sample_number-1)
  real(p_):: bin_npt(sample_number-1)
  real(p_):: interval
  real(p_):: x(sample_number)
  integer:: i,j

  interval=(ending_value-starting_value)/(sample_number-1)
  do i=1,sample_number
     x(i)=starting_value+interval*(i-1)
  enddo

  bin_npt=0
  do i=1,sample_number-1
     do j=1,total_number
        if (v(j).ge.x(i) .and. v(j).le.x(i+1)) bin_npt(i)=bin_npt(i)+1
     enddo
  enddo

  do i=1,sample_number-1
     xcenter(i)=(x(i)+x(i+1))/two
     possibility_density(i)=bin_npt(i)/(total_number)/interval
     write(*,*) xcenter(i), possibility_density(i)
  enddo

end subroutine calculate_possibility_density


subroutine report(t)
  use magnetic_coordinates,only: mtor,nrad
!  use perturbation_field, only:ex=>ex_left,ey=>ey_left,epar=>epar_left

!  use perturbation_field,only:source_e1,source_e2 !,jper_x_i_left,jper_y_i_left
!  use perturbation_field,only:source1,source2,source3,source_faraday1, source_faraday2,source_faraday3
  use constants,only:p_

  real(p_),intent(in):: t
  integer:: i,j

i=mtor/2
j=nrad/2

!write(*,'(7(1pe14.4))') t,ex(i,j),ey(i,j),epar(i,j),mf_x(i,j),mf_y(i,j),mf_par(i,j)
!write(*,*) mf_x(i,j),mf_y(i,j),mf_par(i,j)
!write(*,*) jper_x_i_left(mtor/2,nrad/2),jper_y_i_left(mtor/2,nrad/2)!,source_e1(mtor/2,nrad/2),source_e2(mtor/2,nrad/2)
!write(*,*) source1(mtor/3,nrad/3),source2(mtor/3,nrad/3),source3(mtor/3,nrad/3)
!write(*,*) source_faraday1(mtor/3,nrad/3),source_faraday2(mtor/3,nrad/3),source_faraday3(mtor/3,nrad/3)
end subroutine report


subroutine mode_evolution_analysis(t)
  use constants,only: one
  use constants,only:p_
  use magnetic_coordinates,only: dtor,dradcor,m=>mtor,n=>nrad
!  use perturbation_field,only: ex_left
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


subroutine mode_evolution_analysis2(t) 
  use constants,only: one
  use constants,only:p_
  use magnetic_coordinates,only: m=>mtor,n=>nrad
  use perturbation_field,only: ef_cyl_phi_left,ef_cyl_r_left,ef_cyl_z_left
  use perturbation_field,only: ef_cyl_phi_right,ef_cyl_r_right,ef_cyl_z_right
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


subroutine mode_evolution_analysis4(t, s, m, n, file_unit) 
  use constants,only: one, p_
  use transform_module
  use control_parameters, only : nh_min
  !  use fourn_module,only: twod_fourier_transform_nr
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  real(p_), intent(in) :: t
  integer, intent(in) :: file_unit
  complex(p_) :: spectrum(0:m-1, 0:n-1)
  integer :: ntor, kr

  !  call oned_sine_transform2(s,a_dst,m,n) 
  !  call oned_fourier_transform1(a_dst,a_spectrum,m,n) 

  call dst_dft(s, spectrum, m, n) !using DST for radial direction and DFT for toroidal direction
  spectrum = spectrum/(2*(n+1)*m) !the expansion coeficient is the dst_dft devided by (n+1)*m, see my notes on Fourier analysis

  ntor = nh_min
  write(file_unit,'(20ES18.4)') t, (real(spectrum(ntor,kr)), imag(spectrum(ntor,kr)), kr=0,3)

end subroutine mode_evolution_analysis4

subroutine mode_evolution_analysis5(t, s, m, n, file_unit)
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  use FFTW3, only: plan_toroidal, in1, out1
  use control_parameters, only : nh_min, nh_max
  implicit none
  include 'fftw3.f03'
  real(p_), intent(in) :: t
  integer, intent(in) :: m,n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  integer, intent(in) :: file_unit
  integer :: i, jrad

!  call oned_fourier_transform1(s,s_spectrum,m,n)                          

  jrad=154

  in1(:) = s(:,jrad) !copy in, meanwhile convert real array to complex array                                      
  call fftw_execute_dft(plan_toroidal, in1(:), out1(:))

  write(file_unit,'(2000(1pe20.8))') t, (real(out1(i))/m, imag(out1(i))/m, i=nh_min,nh_max)
end subroutine mode_evolution_analysis5


subroutine mode_evolution_analysis6(t, s, m, n, file_unit) 
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  use FFTW3, only: plan_toroidal, in1, out1
  use control_parameters, only : nh_min, nh_max
  implicit none
  include 'fftw3.f03'
  real(p_), intent(in) :: t
  integer, intent(in) :: m,n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  integer, intent(in) :: file_unit
  integer, parameter ::   jw = 5
  integer :: i, j, j0, jlow, jupp, nr
  complex(p_), allocatable :: spectrum(:,:)

  j0 = n/2
  jlow = j0-jw
  jupp = j0+jw
  nr = jupp - jlow +1
  allocate(spectrum(0:m-1, jlow:jupp))

  do j = jlow, jupp
     in1(:) = s(:,j) !copy in, meanwhile convert real array to complex array
     call fftw_execute_dft(plan_toroidal, in1(:), out1(:))
     spectrum(:,j) = out1(:)
  enddo
  
  i = nh_min
  write(file_unit,'(2000(1pe18.4))') t, (real(spectrum(i,j))/m, imag(spectrum(i,j))/m, j=jlow,jupp)

end subroutine mode_evolution_analysis6


end module diagnosis_mod


module mode_structure
  implicit none
contains
  subroutine mode_structure_on_xy_plane(kt,GCLR,a,partial_file_name)
    use constants,only:p_
    use magnetic_coordinates,only:xgrid,ygrid
!    use domain_decomposition,only:myid
    integer,intent(in)::kt,GCLR
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in)::partial_file_name
    character(100):: full_file_name
    integer:: i,j,m,n,u

    full_file_name='ms/xy_polxxx_txxxxxx'//partial_file_name
    write(full_file_name(10:12),'(i3.3)') GCLR
    write(full_file_name(15:20),'(i6.6)') kt

    open(newunit=u,file=full_file_name)
    m=size(a,1)
    n=size(a,2)
    do j=1,n
       do i=1,m
          write(u,*) xgrid(j),ygrid(i),a(i,j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine mode_structure_on_xy_plane

  subroutine mode_structure_on_xz_plane(kt,a,partial_file_name)
    use constants,only:p_
    use constants,only:pi
    use magnetic_coordinates,only: xgrid,zgrid,dtheta,mtor
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,GCLR_cut,dtheta2
    use constants,only:twopi
    use mpi
    implicit none
    integer,intent(in)::kt
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in)::partial_file_name
    integer:: itor,ipol,j,ierr,m,n 
    real(p_):: my_a_xz_plane(size(a,2)),a_xz_plane(size(a,2),0:numprocs/ntube-1)
    real(p_)::my_theta
    character(100)::full_file_name
    integer:: u !file unit number

    m=size(a,1)
    n=size(a,2)

    itor=mtor/2 !choose a alpha (i.e., y) grid
    do j=1,n
       my_a_xz_plane(j)=a(itor,j)
    enddo
    call MPI_gather(my_a_xz_plane, n, MPI_real8, &
                  & a_xz_plane,    n, MPI_real8, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       full_file_name='ms/xz_txxxxxx'//partial_file_name
       write(full_file_name(8:13),'(i6.6)') kt
       open(newunit=u,file=full_file_name)
       do ipol=0,numprocs/ntube-1 !poloidal direction
          my_theta=-pi+ipol*dtheta2
          !if(ipol.gt.GCLR_cut) my_theta=-pi+(ipol-GCLR_cut-1)*dtheta2
          do j=1,n !radial direction
             write(u,*) xgrid(j), my_theta,a_xz_plane(j,ipol)
          enddo
          write(u,*)
          !if(ipol.eq.GCLR_cut) write(u,*) !to inform gnuplot that this is a new data block, to prevent gnuplot using line connection between the following data and previous data
       enddo
       close(u)
    endif
  end subroutine mode_structure_on_xz_plane


  subroutine mode_structure_on_yz_plane(kt,a,partial_file_name)
    use constants,only: p_, pi, twopi
    use magnetic_coordinates,only: mtor,ygrid,zgrid, dtheta,nrad, mpol2,&
         & r_mc, z_mc, tor_shift_mc
    use domain_decomposition,only: myid,tube_comm,dtheta2, multi_eq_cells
    use mpi
    implicit none
    integer,intent(in)::kt
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in)::partial_file_name
    integer:: itor,ipol,ipol_eq,jrad,ierr,m,n
    real(p_):: my_a_yz_plane(size(a,1)),a_yz_plane(size(a,1),0:mpol2)
    real(p_):: phi
    character(100)::full_file_name
    integer:: u !file unit number
    m=size(a,1)
    n=size(a,2)
    jrad=nrad/2 !choose a radial index
    do itor=1,mtor
       my_a_yz_plane(itor)=a(itor,jrad)
    enddo
    call MPI_gather(my_a_yz_plane, m, MPI_real8,&
         &          a_yz_plane,    m, MPI_real8, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       full_file_name='ms/yz_txxxxxx'//partial_file_name
       write(full_file_name(8:13),'(i6.6)') kt
       open(newunit=u,file=full_file_name)
       do ipol=0,mpol2 !poloidal direction
          ipol_eq=multi_eq_cells*ipol + 1 !index in the equilibrium grids
          do itor=1,mtor !toroidal direction
             phi=ygrid(itor)+tor_shift_mc(ipol_eq,jrad)
             write(u,'(19ES16.5E3)') ygrid(itor), zgrid(ipol_eq), a_yz_plane(itor,ipol), &
                  & r_mc(ipol_eq,jrad), z_mc(ipol_eq,jrad), phi
          enddo
          write(u,*) 
          !if(ipol.eq.GCLR_cut) write(u,*) !to inform gnuplot that this is a new data block, to prevent gnuplot using line connection between the following data and previous data
       enddo
       close(u)
    endif
  end subroutine mode_structure_on_yz_plane


  subroutine mode_structure_on_poloidal_plane(kt, a, str)
    use constants,only : p_, one
    use transform_module,only:oned_fourier_transform1, oned_backward_fourier_transform1
    use domain_decomposition,only:numprocs,myid,ntube, theta_start, dtheta2
    implicit none
    integer, intent(in) :: kt
    character(*), intent(in) :: str
    real(p_), intent(in) :: a(:,:, :)
    integer, parameter :: nz = 40 !along field line
    real(p_) :: a0(size(a,1), size(a,2), 0:nz-1), a1(size(a,1),size(a,2), 0:nz-1)
    real(p_) :: theta(0:nz-1), c
    complex(p_) :: a_dft(size(a,1), size(a,2), 0:nz-1)
    real(p_) :: perturb(size(a,1), size(a,2), 0:nz-1)
    character(len=100) :: file_name
    integer :: m, n, j, i

    m=size(a,1) !toroidal
    n=size(a,2) !radial

    do i = 0, nz-1 !interpolate to more z gridpoints 
       theta(i) = theta_start + dtheta2/nz*i
       c = (theta(i) - theta_start)/dtheta2
       perturb(:,:,i) = a(:,:,1)*(one-c) + a(:,:,2)*c
    enddo

    do i = 0, nz-1
       call oned_fourier_transform1(perturb(:,:,i), a_dft(:,:,i), m,n) !1d DFT long the first dimension
       do j=1,n
          a0(:,j, i) = real(a_dft(1,j, i)) !space distribution of the n=0 harmonic
       enddo
    enddo
    file_name='ms/poloidal_plane_txxxxxx'//str//'_zonal'
    write(file_name(20:25),'(i6.6)') kt
    call mode_structure_on_poloidal_plane0(a0, m, n, nz, theta, file_name)

    a_dft(1,:,:) = 0 !remove the n=0 component for ploting, this does not affect the simulation
    do i = 0, nz-1
       call oned_backward_fourier_transform1(a_dft(:,:, i), a1(:,:,i), m,n)
    enddo
    
    file_name='ms/poloidal_plane_txxxxxx'//str//'_nonzonal'
    write(file_name(20:25),'(i6.6)') kt
    call mode_structure_on_poloidal_plane0(a1, m, n, nz, theta, file_name)

  end subroutine mode_structure_on_poloidal_plane

  subroutine mode_structure_on_poloidal_plane0(yxz, ny, nx, nz, theta, file_name)
    use constants, only : p_, pi, twopi, one
    use magnetic_coordinates, only: ygrid, r_mc, z_mc, tor_shift_mc, mpol, nsegment, &
         & dtheta, zgrid, pfn, toroidal_range
    use domain_decomposition, only : numprocs,myid,ntube,tube_comm, theta_start, dtheta2, ipol_eq
    use interpolate_module
    use math,only: shift_toroidal
    use mpi
    implicit none
    integer, intent(in) :: ny, nx, nz
    real(p_),intent(in) :: yxz(ny, nx, 0:nz-1), theta(0:nz-1)
    character(len=*), intent(in) :: file_name
    real(p_), allocatable :: field(:,:), field0(:,:)
    real(p_) :: phi, my_theta, alpha, r, z, th, tor_shift
    integer :: j, ierr, i, jeq, u, npol, iz

    allocate(field(nx, 0:nz-1))
    npol = nz*numprocs/ntube
    allocate(field0(nx, 0:npol-1))

    !choose a cylindrical toroidal angle, for which the mode structure is computed
    phi = 0._p_ + 0.5_p_*twopi/nsegment

    do j = 1, nx
       jeq = j
       do i = 0, nz-1
          call linear_1d_interpolate(mpol, zgrid, tor_shift_mc(:,jeq), theta(i), tor_shift)
          alpha = phi - tor_shift
          call shift_toroidal(alpha, toroidal_range)
          !to the same cylindrical toroidal angle:
          call linear_1d_interpolate(ny, ygrid, yxz(:,j,i), alpha, field(j, i) ) 
       enddo
    enddo

    call MPI_gather(field, nx*nz, MPI_real8, field0, nx*nz, MPI_real8, 0, tube_COMM, ierr)

    if(myid.eq.0) then
       open(newunit=u,file=file_name)
       do j=1,nx
          jeq = j
          do i = 0, numprocs/ntube-1
             my_theta = -pi + i*dtheta2
             do iz = 0, nz-1
                th = my_theta + dtheta2/nz*iz
                call linear_1d_interpolate(mpol, zgrid, r_mc(:,jeq), th, r)
                call linear_1d_interpolate(mpol, zgrid, z_mc(:,jeq), th, z)
                write(u,'(19ES16.5E3)') r, z, field0(j, i*nz+iz), pfn(jeq), th
             enddo
          enddo
          write(u,*)
          write(u,*) 
       enddo
       close(u)
    endif
  end subroutine mode_structure_on_poloidal_plane0

end module mode_structure


module spectrum_diagnostic
  use constants,only:p_
  implicit none
  save
  real(p_):: spectrum_amplitude_old
contains
  subroutine spectrum_diagnostic_routine(iter,t,a,file_unit)
    use constants,only:pi
    use magnetic_coordinates,only: xgrid,zgrid,dtheta,mtor
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,GCLR_cut
    use transform_module,only: twod_fourier_transform_xz,oned_fourier_transform1
    use constants,only:twopi
    use mpi
    integer,intent(in):: iter
    real(p_),intent(in):: a(:,:),t
    integer,intent(in):: file_unit
    complex(p_):: a_dft(size(a,1),size(a,2))
    complex(p_):: my_a_xz_plane(size(a,2)),a_xz_plane(size(a,2),numprocs/ntube)
    complex(p_):: a_xz_plane_fft(size(a,2),numprocs/ntube)
    integer:: itor,ipol,j,ierr,m,n,mz
    real(p_):: spectrum_amplitude

    m=size(a,1)
    n=size(a,2)
    mz=numprocs/ntube
!!$    itor=mtor/2 !choose a y grid
!!$    do j=1,n
!!$       my_a_xz_plane(j)=a(itor,j)
!!$    enddo
    call oned_fourier_transform1(a,a_dft,m,n) !Fourier transform along toroidal direction
    do j=1,n
       my_a_xz_plane(j)=a_dft(2,j) !select the fundamental harmonic of the toroidal expansion
    enddo

    call MPI_gather(my_a_xz_plane, n, MPI_complex8, &
         & a_xz_plane,    n, MPI_complex8, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       call twod_fourier_transform_xz(a_xz_plane,a_xz_plane_fft,n,mz)
       write(file_unit,'(20(1pe20.8))') t,(real(a_xz_plane(j,j)),imag(a_xz_plane(j,j)),j=1,3)
!!$       spectrum_amplitude=abs(a_xz_plane_fft(4,4))
!!$       write(*,*) iter,'relative_error_in_low_kx_low_kz_spectrum_amplitude=',&
!!$            & abs(spectrum_amplitude_old-spectrum_amplitude)/spectrum_amplitude
!!$       spectrum_amplitude_old=spectrum_amplitude !prepare for the next convergence checking
    endif

  end subroutine spectrum_diagnostic_routine
end module spectrum_diagnostic


subroutine visualize_grid()
  use constants, only : p_, twopi, zero
  use magnetic_coordinates, only : r_mc, z_mc, mpol,nrad, mtor, ygrid, tor_shift_mc, toroidal_range
  use magnetic_field, only : psi_func, qfunc0
  use math, only : shift_toroidal
  implicit none
  integer, parameter :: npt=2000
  integer :: u,i,j,k
  real(p_) :: phi, phi0
  real(p_) :: rf(npt), zf(npt), phif(npt), dphi, qval

  open(newunit=u,file='grid.txt')
  do k=1, mtor
     do i=1, mpol
        do j=1,nrad
           phi = ygrid(k) + tor_shift_mc(i,j)
           !call shift_toroidal(phi,toroidal_range) 
           write(u,*) r_mc(i,j), z_mc(i,j), phi
        enddo
     enddo
  enddo
  close(u)

  i = (mpol+1)/2
  j = nrad/2
  !phi0 = ygrid(1) + tor_shift_mc(i,j)
  phi0 = 0
  qval = qfunc0(psi_func(r_mc(i,j), z_mc(i,j)))
  dphi = twopi/2*abs(qval)/(npt-1)
  call field_line_tracing0(r_mc(i,j), z_mc(i,j), phi0, npt, dphi, rf,zf,phif)
  open(newunit=u,file='grid_z.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  call field_line_tracing0(r_mc(i,j),z_mc(i,j), phi0, npt, -dphi, rf,zf,phif)
  open(newunit=u,file='grid_z2.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)

  !-------------------------------
  j = 1
  phi0=ygrid(1) + tor_shift_mc(i,j)
  qval = qfunc0(psi_func(r_mc(i,j), z_mc(i,j)))
  dphi = twopi*abs(qval)/(npt-1)/2
  call field_line_tracing0(r_mc(i,j), z_mc(i,j), phi0, npt, dphi, rf,zf,phif)
  open(newunit=u,file='grid_z3.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  call field_line_tracing0(r_mc(i,j),z_mc(i,j), phi0, npt, -dphi, rf,zf,phif)
  open(newunit=u,file='grid_z4.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  !-------------------------------
  j = nrad
  phi0=ygrid(1) + tor_shift_mc(i,j)
  qval = qfunc0(psi_func(r_mc(i,j), z_mc(i,j)))
  dphi = twopi*abs(qval)/(npt-1)/2
  call field_line_tracing0(r_mc(i,j), z_mc(i,j), phi0, npt, dphi, rf,zf,phif)
  open(newunit=u,file='grid_z5.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  call field_line_tracing0(r_mc(i,j),z_mc(i,j), phi0, npt, -dphi, rf,zf,phif)
  open(newunit=u,file='grid_z6.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  
  open(newunit=u,file='inner_bdry.txt')
  do i =1, mpol
     write(u,*) r_mc(i, 1), z_mc(i, 1)
  enddo
  close(u)

  open(newunit=u,file='outer_bdry.txt')
  do i =1, mpol
     write(u,*)  r_mc(i, nrad), z_mc(i, nrad)
  enddo
  close(u)

end subroutine visualize_grid


