module diagnosis_mod
contains

  subroutine compute_particle_and_heat_flux(time, ns, lost, mu_gk, vpar_gk, w_gk,&
       &  xgc, zgc, xdrift1, xdrift0)
    use constants, only: p_, two, kev, three
    use magnetic_coordinates, only: radial_width, dradcor, dtheta, mpol, toroidal_range, xgrid, &
         & jacobian, vol, grad_psi
    use func_in_mc, only: b_mc_func, grad_psi_func
    use gk_module,only: mass_gk, charge_gk, nm_gk, vn_gk, w_unit, nsm
    use gk_profile_funcs, only : gkt_func, gkn_func, gkdtdx_func
    use radial_module, only : baxis, minor_a
    use domain_decomposition, only: myid
    use mpi
    implicit none
    real(p_), intent(in) :: time
    integer, intent(in) :: ns
    logical, intent(in) :: lost(:)
    real(p_), intent(in) :: mu_gk(:), vpar_gk(:), w_gk(:), xgc(:), zgc(:), xdrift1(:), xdrift0(:)
    integer, parameter :: nsub = 50
    real(p_) :: tmp_xa(nsub), dv(nsub)
    real(p_) :: myheat_flux(nsub), heat_flux(nsub), myheat_flux0(nsub), heat_flux0(nsub)
    real(p_) :: myptcl_flux(nsub), ptcl_flux(nsub), mydensity_perturb(nsub), density_perturb(nsub)
    real(p_) :: myptcl_flux0(nsub), ptcl_flux0(nsub)
    real(p_) :: diffusivity(nsub), gyro_bohm(nsub)
    real(p_) :: x, dx, kinetic, gradient, trash
    integer :: nm, i, j, k, jeq, ierr
    character(len=64) :: fn1, fn2,  fn1b, fn2b, fn3
    logical, save :: is_first = .true.
    integer, allocatable, save :: u1(:), u2(:), u1b(:), u2b(:), u3(:)

    nm = nm_gk(ns)

    if ((is_first .eqv. .true.) .and. (myid==0)) then 
       is_first=.false.
       allocate(u1(nsm), u1b(nsm))
       allocate(u2(nsm), u2b(nsm))
       allocate(u3(nsm))       
       do i = 1, nsm
          fn1 = 'heat_flux_nsx.txt'; fn1b = 'heat_flux0_nsx.txt'
          fn2 = 'ptcl_flux_nsx.txt'; fn2b = 'ptcl_flux0_nsx.txt'
          fn3 = 'dens_pert_nsx.txt'
          write(fn1(13:13),'(i1.1)') i
          write(fn2(13:13),'(i1.1)') i
          write(fn1b(14:14),'(i1.1)') i
          write(fn2b(14:14),'(i1.1)') i
          write(fn3(13:13),'(i1.1)') i
          open(newunit=u1(i), file=fn1, position="append")
          open(newunit=u2(i), file=fn2, position="append")
          open(newunit=u1b(i), file=fn1b, position="append")
          open(newunit=u2b(i), file=fn2b, position="append")
          open(newunit=u3(i), file=fn3, position="append")
       enddo
    endif

    dx = radial_width/nsub
    do j = 1, nsub
       tmp_xa(j) = xgrid(1) + dx*(j-1)
    enddo

    myheat_flux = 0._p_; myheat_flux0 = 0._p_
    myptcl_flux = 0._p_; myptcl_flux0 = 0._p_
    mydensity_perturb = 0._p_
    do k = 1, nm
       if(lost(k) .eqv. .true.) cycle
       j = floor((xgc(k)- tmp_xa(1))/dx)+1
       j = min(j, nsub)
       j = max(j, 1)
       kinetic = vpar_gk(k)**2/two + mu_gk(k)*b_mc_func(zgc(k), xgc(k))
       myheat_flux(j) = myheat_flux(j) + w_gk(k)*kinetic*xdrift1(k)/grad_psi_func(zgc(k), xgc(k))
       myptcl_flux(j) = myptcl_flux(j) + w_gk(k)        *xdrift1(k)/grad_psi_func(zgc(k), xgc(k))
       myheat_flux0(j) = myheat_flux0(j) + w_gk(k)*kinetic*xdrift0(k)/grad_psi_func(zgc(k), xgc(k))
       myptcl_flux0(j) = myptcl_flux0(j) + w_gk(k)        *xdrift0(k)/grad_psi_func(zgc(k), xgc(k))
       mydensity_perturb(j) = mydensity_perturb(j) + w_gk(k)
    enddo

    call MPI_Reduce(myheat_flux, heat_flux, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(myptcl_flux, ptcl_flux, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(myheat_flux0, heat_flux0, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(myptcl_flux0, ptcl_flux0, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(mydensity_perturb, density_perturb, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)    

    if(myid==0) then
       do j = 1, nsub
          jeq = 1 + (j-1)*int(dx/dradcor)
          x =  tmp_xa(j)
          dv(j) = sum(abs(jacobian(1:mpol-1, jeq)))*dx*dtheta*toroidal_range
          ptcl_flux(j) = ptcl_flux(j)*w_unit            *vn_gk(ns)   /dv(j) 
          heat_flux(j) = heat_flux(j)*w_unit*mass_gk(ns)*vn_gk(ns)**3/dv(j) !in unit m^2/s
          ptcl_flux0(j) = ptcl_flux0(j)*w_unit            *vn_gk(ns)   /dv(j) 
          heat_flux0(j) = heat_flux0(j)*w_unit*mass_gk(ns)*vn_gk(ns)**3/dv(j) !in unit m^2/s

          density_perturb(j) = density_perturb(j)*w_unit/dv(j) 
          gyro_bohm(j) = sqrt(mass_gk(ns))*sqrt(gkt_func(x,ns)*kev)**3/(minor_a*Baxis**2*charge_gk(ns)**2) !in unit m^2/s
          gradient = -gkn_func(x,ns)*gkdtdx_func(x,ns)*kev*sum(grad_psi(:,jeq))/mpol !assume constant density
          diffusivity(j) = 2./3.0*heat_flux(j)/gradient
       enddo
       !write(u1(ns),'(999(1pe20.8))') time, sum(diffusivity)/nsub, sum(diffusivity/gyro_bohm)/nsub, &
       write(u1(ns),'(999(1pe20.8))')  time, sum(heat_flux)/nsub, (heat_flux(j), j=1, nsub)
       write(u2(ns),'(999(1pe20.8))')  time, sum(ptcl_flux)/nsub, (ptcl_flux(j), j=1, nsub)
       write(u1b(ns),'(999(1pe20.8))') time, sum(heat_flux0)/nsub, (heat_flux0(j), j=1, nsub)
       write(u2b(ns),'(999(1pe20.8))') time, sum(ptcl_flux0)/nsub, (ptcl_flux0(j), j=1, nsub)
       write(u3(ns), '(999(1pe20.8))') time, sum(density_perturb)/nsub, (density_perturb(j), j=1, nsub)
    endif

  end subroutine compute_particle_and_heat_flux


  subroutine compute_entropy(time, ns, lost, w_gk, xgc, v)
    use constants, only: p_, two
    use gk_module,only: nm_gk, w_unit, nsm, ps_vol_gk
    use load_gk_mod, only: f0
    use domain_decomposition, only: myid, numprocs
    use mpi
    implicit none
    real(p_), intent(in) :: time
    integer, intent(in) :: ns
    logical, intent(in) :: lost(:)
    real(p_), intent(in) :: w_gk(:), xgc(:), v(:)
    real(p_) :: my_entropy, entropy, my_f1f0, f1f0, tmp
    integer :: nm, i, k, ierr, nptcl
    character(len=64) :: fn1, fn2
    logical, save :: is_first = .true.
    integer, allocatable, save :: u1(:), u2(:)

    if ((is_first .eqv. .true.) .and. (myid==0)) then 
       is_first=.false.
       allocate(u1(nsm), u2(nsm))
       do i = 1, nsm
          fn1 = 'entropy_nsx.txt'
          fn2 = 'f1f0_nsx.txt'
          write(fn1(11:11),'(i1.1)') i
          write(fn2(8:8),'(i1.1)') i
          open(newunit=u1(i), file=fn1, position="append")
          open(newunit=u2(i), file=fn2, position="append")
       enddo
    endif

    nm = nm_gk(ns)
    
    my_entropy = 0.
    my_f1f0 = 0.
    nptcl = 0
    do k = 1, nm
       if(lost(k) .eqv. .true.) cycle
       nptcl = nptcl +1
       tmp = w_gk(k)/(f0(xgc(k), v(k), ns)*ps_vol_gk(k,ns)/w_unit) ! delta_f/f0
       my_entropy = my_entropy + tmp**2
       my_f1f0 = my_f1f0 + tmp
    enddo
    my_entropy = my_entropy / nptcl
    my_f1f0 = my_f1f0 / nptcl
    call MPI_Reduce(my_entropy, entropy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(my_f1f0,    f1f0,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)    
    
    if(myid==0) write(u1(ns),'(50(1pe20.8))') time, entropy/numprocs
    if(myid==0) write(u2(ns),'(50(1pe20.8))') time, f1f0/numprocs
  end subroutine compute_entropy
 
  
subroutine count_lost_markers_gk(ns, radcor)  !diagnosis
  use constants, only : p_
    use gk_module, only: nm_gk
    use domain_decomposition, only: myid
    use magnetic_coordinates, only: xlow, xupp
    implicit none
    integer, intent(in) :: ns
    real(p_), intent(in) :: radcor(:)
    integer :: nlost, k

    nlost = 0
    do k = 1, nm_gk(ns)
       if( (radcor(k) >= xupp) .or.  (radcor(k) <= xlow) ) nlost = nlost + 1
    enddo
    write(*,'(A50, 20i15)') 'myid, species, total_makers, lost_markers =', myid, ns, nm_gk(ns), nlost
  end subroutine count_lost_markers_gk

  
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
!write(*,*) jper_x_i_left(mtor/2,nrad/2),jper_y_i_left(mtor/2,nrad/2)!,source_e1(mtor/2,nrad/2),source_e2(mtor/2,nrad/2)
!write(*,*) source1(mtor/3,nrad/3),source2(mtor/3,nrad/3),source3(mtor/3,nrad/3)
!write(*,*) source_faraday1(mtor/3,nrad/3),source_faraday2(mtor/3,nrad/3),source_faraday3(mtor/3,nrad/3)
end subroutine report


subroutine bperp_perturbation(t, u)
  use constants, only: p_, kev
  use normalizing, only: tu, qu, vu
  use perturbation_field, only: apara
  use magnetic_coordinates, only: nrad, mtor, grad_psi_dot_grad_alpha, grad_psi, grad_alpha
  use domain_decomposition, only: ipol_eq
  use derivatives_in_xyz, only: x_derivative, y_derivative
  implicit none
  real(p_), intent(in) :: t
  integer, intent(in) :: u
  real(p_) :: bperp(nrad), ax(mtor, nrad), ay(mtor, nrad), gx(nrad), gy(nrad), gxdgy(nrad)
  integer :: i
  
  call x_derivative(apara(:,:,1), ax(:,:))
  call y_derivative(apara(:,:,1), ay(:,:))
  ax = ax*tu*kev/(qu*vu) ! change to SI unit
  ay = ay*tu*kev/(qu*vu)
  gx(:) = grad_psi(ipol_eq, :)
  gy(:) = grad_alpha(ipol_eq, :)
  gxdgy(:) = grad_psi_dot_grad_alpha(ipol_eq, :)
  i = mtor/2
  bperp = sqrt(ax(i,:)**2*gx(:)**2 + ay(i,:)**2*gy(:)**2 + 2*ax(i,:)*ay(i,:)*gxdgy(:))
  write(u, '(3000(1pe18.4))') t, bperp(1:nrad:2)
end subroutine bperp_perturbation


subroutine mode_evolution(t, s, m, n, file_unit) 
  use constants, only: p_
  implicit none
  real(p_), intent(in) :: t
  integer, intent(in) :: m, n
  real(p_), intent(in) :: s(0:m-1, n)
  integer, intent(in) :: file_unit
  integer :: i0, j0, j

  i0 = m/2 ! Toroidal location
  !j0 = n/2 ! Radial location
  write(file_unit,'(3000(1pe18.4))') t, (s(i0, j),  j = 1, n, 4)
  !write(file_unit,'(3000(1pe18.4))') t, s(i0, j0)

end subroutine mode_evolution

subroutine mode_evolution2(t) 
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

end subroutine mode_evolution2


subroutine mode_evolution3(t,a,m,n,file_unit) 
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

end subroutine mode_evolution3


subroutine mode_evolution4(t, s, m, n, file_unit) 
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

end subroutine mode_evolution4

subroutine mode_evolution5(t, s, m, n, file_unit)
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  use my_FFTW3, only: plan_toroidal, in1, out1
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

  write(file_unit,'(2000(1pe20.8))') t, (real(out1(i))/m, imag(out1(i))/m, i = nh_min, nh_max)
end subroutine mode_evolution5


subroutine nharmonic_evolution(t, sdft, file_unit) 
  use constants, only: p_
  use magnetic_coordinates, only : nrad, mtor
  use control_parameters, only: nh_min, nh_max
  implicit none
  real(p_), intent(in) :: t
  complex(p_), intent(in) :: sdft(0:,:)
  integer, intent(in) :: file_unit
  integer :: i, j

  !j0 = nrad/2 !choose a radial location
  write(file_unit,'(900(1pe18.4))') t, ((real(sdft(i, j))/mtor, imag(sdft(i, j))/mtor, i = nh_min, nh_max), j = 1, nrad-2, 4)

end subroutine nharmonic_evolution



subroutine mode_evolution7(t, s, m, n, file_unit) 
  use constants, only: p_
  use, intrinsic :: iso_c_binding
  use my_FFTW3, only: plan_toroidal, in1, out1
  use control_parameters, only: nh_min, nh_max
  implicit none
  include 'fftw3.f03'
  real(p_), intent(in) :: t
  integer, intent(in) :: m, n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  integer, intent(in) :: file_unit
  integer :: i, j
  complex(p_) :: spectrum(0:m-1)

  j = n/2
  in1(:) = s(:,j) !copy in, meanwhile convert real array to complex array
  call fftw_execute_dft(plan_toroidal, in1(:), out1(:))
  spectrum(:) = out1(:)/m

  write(file_unit,'(2000(1pe18.4))') t, (real(spectrum(i)), imag(spectrum(i)), i = nh_min, nh_max)

end subroutine mode_evolution7



end module diagnosis_mod


module mode_structure
  implicit none
contains
  subroutine mode_structure_in_xy_plane(kt,GCLR,a,partial_file_name)
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
          write(u,*) xgrid(j), ygrid(i), a(i,j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine mode_structure_in_xy_plane

  subroutine mode_structure_in_xz_plane(kt, a, partial_file_name)
    use constants,only:p_
    use constants,only:pi
    use magnetic_coordinates,only: xgrid,zgrid,dtheta,mtor
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,dtheta2
    use constants,only:twopi
    use mpi
    implicit none
    integer, intent(in) :: kt
    real(p_), intent(in) :: a(:,:)
    character(len=*), intent(in) :: partial_file_name
    integer :: itor, ipol, j, ierr, m, n 
    real(p_) :: my_a_xz_plane(size(a,2)),a_xz_plane(size(a,2),0:numprocs/ntube-1)
    real(p_) :: my_theta
    character(100) :: full_file_name
    integer :: u !file unit number

    m=size(a,1)
    n=size(a,2)

    itor=mtor/2 !choose a alpha (i.e., y) gridpoint
    do j=1,n
       my_a_xz_plane(j)=a(itor,j)
    enddo
    call MPI_gather(my_a_xz_plane, n, MPI_real8, &
                  & a_xz_plane,    n, MPI_real8, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       full_file_name='ms/xz_txxxxxx'//partial_file_name
       write(full_file_name(8:13),'(i6.6)') kt
       open(newunit=u,file=full_file_name)
       do ipol = 0, numprocs/ntube-1 !poloidal direction
          my_theta=-pi+ipol*dtheta2
          do j=1,n !radial direction
             write(u,*) xgrid(j), my_theta, a_xz_plane(j,ipol)
          enddo
          write(u,*)
       enddo
       close(u)
    endif
  end subroutine mode_structure_in_xz_plane


  subroutine mode_structure_in_yz_plane(kt, ayxz, partial_file_name)
    use constants, only: p_, pi, twopi
    use magnetic_coordinates, only: mtor, ygrid, zgrid, nrad, mpol2, &
         & r_mc, z_mc, tor_shift_mc
    use domain_decomposition, only: myid, tube_comm, multi_eq_cells, gclr
    use mpi
    implicit none
    integer, intent(in) :: kt
    real(p_), intent(in) :: ayxz(:, :, :)
    character(len=*), intent(in) :: partial_file_name
    real(p_), allocatable :: ay(:), ayz(:,:)
    real(p_) :: phi
    character(100) :: full_file_name
    integer :: itor, ipol, ipol_eq, jrad, ierr, m, n, u
    INTEGER :: status(MPI_STATUS_SIZE)
    
    m = size(ayxz, 1) !toroidal
    n = size(ayxz, 2) !radial
    allocate(ay(m))
    allocate(ayz(m, 0:mpol2))
    jrad = nrad/2 !choose a radial index
    
    ay(:) = ayxz(:, jrad, 1)
    call MPI_gather(ay,  m, MPI_real8, &
         &          ayz, m, MPI_real8, 0, tube_COMM, ierr)

    ! z boundary
    if(gclr==mpol2-1) call mpi_send(ayxz(:,jrad, 2), m, MPI_real8, 0, 123, tube_comm, ierr)
    if(gclr==0) call mpi_recv(ayz(:, mpol2), m, MPI_real8, mpol2-1, 123, tube_comm, status, ierr)

    if(myid==0) then
       full_file_name='ms/yz_txxxxxx'//partial_file_name
       write(full_file_name(8:13),'(i6.6)') kt
       open(newunit=u, file=full_file_name)
       do ipol = 0, mpol2 !z : poloidal
          ipol_eq=multi_eq_cells*ipol + 1 !index in the equilibrium grids
          do itor =1, mtor !y : toroidal 
             phi = ygrid(itor) + tor_shift_mc(ipol_eq, jrad)
             write(u,'(19ES16.5E3)') ygrid(itor), zgrid(ipol_eq), ayz(itor,ipol), &
                  & r_mc(ipol_eq,jrad), z_mc(ipol_eq,jrad), phi
          enddo
          write(u,*) 
          !if(ipol.eq.GCLR_cut) write(u,*) !to inform gnuplot that this is a new data block
       enddo
       close(u)
    endif
  end subroutine mode_structure_in_yz_plane


  subroutine nharmonics_in_poloidal_plane(kt, a, str)
    use constants, only: p_, one, zero, pi
    use domain_decomposition, only: dtheta2, myid, tube_comm, multi_eq_cells
    use control_parameters, only: nh_min, nh_max
    use magnetic_coordinates, only: mtor, nrad, nz=>mpol2, pfn, tfn, r=>r_mc, z=>z_mc
    use mpi
    implicit none
    integer, intent(in) :: kt
    character(*), intent(in) :: str
    complex(p_), intent(in) :: a(0:mtor-1, 1:nrad-2)
    complex(p_), allocatable :: atmp(:,:), atmp0(:,:,:)
    character(len=100) :: file_name
    real(p_) :: theta
    integer :: j, i, ieq, nh, ierr, u

    allocate(atmp(0:mtor-1, nrad), source=(zero,zero))    
    atmp(:, 2:nrad-1) = a(:,:)
    allocate(atmp0(0:mtor-1, nrad, 0:nz-1))
    
    call MPI_gather(atmp, mtor*nrad, MPI_complex16, &
         & atmp0, mtor*nrad, MPI_complex16, 0, tube_COMM, ierr)

    if(myid == 0) then
       do nh = nh_min, nh_max
          file_name='ms/poloidal_plane_txxxxxx_nhxxx'//str
          write(file_name(20:25),'(i6.6)') kt
          write(file_name(29:31),'(i3.3)') nh
          open(newunit=u, file=file_name)
          do j = 1, nrad
             do i = 0, nz-1
                ieq = 1 + i*multi_eq_cells
                theta = -pi + i*dtheta2
                write(u,'(19ES16.5E3)') r(ieq, j), z(ieq, j), 2*abs(atmp0(nh, j, i))/mtor, pfn(j), theta, tfn(j)
             enddo
          enddo
          close(u)
       enddo
    endif
  end subroutine nharmonics_in_poloidal_plane
  
  subroutine mode_structure_in_poloidal_plane(kt, a, str)
    use constants, only: p_, one
    use domain_decomposition, only: theta_start, dtheta2
    implicit none
    integer, intent(in) :: kt
    character(*), intent(in) :: str
    real(p_), intent(in) :: a(:, :, :)
    integer, parameter :: nz = 20 !along field line
    real(p_) :: theta(0:nz-1), c
    real(p_) :: perturb(size(a,1), size(a,2), 0:nz-1)
    real(p_) :: a0(size(a,1), size(a,2), 0:nz-1), a1(size(a,1),size(a,2), 0:nz-1)
    character(len=100) :: file_name
    integer :: m, n, j, i

    m = size(a,1) !toroidal
    n = size(a,2) !radial

    do i = 0, nz-1 !interpolate to more z gridpoints 
       theta(i) = theta_start + dtheta2/nz*i
       c = (theta(i) - theta_start)/dtheta2
       perturb(:,:,i) = a(:,:,1)*(one-c) + a(:,:,2)*c
    enddo

    do i = 0, nz-1
       do j = 1, n
          a0(:, j, i) = sum(perturb(:,j, i))/m
       enddo
    enddo

    file_name='ms/poloidal_plane_txxxxxx'//str//'_neq0'
    write(file_name(20:25),'(i6.6)') kt
    call to_poloidal_plane(a0, m, n, nz, theta, file_name)

    a1 = perturb - a0
    file_name='ms/poloidal_plane_txxxxxx'//str//'_nneq0'
    write(file_name(20:25),'(i6.6)') kt
    call to_poloidal_plane(a1, m, n, nz, theta, file_name)

  end subroutine mode_structure_in_poloidal_plane

  subroutine to_poloidal_plane(yxz, ny, nx, nz, theta, file_name)
    use constants, only: p_, pi, twopi
    use magnetic_coordinates, only: ygrid, r_mc, z_mc, tor_shift_mc, mpol, mpol2, nsegment, &
         & zgrid, pfn, tfn, toroidal_range
    use domain_decomposition, only: myid, tube_comm, dtheta2
    use interpolate_module
    use math, only: shift_toroidal
    use mpi
    implicit none
    integer, intent(in) :: ny, nx, nz
    real(p_), intent(in) :: yxz(ny, nx, 0:nz-1), theta(0:nz-1)
    character(len=*), intent(in) :: file_name
    real(p_), allocatable :: field(:,:), field0(:,:)
    real(p_) :: tmp(ny+1)
    real(p_) :: phi, my_theta, alpha, r, z, th, tor_shift
    integer :: j, i, iz,  u, ierr

    allocate(field(nx, 0:nz-1))
    allocate(field0(nx, 0:nz*mpol2-1))

    !choose a cylindrical toroidal angle, for which the mode structure is computed
    phi = 0.5_p_*twopi/nsegment

    do j = 1, nx ! radial
       do i = 0, nz-1 ! local theta range
          call linear_1d_interpolate(mpol, zgrid, tor_shift_mc(:,j), theta(i), tor_shift)
          alpha = phi - tor_shift ! alpha of the fixed phi
          call shift_toroidal(alpha, toroidal_range)
          !Interpolate in alpha, to the same cylindrical toroidal angle:
          tmp(1:ny) = yxz(:, j, i)
          tmp(ny+1) = yxz(1, j, i) !periodic boundary condition
          call linear_1d_interpolate(ny+1, ygrid, tmp, alpha, field(j, i)) 
       enddo
    enddo

    call MPI_gather(field, nx*nz, MPI_real8, field0, nx*nz, MPI_real8, 0, tube_COMM, ierr)

    if(myid == 0) then
       open(newunit=u, file=file_name)
       do j = 1, nx
          do i = 0, mpol2 - 1
             my_theta = -pi + i*dtheta2
             do iz = 0, nz-1
                th = my_theta + dtheta2/nz*iz
                call linear_1d_interpolate(mpol, zgrid, r_mc(:,j), th, r)
                call linear_1d_interpolate(mpol, zgrid, z_mc(:,j), th, z)
                write(u,'(19ES16.5E3)') r, z, field0(j, i*nz+iz), pfn(j), th, tfn(j)
             enddo
          enddo
       enddo
       close(u)
    endif
  end subroutine to_poloidal_plane

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
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR
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

    call MPI_gather(my_a_xz_plane, n, MPI_complex16, &
         & a_xz_plane,    n, MPI_complex16, 0, tube_COMM, ierr)
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


subroutine field_lines_analyse()
  use constants,only:p_
  use constants,only: two
  use radial_module,only: z_axis
  implicit none
  
  integer:: n_tor_loop,max_npt_along_field_line,krad !used in field line tracing module
!  namelist /field_line_tracing_nl/  n_tor_loop,max_npt_along_field_line,krad
  real(p_),allocatable:: r_start(:),z_start(:),phi_start(:) !starting point of field lines
  real(p_),allocatable:: r_poincare(:,:),z_poincare(:,:),phi_poincare(:,:) !pointcare points
  integer:: j,k
  integer,allocatable::nloop_actual(:)

n_tor_loop=10
max_npt_along_field_line=8000
krad=1

  
!!$  open(31,file='input.nmlt')
!!$  read(31,field_line_tracing_nl)
!!$  close(31)
!!$  write(*,field_line_tracing_nl)

  allocate(r_start(krad))
  allocate(z_start(krad))
  allocate(phi_start(krad))
  allocate(r_poincare(n_tor_loop+1,krad))
  allocate(z_poincare(n_tor_loop+1,krad))
  allocate(phi_poincare(n_tor_loop+1,krad))
  allocate(nloop_actual(krad))

!  !$omp parallel do
  do k=1,krad
     r_start(k)=2.15_p_+(k-1)*0.15d0/(20-1)
     z_start(k)=z_axis
     phi_start(k)=0._p_
     call field_line_tracing(r_start(k),z_start(k),phi_start(k), max_npt_along_field_line,n_tor_loop,&
          & r_Poincare(:,k),z_Poincare(:,k),phi_Poincare(:,k),nloop_actual(k))
  enddo
!  !$omp end parallel do

  open(26,file='poincare.txt')
  do k=1,krad
     do j=1,nloop_actual(k)
        write(26,*) r_Poincare(j,k),z_Poincare(j,k),phi_Poincare(j,k)
     enddo
     write(26,*)
     write(26,*)
  enddo
  close(26)

!!$  do k=1,krad
!!$     call draw_magnetic_surface(r_start(k),z_start(k),'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
!!$  enddo


end subroutine field_lines_analyse



subroutine field_line_tracing(r0,z0,phi0,npt,n_tor_loop,r_poincare,z_poincare,phi_poincare,nloop_actual)
  !given coordinates (R,Z,phi), this subroutine follows the field lines passing through this point until it has finish n_tor_loop toroidal loop or exceeds the specifed maximum number of points along the field-line, npt. This subroutine also calculates the safety factor of the field-line found.
  use constants,only:p_,two,twopi,one_half
  use boundary,only: np_lcfs,x_lcfs,z_lcfs,nlim,rlim,zlim !use to check whether field line touch the boundary
 use magnetic_field, only : br,bz,bphi
 implicit none

  real(p_),intent(in):: r0,z0,phi0
  integer,intent(in)::npt,n_tor_loop
  real(p_):: r(npt),z(npt),phi(npt)
  real(p_),intent(out):: r_poincare(n_tor_loop+1),z_poincare(n_tor_loop+1),phi_poincare(n_tor_loop+1)
  integer,intent(out):: nloop_actual
  real(p_),parameter:: step=1d-3 !meter, trial of dr or dz step
  real(p_):: brval,bzval,bphival,bpolval,dr,dz,dphi
  
  real(p_):: r_mid,z_mid,dl_pol,qval
  logical:: loss
  integer:: j,k,jj

  k=1 !Poincare points
  r_poincare(k)=r0 !Poincare points
  z_poincare(k)=z0
  phi_poincare(k)=phi0

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  loss=.false.
  do j=1,npt-1
     !2nd Runge-Kutta
     brval=    br(r(j),z(j))
     bzval=    bz(r(j),z(j))
     bphival=bphi(r(j),z(j))
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr=step*one_half
        if(brval.lt.0._p_) dr=-step*one_half
        dz=bzval/brval*dr
     else
        dz=step*one_half
        if(bzval.lt.0._p_) dz=-step*one_half
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/bpolval*dl_pol/r(j)

     !first step:
     r_mid=r(j)+dr
     z_mid=z(j)+dz
     brval=    br(r_mid,z_mid)
     bzval=    bz(r_mid,z_mid)
     bphival=bphi(r_mid,z_mid)
     bpolval=sqrt(brval**2+bzval**2)


     if(abs(bzval).lt.abs(brval)) then
        dr=step
        if(brval.lt.0._p_) dr=-step
        dz=bzval/brval*dr
     else
        dz=step
        if(bzval.lt.0._p_) dz=-step
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/(bpolval*r_mid)*dl_pol

     r(j+1)=r(j)+dr
     z(j+1)=z(j)+dz
     phi(j+1)=phi(j)+dphi

     call  check_whether_field_line_touch_boundary(r(j+1),z(j+1),phi(j+1),x_lcfs,z_lcfs,np_lcfs,loss)
     if (loss.eqv. .true.) exit

     if(abs(floor(abs(phi(j)-phi0)/twopi)-floor(abs(phi(j+1)-phi0)/twopi)).eq.1) then ! finish one toroidal turn
        !   write(*,*) 'j=',j,'k=',k, phi(j),phi(j+1)
        k=k+1
        r_poincare(k)=(r(j)+r(j+1))/two
        z_poincare(k)=(z(j)+z(j+1))/two
        phi_poincare(k)=((phi(j)+phi(j+1))/two)/twopi
     endif

     if(abs(phi0-phi(j+1))/twopi.ge.n_tor_loop) exit

  enddo

 if(j.eq.npt) then
     open(76,file='bad_line.txt')
     do jj=1,j-1
     write(76,*) r(jj),z(jj),phi(jj)
     enddo
     close(76)
     call safety_factor_a_field_line(r,z,phi,j,qval)
     call draw_magnetic_surface(r0,z0,'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
     stop 'max number of tracing steps of field line is exceeded before achiving the specified number of toroidal loop'
endif

  nloop_actual=k
  write(*,*) 'nloop_actual=',nloop_actual, 'actual step along field line=',j

  call safety_factor_a_field_line(r,z,phi,j+1,qval)

!!$  open(76,file='field_line.txt')
!!$  do jj=1,j+1
!!$     write(76,*) r(jj),z(jj),phi(jj)
!!$  enddo
!!$  close(76)


!call check_field_line_in_field_aligned_coordinates(r,z,phi,j,qval)

end subroutine field_line_tracing


subroutine field_line_tracing_simplified(r0,z0,phi0) !see the comments in subroutine "field_line_tracing"
  use constants,only:p_,two,twopi,one_half
use magnetic_field, only :  br,bz,bphi
  implicit none
  real(p_),intent(in):: r0,z0,phi0
  integer,parameter::npt=3000
  real(p_):: r(npt),z(npt),phi(npt)
  real(p_),parameter:: step=1d-3 !meter, trial of dr or dz step
  real(p_):: brval,bzval,bphival,bpolval,dr,dz,dphi
  real(p_):: r_mid,z_mid,dl_pol,qval
  integer:: j,u

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  do j=1,npt-1 !2nd Runge-Kutta
     brval=    br(r(j),z(j))
     bzval=    bz(r(j),z(j))
     bphival=bphi(r(j),z(j))
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr=step*one_half
        if(brval.lt.0._p_) dr=-step*one_half
        dz=bzval/brval*dr
     else
        dz=step*one_half
        if(bzval.lt.0._p_) dz=-step*one_half
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/bpolval*dl_pol/r(j)

     !first step of 2nd Runge-Kutta:
     r_mid=r(j)+dr
     z_mid=z(j)+dz
     brval=    br(r_mid,z_mid)
     bzval=    bz(r_mid,z_mid)
     bphival=bphi(r_mid,z_mid)
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr=step
        if(brval.lt.0._p_) dr=-step
        dz=bzval/brval*dr
     else
        dz=step
        if(bzval.lt.0._p_) dz=-step
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/(bpolval*r_mid)*dl_pol

     r(j+1)=r(j)+dr
     z(j+1)=z(j)+dz
     phi(j+1)=phi(j)+dphi

  enddo

  open(newunit=u,file='field_line')
  do j=1,npt
     write(u,*) phi(j),z(j),r(j)
  enddo
  close(u)
end subroutine field_line_tracing_simplified


subroutine field_line_tracing0(r0,z0,phi0, npt, dphi, r,z,phi)
  use constants,only : p_, two, twopi, one_half
  use magnetic_field, only :  br, bz, bphi
  implicit none
  real(p_),intent(in) :: r0,z0,phi0
  real(p_),intent(in) :: dphi
  integer, intent(in) :: npt
  real(p_), intent(out) :: r(npt),z(npt),phi(npt)
  real(p_) :: brval,bzval,bphival
  real(p_) :: dr, dz, r_mid, z_mid
  integer :: j

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  do j=1,npt-1 !2nd Runge-Kutta
     brval=    br(r(j),z(j))
     bzval=    bz(r(j),z(j))
     bphival=bphi(r(j),z(j))

     dr= brval/bphival*(r(j)*dphi)
     dz= bzval/bphival*(r(j)*dphi)

     r_mid=r(j)+dr*one_half !first step of 2nd Runge-Kutta
     z_mid=z(j)+dz*one_half
     brval=    br(r_mid,z_mid)
     bzval=    bz(r_mid,z_mid)
     bphival=bphi(r_mid,z_mid)

     dr= brval/bphival*(r_mid*dphi)
     dz= bzval/bphival*(r_mid*dphi)
    
     r(j+1)=r(j)+dr !second step of 2nd Runge-Kutta
     z(j+1)=z(j)+dz
     phi(j+1)=phi(j)+dphi
  enddo
end subroutine field_line_tracing0

subroutine check_field_line_in_field_aligned_coordinates(r,z,phi,npt,qval)
  use constants, only: p_, two,twopi
  use radial_module,only:z_axis
  use mapping_module,only: nx_mapping ,j0,r_cyl,tor_shift_b
  use map_to_mc, only : interpolate_from_cylindrical_to_magnetic_coordinates
  use magnetic_field, only : psi_func, pfn_func, radcor_as_func_of_pfn
  implicit none
  integer,intent(in):: npt
  real(p_),intent(in):: r(npt),z(npt),phi(npt),qval
  real(p_):: radcor(npt),theta(npt),alpha(npt),tor_shift(npt)
  real(p_)::x1,x2,y1,y2,z1,z2,dx,dy,dz,dl(npt)
  real(p_):: sum=0._p_,real_shift,total_shift,tmp_array(nx_mapping)
  integer:: j,kk



  do j=1,npt
     radcor(j)=radcor_as_func_of_pfn(pfn_func(r(j),z(j))) !get radial coordinate
     call interpolate_from_cylindrical_to_magnetic_coordinates(r(j),z(j),theta(j),tor_shift(j))
  enddo

  dl(1)=0._p_
  do j=2,npt
     x1=r(j-1)*cos(phi(j-1))
     x2=r(j)*cos(phi(j))
     y1=r(j-1)*sin(phi(j-1))
     y2=r(j)*sin(phi(j))
     z1=z(j-1)
     z2=z(j)
     dx=x2-x1
     dy=y2-y1
     dz=z2-z1
     dl(j)=dl(j-1)+sqrt(dx*dx+dy*dy+dz*dz)
  enddo

  open(11,file='field_line_in_field_aligned_co.txt')
  do j=2,npt

!!$     if(abs(theta(j)-theta(j-1)) .ge. twopi*0.9) then !indecate finishing one poloidal loop
!!$        !total_shift=tor_shift(j)*(z_axis-z(j-1))/(z(j)-z_axis)+tor_shift(j-1)
!!$        !total_shift=qval*twopi*1.0004
!!$        !total_shift=2.35227*twopi
!!$        do kk=1,nx_mapping
!!$           tmp_array(kk)= tor_shift_b(kk,j0)
!!$        enddo
!!$        call linear_1d_interpolate(nx_mapping,r_cyl,tmp_array,r(j),total_shift)
!!$        sum=sum+total_shift
!!$     endif
!!$     alpha(j)=phi(j)-(tor_shift(j)+sum)
   call accumulate_tor_shift(theta(j-1),theta(j),r(j),tor_shift(j),real_shift)
    alpha(j)=phi(j)-real_shift 
     write(11,*) dl(j),radcor(j),theta(j),alpha(j), tor_shift(j), phi(j),r(j),z(j)
  enddo
  close(11)
end subroutine check_field_line_in_field_aligned_coordinates


subroutine safety_factor_a_field_line(r,z,phi,npt,qval)
  !given a field line, calculate its safety factor
  use constants, only: p_, two, twopi
  use magnetic_field, only : psi_func
  use magnetic_field, only : qfunc0
  implicit none
  integer, intent(in) :: npt
  real(p_), intent(in) :: r(npt), z(npt), phi(npt)
  real(p_), intent(out) :: qval
  real(p_) :: phi_old
  integer :: j, npass

  npass = 0
  do j=1,npt-1
     if(z(j)*z(j+1).lt.0) then !indicates one midplane-crossing
        npass = npass+1
        if(npass == 1) phi_old = (phi(j)+phi(j+1))/two
     endif
     if(npass == 3) then !indicates that the line has finished one poloidal period
        qval=abs(phi_old-phi(j))/twopi
        write(*,*) 'safety factor of field line passing (r,z)',r(1),z(1),'is', qval,&
             & 'q value specified in gfile =', qfunc0(psi_func(r(1),z(1)))
        exit
     endif
  enddo

  write(*,*) 'toroidal loops the field line travels=',(phi(npt)-phi(1))/twopi

end subroutine safety_factor_a_field_line


subroutine check_whether_field_line_touch_boundary(r,z,phi,rlim,zlim,nlim,loss)
  use constants,only:p_
  use math, only : pnpoly
!  use boundary,only: nlim,rlim,zlim 
  implicit none
  real(p_),intent(in):: r,z,phi
  integer,intent(in):: nlim
  real(p_),intent(out):: rlim(nlim),zlim(nlim)
  logical,intent(out):: loss
  integer:: inout

  call PNPOLY(r,z,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter
  !        if (inout.eq.1) then !within the LCFS
  if (inout.eq.-1 .or.inout.eq.0) then !the particle is out of the limiter
     write(*,*) '==>This field line touches the limiter at (R,Z,phi)=', r,z,phi
     loss=.true.
     !stop
  else
     loss=.false.
  endif

end subroutine


  subroutine accumulate_tor_shift(theta_old,theta_new,r,tor_shift,real_shift) !,kt)
    use constants,only:p_
    use constants,only: twopi
    use mapping_module,only: nx_mapping,tor_shift_b,j0,r_cyl
  use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_),intent(in):: theta_old,theta_new,r,tor_shift
    real(p_),intent(out)::real_shift
    real(p_),save::sum=0._p_
    real(p_):: tmp_array(nx_mapping),twopi_q
!integer,intent(in):: kt
    integer::kk
    if(abs(theta_old-theta_new) .ge. twopi*0.9) then !indicate finishing one poloidal loop
       do kk=1,nx_mapping
          tmp_array(kk)= tor_shift_b(kk,j0)
       enddo
       call linear_1d_interpolate(nx_mapping,r_cyl,tmp_array,r,twopi_q) !the result is twopi*q, I use this instead of directly using twopi*q because the latter may cause some cancellation problem
       sum=sum+twopi_q*sign(1._p_,theta_old-theta_new)
!write(*,*) 'sum=',sum, 'twopi_q=',twopi_q,'tor_shift=',tor_shift,'sum+tor_shift/ real_shift',sum+tor_shift, 'kt=',kt
    endif

    real_shift=sum+tor_shift


    !real_shift=sum-(twopi_q-tor_shift)
  end subroutine accumulate_tor_shift


  subroutine draw_magnetic_surface(r0,z0,filename) !draw the magnetic surface which passes through (r0,z0)
  use constants,only:p_
  use constants,only:zero,one,two,twopi
  use boundary,only:np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: r_axis,z_axis
  use contour_mod,only : find_contour
  use magnetic_field, only : psi_func
  implicit none

  real(p_),intent(in):: r0,z0
  character(*),intent(in)::  filename
  real(p_):: psival
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  integer:: i,u
  
  psival=psi_func(r0,z0)

  call find_contour(psival,x_contour,z_contour)

  open(newunit=u,file=filename)
  do i=1,np_lcfs
     write(u,*) x_contour(i),z_contour(i)
  enddo
  close(u)
end subroutine draw_magnetic_surface
