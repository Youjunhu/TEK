module diagnostic_mod
contains

subroutine parallel_current_diagnostic(niter,iter)
use normalizing,only: vn_i
use fk_module,only: ni0,charge_i
use gk_module,only: ne0
!use perturbation_field_matrix,only: jpar_i_old, jpar_e_old

integer,intent(in):: niter,iter

!if(iter.eq.niter) write(*,*) iter, jpar_e_old(5,5)*vn_e*ne0*(-charge_i),jpar_i_old(5,5)*vn_i*ni0*charge_i,&
!     & jpar_e_old(5,5)*vn_e*ne0*(-charge_i)/(jpar_i_old(5,5)*vn_i*ni0*charge_i)

 !write(*,*) jpar_e_old(15,15)*vn_e*ne0*(-charge_i),jpar_i_old(15,15)*vn_i*ni0*charge_i,&
  !   & jpar_e_old(15,15)*vn_e*ne0*(-charge_i)/(jpar_i_old(15,15)*vn_i*ni0*charge_i)

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

subroutine compute_heat_flux(time,ns,nm, touch_bdry_e, mu_e, vpar_e, w_e,&
     &  radcor_e, bval_e,grad_psi_e, perturbed_radial_drift)
  use constants,only:p_
  use mpi
  use constants,only: two,kev
  use control_parameters,only: kstart, kend
  !  use gk_module,only:ps_vol_e
  use magnetic_coordinates,only: radial_width, dradcor, dtheta, j_low2, mpol, toroidal_range, radcor_1d_array2
  use magnetic_coordinates,only: jacobian, vol
  use gk_module,only: te0,ne0, kappa_te,mass_e, rho_e, vn_e
  use domain_decomposition,only: myid
  implicit none
  real(p_),intent(in):: time
  integer,intent(in):: ns, nm
  logical,intent(in):: touch_bdry_e(:)
  real(p_),intent(in):: mu_e(:), vpar_e(:), w_e(:),&
       &  radcor_e(:), bval_e(:),grad_psi_e(:), perturbed_radial_drift(:)
  integer,parameter:: nsub=10
  real(p_):: tmp_radial_array(nsub), dv(nsub), myheat_flux(nsub), heat_flux(nsub),heat_flux_vol_av
  real(p_):: dx, kinetic_energy, factor, gyro_bohm_diffusivity
  integer:: i,j,k, jeq,ierr
  character(len=64)::filename3
  logical,save:: is_first=.true.
  integer,save:: unit_flux
  if ((is_first.eqv..true.).and.(myid.eq.0)) then 
!!$     filename3 = 'heat_fluxxxxxxx_xxxxxx'
!!$     write(filename3(10:22),'(i6.6,a1,i6.6)') kstart,'_',kend
     filename3 = 'heat_flux.txt'
     open(newunit=unit_flux, file=filename3)
     is_first=.false.
  endif

  factor=mass_e(ns)*vn_e(ns)**3/(ne0(ns)*te0(ns)*kev*kappa_te(ns)) !to transform heat flux to heat diffusivity in SI units
  gyro_bohm_diffusivity=rho_e(ns)**2*sqrt(te0(ns)*kev/mass_e(ns))*kappa_te(ns) !in SI units

  dx=radial_width/nsub
  do j=1,nsub
     tmp_radial_array(j)=radcor_1d_array2(1)+dx*(j-1)
  enddo

  myheat_flux=0._p_
  do k=1,nm
     if(touch_bdry_e(k).eqv..false.) cycle
     j=floor((radcor_e(k)-tmp_radial_array(1))/dx)+1
     kinetic_energy=vpar_e(k)**2/two+mu_e(k)*bval_e(k)
     !myheat_flux(j)=myheat_flux(j)+w_e(k)*ps_vol_e(k)*kinetic_energy*perturbed_radial_drift(k) !a bug, ps_vol_e is not needed
     myheat_flux(j)=myheat_flux(j)+w_e(k)*kinetic_energy*perturbed_radial_drift(k)/grad_psi_e(k) !corrected
  enddo

  call MPI_Reduce(myheat_flux, heat_flux, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr)

  if(myid.eq.0) then
     heat_flux_vol_av=sum(heat_flux)/vol*factor/gyro_bohm_diffusivity
     dv(:)=0._p_
     do j=1,nsub
        jeq=j_low2+(j-1)*int(dx/dradcor)
        do i=1,mpol-1
           dv(j)=dv(j)+abs(jacobian(i,jeq))*dx*dtheta*toroidal_range
        enddo
     enddo
     heat_flux(:)=heat_flux(:)/dv(:)*factor/gyro_bohm_diffusivity
     write(unit_flux,'(20(1pe20.8))') time, heat_flux_vol_av,sum(heat_flux)/nsub,(heat_flux(j),j=1,nsub)
  endif

end subroutine compute_heat_flux




subroutine draw_grids_on_theta_isosurface(mpol,nflux,tor_shift_mc,r_mc,z_mc) !on theta=constant surface, the subroutine name is wrong, not necessarily top view
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nflux
  real(p_),intent(in)::tor_shift_mc(mpol,nflux),r_mc(mpol,nflux),z_mc(mpol,nflux)
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
     !     do j=1,nflux,10
     !    do j=1,nflux,1
     do j=1,1
        phi=alpha+tor_shift_mc(i,j) !phi is changing due to the radial dependene of tor_shift
        write(113,*) phi,r_mc(i,j),z_mc(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

  open(113,file='grids_on_theta_isosurface2.txt')
  do j=1,nflux,10
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


subroutine draw_alpha_isosurface(mpol,nflux,tor_shift_mc,r_mc,z_mc) !on a nonzero radial range
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nflux
  real(p_),intent(in)::tor_shift_mc(mpol,nflux),r_mc(mpol,nflux),z_mc(mpol,nflux)
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


subroutine draw_alpha_contours_on_a_magnetic_surface(mpol,nflux,tor_shift_mc,r_mc,z_mc) !field lines on a magnetic surface
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nflux
  real(p_),intent(in)::tor_shift_mc(mpol,nflux),r_mc(mpol,nflux),z_mc(mpol,nflux)
  real(p_):: phi,alpha0,tor_range
  integer:: i,j,ialpha,ishift,u,iphi
  integer,parameter::nalpha=10,nphi=20

  open(newunit=u,file='alpha_contours_on_magnetic_surface.txt')
  !do j=1,30,2
  j=nflux/2 !select a radial location, i.e., a magnetic surface 
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


end module diagnostic_mod


module mode_structure
  implicit none
contains
  subroutine mode_structure_on_xy_plane(kt,GCLR,a,partial_file_name)
    use constants,only:p_
    use magnetic_coordinates,only:radcor_1d_array2,tor_1d_array
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
          write(u,*) radcor_1d_array2(j),tor_1d_array(i),a(i,j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine mode_structure_on_xy_plane

  subroutine mode_structure_on_xz_plane(kt,a,partial_file_name)
    use constants,only:p_
    use constants,only:pi
    use magnetic_coordinates,only: radcor_1d_array2,theta_1d_array,dtheta,mtor
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
             write(u,*) radcor_1d_array2(j), my_theta,a_xz_plane(j,ipol)
          enddo
          write(u,*)
          !if(ipol.eq.GCLR_cut) write(u,*) !to inform gnuplot that this is a new data block, to prevent gnuplot using line connection between the following data and previous data
       enddo
       close(u)
    endif
  end subroutine mode_structure_on_xz_plane


  subroutine mode_structure_on_yz_plane(kt,a,partial_file_name)
    use constants,only: p_, pi, twopi
    use magnetic_coordinates,only: mtor,tor_1d_array,theta_1d_array, dtheta,nflux2, mpol2,&
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
    jrad=nflux2/2 !choose a radial index
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
             phi=tor_1d_array(itor)+tor_shift_mc(ipol_eq,jrad)
             write(u,'(19ES16.5E3)') tor_1d_array(itor), theta_1d_array(ipol_eq), a_yz_plane(itor,ipol), &
                  & r_mc(ipol_eq,jrad), z_mc(ipol_eq,jrad), phi
          enddo
          write(u,*) 
          !if(ipol.eq.GCLR_cut) write(u,*) !to inform gnuplot that this is a new data block, to prevent gnuplot using line connection between the following data and previous data
       enddo
       close(u)
    endif
  end subroutine mode_structure_on_yz_plane


  subroutine mode_structure_on_poloidal_plane(kt,a)
    use constants,only : p_
    use transform_module,only:oned_fourier_transform1, oned_backward_fourier_transform1
    implicit none
    integer,intent(in)::kt
    real(p_),intent(in):: a(:,:)
    complex(p_):: a_dft(size(a,1),size(a,2))
    real(p_) :: a0(size(a,1),size(a,2)), a1(size(a,1),size(a,2))
    integer:: m,n, j
    character(len=100) :: file_name
    
    m=size(a,1) !toroidal
    n=size(a,2) !radial
    call oned_fourier_transform1(a,a_dft,m,n) !calculating 1d DFT of s(:,:) along the first dimension

    do j=1,n
       a0(:,j)=real(a_dft(1,j)) !space distribution of the n=0 harmonic
    enddo
    file_name='ms/poloidal_plane_txxxxxxneq0'
    write(file_name(20:25),'(i6.6)') kt
    call mode_structure_on_poloidal_plane0(a0, file_name)

    a_dft(1,:)=0 !remove the n=0 component for ploting n.neq.0 mode structure, this removing does not affect the simulation
    call oned_backward_fourier_transform1(a_dft,a1,m,n)
    file_name='ms/poloidal_plane_txxxxxxnneq0'
    write(file_name(20:25),'(i6.6)') kt
    call mode_structure_on_poloidal_plane0(a1,file_name)

  end subroutine mode_structure_on_poloidal_plane

  subroutine mode_structure_on_poloidal_plane0(a, file_name)
    use constants,only:p_, pi, twopi
    use magnetic_coordinates,only: tor_1d_array,r_mc, z_mc,tor_shift_mc,mpol,nsegment,&
         & dtheta,theta_1d_array, pfn, j_low2
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,theta_start,dtheta2,GCLR_cut
    use interpolate_module
    use math,only: shift_to_specified_toroidal_range
    use mpi
    implicit none
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in):: file_name
    integer:: itor,ipol,ipol_eq,j,ierr,m,n
    real(p_):: my_a_poloidal_plane(size(a,2)),a_poloidal_plane(size(a,2),0:numprocs/ntube-1)
    real(p_)::phi_val,my_theta, alpha
    integer:: jeq, u !file unit number

    m=size(a,1) !toroidal
    n=size(a,2) !radial

    phi_val=0._p_+0.5_p_*twopi/nsegment !choose a fixed cylindrical toroidal angle, this is the poloidal plane on which the mode structure is computed
    ipol_eq=1+nint((theta_start-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
    do j=1,n
       jeq=j+j_low2-1
       alpha = phi_val - tor_shift_mc(ipol_eq,jeq)
       call shift_to_specified_toroidal_range(alpha)
       call linear_1d_interpolation(m, tor_1d_array, a(:,j), alpha ,my_a_poloidal_plane(j)) !interpolated to the same cylindrical toroidal angle
    enddo
    call MPI_gather(my_a_poloidal_plane, n, MPI_real8,&
         a_poloidal_plane, n, MPI_real8, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       open(newunit=u,file=file_name)
       do j=1,n
          jeq=j-1+j_low2
          do ipol=0,numprocs/ntube-1
             my_theta=-pi+ipol*dtheta2
             ipol_eq=1+nint((my_theta-theta_1d_array(1))/dtheta) !index of begining theta angle of this subdomain in the original magnetic coordinate grids
             write(u,'(19ES16.5E3)') r_mc(ipol_eq,jeq), z_mc(ipol_eq,jeq), a_poloidal_plane(j,ipol), pfn(jeq), my_theta
          enddo
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
    use magnetic_coordinates,only: radcor_1d_array2,theta_1d_array,dtheta,mtor
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
  use magnetic_coordinates, only : r_mc, z_mc, mpol,nflux, mtor, tor_1d_array, tor_shift_mc
  use magnetic_field, only : psi_func, qfunc0
  use math, only : shift_to_specified_toroidal_range
  implicit none
  integer, parameter :: npt=2000
  integer :: u,i,j,k
  real(p_) :: phi, phi0
  real(p_) :: rf(npt), zf(npt), phif(npt), dphi, qval

  open(newunit=u,file='grid.txt')
  do k=1, mtor
     do i=1, mpol
        do j=1,nflux
           phi = tor_1d_array(k) + tor_shift_mc(i,j)
           !call shift_to_specified_toroidal_range(phi) 
           write(u,*) r_mc(i,j), z_mc(i,j), phi
        enddo
     enddo
  enddo
  close(u)

  i = (mpol+1)/2
  j = nflux/2
  !phi0 = tor_1d_array(1) + tor_shift_mc(i,j)
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
  phi0=tor_1d_array(1) + tor_shift_mc(i,j)
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
  j = nflux
  phi0=tor_1d_array(1) + tor_shift_mc(i,j)
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
     write(u,*)  r_mc(i, nflux), z_mc(i, nflux)
  enddo
  close(u)

end subroutine visualize_grid


