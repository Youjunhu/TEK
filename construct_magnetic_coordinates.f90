subroutine construct_magnetic_coordinates()
  use constants,only : p_, zero, twopi,pi
  use control_parameters,only: diagnosis
  use boundary,only: x_lcfs, z_lcfs, np_lcfs
  use radial_module,only: r_axis,z_axis,psi_lcfs,psi_axis
  use magnetic_field,only: radcor_minor_r, minor_r_radcor,minor_r_prime  !function name
  use magnetic_coordinates,only: mpol, mtor, nsegment, nflux, pfn, GSpsi_array, &
       & r_mc,z_mc, & !as output
       & theta_1d_array, dtheta, & !as output
       & minor_r_array,minor_r_prime_array,i_theta_zero, dl_mc, & !output
       & toroidal_range, dtor, tor_1d_array !output
  use domain_decomposition,only: myid
  use contour_mod,only : contour
  use splines, only : spline3ders
  use math,only:  one_dimensional_derivative, arc_length
  implicit none
  real(p_), allocatable :: r_mag_surf0(:,:), z_mag_surf0(:,:)
  integer :: i, j

  call choose_radial_grids()
  allocate(r_mag_surf0(np_lcfs,nflux), z_mag_surf0(np_lcfs,nflux))
  do j=1,nflux
     call contour(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis,GSpsi_array(j),r_mag_surf0(:,j),z_mag_surf0(:,j))
  enddo

  if((myid.eq.0) .and. (diagnosis.eqv..true.)) call diagnostic1()

  allocate(minor_r_array(nflux))
  allocate(minor_r_prime_array(nflux))
  do j=1,nflux
     minor_r_array(j) = r_axis - r_mag_surf0(1,j) !defined on the high-field-side midplane
  enddo
  call one_dimensional_derivative(nflux,pfn,minor_r_array,minor_r_prime_array) !to be used as interpolating table
  !call spline3ders(pfn(:), minor_r_array(:), pfn(:), dynew=minor_r_prime_array(:))

  dtheta = twopi/(mpol-1) !poloidal grid spacing for equilibrium
  allocate(theta_1d_array(mpol))
  theta_1d_array = [ (-pi+dtheta*(i-1), i=1,mpol) ]
  i_theta_zero=(mpol+1)/2 !poloidal index corresponding to theta=0 (mpol is assumed odd)

  allocate(r_mc(mpol,nflux))
  allocate(z_mc(mpol,nflux))

  do j=1,nflux
     call construct_poloidal_coordinate(r_mag_surf0(:,j),z_mag_surf0(:,j),np_lcfs, &
          & mpol,theta_1d_array, r_mc(:,j), z_mc(:,j))
  enddo

  allocate(dl_mc(mpol-1,nflux))
  do j=1,nflux
     call arc_length(r_mc(:,j), z_mc(:,j), mpol, dl_mc(:,j))
  enddo

  call calculate_metric() 

  toroidal_range = twopi/nsegment
  dtor = toroidal_range/mtor
  allocate(tor_1d_array(mtor+1)) 
  tor_1d_array = [ (zero+dtor*(i-1), i=1,mtor+1) ]

  call plasma_volume_of_computational_region()
  if ((myid==0) .and. (diagnosis .eqv. .true.)) call diagnostic2()
  if((myid==0) .and. (diagnosis.eqv..true.)) call diagnostic3()
  deallocate(r_mag_surf0, z_mag_surf0)

contains

  subroutine diagnostic1()
    integer:: u
    open(newunit=u,file='mag_surf_shape0.txt')
    do j=1,nflux
       do i=1,np_lcfs
          write(u,*) r_mag_surf0(i,j),z_mag_surf0(i,j)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)
  end subroutine diagnostic1

  subroutine diagnostic2()
    integer:: u
    open(newunit=u,file='theta_line.txt')
    do i=1,mpol
       do j=1,nflux
          write(u,*) r_mc(i,j),z_mc(i,j),theta_1d_array(i)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)

    open(newunit=u,file='mag_surf_shape.txt')
    do j=1,nflux
       do i=1,mpol
          write(u,*) r_mc(i,j),z_mc(i,j)
       enddo
       write(u,*);  write(u,*)
       write(u,*);  write(u,*)
    enddo
    close(u)
  end subroutine diagnostic2

  subroutine diagnostic3()
    integer:: u
    open(newunit=u,file='minor_r.txt')
    do j=1,nflux
       write(u,*) pfn(j), minor_r_array(j),  minor_r_prime_array(j)
    enddo
    close(u)
  end subroutine diagnostic3

end subroutine construct_magnetic_coordinates



subroutine choose_radial_grids()
  use constants,only:p_
  use radial_module,only: psi_axis,psi_lcfs
  use magnetic_field,only: radcor_as_func_of_pfn
  use magnetic_coordinates,only: nflux2, pfn_inner, pfn_bdry, &
       & nflux, pfn,GSpsi_array,GSpsi_prime, radcor_1d_array,dradcor, & !output
       & radcor_low0,radcor_upp0, &  !as output
       & radcor_low2,radcor_upp2, j_low2,j_upp2, radcor_1d_array2, & !as output
       & radial_width !as output
  use flux_tube_model,only: j_fixed, radcor_fixed
  implicit none
  integer,parameter:: points_in_buffer=0 !additional grid points in the buffer region
  real(p_):: dpfn
  integer:: j

  nflux=nflux2+2*points_in_buffer
  allocate(pfn(nflux))  
  allocate(GSpsi_array(nflux))
  allocate(radcor_1d_array(nflux))
  dpfn=(pfn_bdry-pfn_inner)/real(nflux2-1)
  if((pfn_inner-points_in_buffer*dpfn<0.01) .or. (pfn_bdry+points_in_buffer*dpfn>0.99)) then
     stop "radial buffer region is too near the magnetic axis or the LCFS, decrease the points_in_buffer to correct this"
  endif
  do j=1,nflux !select some flux surfaces (labeld by GSpsi_array)
     pfn(j)=pfn_inner-points_in_buffer*dpfn+dpfn*(j-1)
     GSpsi_array(j)=psi_axis+pfn(j)*(psi_lcfs-psi_axis) !GSpsi=Aphi*R, this array is used in finding magnetic surfaces
     radcor_1d_array(j)=radcor_as_func_of_pfn(pfn(j))
  enddo

  radcor_low0=radcor_1d_array(1)
  radcor_upp0=radcor_1d_array(nflux)
  j_low2=points_in_buffer+1
  j_upp2=points_in_buffer+nflux2
  radcor_low2=radcor_1d_array(j_low2)
  radcor_upp2=radcor_1d_array(j_upp2)
  allocate(radcor_1d_array2(nflux2))
  do j=1,nflux2
     radcor_1d_array2(j)=radcor_1d_array(j_low2+j-1)
  enddo
  radial_width=radcor_1d_array2(nflux2)-radcor_1d_array2(1)
  dradcor = radcor_1d_array2(2) - radcor_1d_array2(1) !radial grid interval
  GSpsi_prime = psi_lcfs - psi_axis !dGSpsi/dx, x is the normalized poloidal magnetic flux

  j_fixed = (j_low2+j_upp2)/2
  radcor_fixed=radcor_1d_array(j_fixed) !the radcor of the center of computational region, used in flux tube model


end subroutine choose_radial_grids


subroutine construct_poloidal_coordinate(r_old,z_old,mpol_old, mpol, theta_new, r_new, z_new) !on a magnetic surface
  use constants,only: p_, two,pi,twopi
  use control_parameters,only: poloidal_angle_type
  use magnetic_field, only : psi_gradient_func, b
  use math, only: arc_length
  implicit none
  integer, intent(in) :: mpol_old, mpol
  real(p_), intent(in) :: r_old(mpol_old), z_old(mpol_old)
  real(p_), intent(in) :: theta_new(mpol)
  real(p_), intent(out) :: r_new(mpol), z_new(mpol)
  real(p_) :: theta_old(mpol_old), dl(mpol_old-1)
  real(p_) :: rmid, zmid, y2(mpol_old)
  integer :: i

  call arc_length(r_old,z_old,mpol_old,dl)
  theta_old(1)=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,mpol_old
        theta_old(i)=theta_old(i-1)+dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'straight-field-line') then
     do i=2,mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)/(rmid*psi_gradient_func(rmid,zmid)) !straight-field-line poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'Boozer') then
     do i=2,mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)*b(rmid,zmid)**2*rmid/psi_gradient_func(rmid,zmid)
     enddo
  else
     stop 'please choose poloidal angle type among equal-arc/equal-volume/Boozer/straight-field-line'
  endif

  theta_old=theta_old*twopi/theta_old(mpol_old) - pi !normalized to the range [-pi:pi]

  !interpolate R  to theta_new gridpoints
  call spline(theta_old, r_old, mpol_old, 2.d30, 2.d30, y2) !prepare the second order derivative needed in the cubic spline interpolation
  do i=2,mpol-1
     call splint(theta_old,r_old,y2,mpol_old,theta_new(i), r_new(i)) !to get R
  enddo
  !interpolate Z  to theta_new gridpoints
  call spline(theta_old, z_old, mpol_old, 2.d30, 2.d30, y2) 
  do i=2,mpol-1
     call splint(theta_old,z_old,y2,mpol_old,theta_new(i), z_new(i))
  enddo

  r_new(1)=r_old(1) !ending points are not included in the above interpolation
  z_new(1)=z_old(1)
  r_new(mpol)=r_new(1)
  z_new(mpol)=z_new(1)

end subroutine construct_poloidal_coordinate



subroutine plasma_volume_of_computational_region()
  !to calculate the spatial volume of computational region
  use constants,only: p_, two, twopi
  use domain_decomposition,only: myid
  use magnetic_coordinates,only: j_low2,j_upp2,mpol,nflux, dradcor, dtheta, toroidal_range, jacobian, &
       & vol !as output
  implicit none
  real(p_) :: dv(mpol, nflux)
  integer :: i, j

  do i=1,mpol
     do j=1,nflux
        dv(i,j)=abs(jacobian(i,j))*dradcor*dtheta*toroidal_range
     enddo
  enddo

  vol=0._p_
  do i=1,mpol-1
     do j=j_low2,j_upp2
        vol = vol + dv(i,j)
     enddo
  enddo

  if(myid.eq.0) write(*,*) 'volume of the toroidal wedge=', vol
  if(myid.eq.0) write(*,*) 'volume of full torus of the computational domain=', vol*twopi/toroidal_range

end subroutine plasma_volume_of_computational_region



subroutine magnetic_coordinates_to_cylindrical_coordinates(theta,radcor,r,z) !given (theta,radcor), return (R,Z)
  use constants,only:p_
  use magnetic_coordinates,only: mpol,nflux,r_mc,z_mc,theta_1d_array,radcor_1d_array
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_),intent(in):: theta,radcor
  real(p_),intent(out):: r,z
  call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,r_mc,theta,radcor,R)  !uniform 1darray is assumed
  call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,z_mc,theta,radcor,Z)  !uniform 1darray is assumed
end subroutine magnetic_coordinates_to_cylindrical_coordinates


function abs_jacobian_func(theta,radcor) result (z) !!used in generating non-uniformly distributed random numbers that satisfied a probability density function proportional to the |jacobian|
  use constants,only: p_
  use magnetic_coordinates,only: mpol,nflux,theta_1d_array,radcor_1d_array,jacobian !as input
  use interpolate_module,only: linear_2d_interpolation
  implicit none
  real(p_)::radcor,theta,z
  real(p_):: abs_jacobian_min,abs_jacobian_max

  call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,abs(jacobian),theta,radcor,z)  !uniform 1d array is assumed

!!$  do i=1,mpol
!!$     do j=1,nflux1
!!$        jshift=j-1+j_low1
!!$        jac1(i,j)=jacobian0(i,jshift)
!!$     enddo
!!$  enddo
!!$  abs_jacobian_max=maxval(jac1)
!abs_jacobian_max=maxval(jacobian0(:,j_low2:j_upp2))
!  abs_jacobian_min=minval(jac1)

  !  z=(z-abs_jacobian_min)/(abs_jacobian_max-abs_jacobian_min) !turns out to be wrong! a subtle bug which I spent much time in finding, the shift is wrong. This bug is due to my misunderstanding about the rejection method, not due to programming mistakes.
    !z=z/abs_jacobian_max
end function abs_jacobian_func
