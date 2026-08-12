module magnetic_coordinates
  use constants, only: p_
  implicit none
  save
  integer :: mpol, nrad, mpol2, mtor 
  !(mpol, mtor, nrad) is number of gridpoints along the (poloidal, toroidal, radial) direction, respectively.
  !mpol2 is number of poloidal gridpoints for perturbed field
  real(p_), dimension(:), allocatable :: xgrid, ygrid, zgrid  !(radial, toroidal, poloidal) grid
  real(p_) :: pfn_inner, pfn_bdry, radial_width
  real(p_) :: dtheta, dradcor, dtor
  real(p_) :: xlow, xupp
  integer :: nsegment !computational region is 1/nsegment of the full torus
  real(p_) :: toroidal_range !toroidal_range=twopi/nsegment

  integer :: itheta0 !poloidal index corresponding to theta=0
  real(p_) :: vol, GSpsi_prime !dGSpsi/dx, x is the normalized poloidal magnetic flux
  real(p_) :: sign_of_jacobian, sign_of_GSpsi_prime
  real(p_), dimension(:,:), allocatable :: r_mc, z_mc !SI units
  real(p_), dimension(:,:), allocatable :: tor_shift_mc, jacobian, abs_jacobian, qhat
  real(p_), dimension(:), allocatable :: tor_shift_left_bdry_minus_one 
  real(p_), dimension(:), allocatable :: jacobian_av
  !GSpsi_array is the poloidal flux appearing in the GS equation (i.e. poiloidal_magnetic_flux/twopi) in SI units (Web/rad)
  real(p_), dimension(:), allocatable :: GSpsi_array, pfn, tfn, minor_r_array, minor_r_prime_array

  real(p_), dimension(:,:),allocatable :: dl_mc, &
       & rth, zth, rpsi, zpsi, &
       & grad_psi_r, grad_psi_z, &
       & grad_theta_r, grad_theta_z, &
       & grad_alpha_r, grad_alpha_z, grad_alpha_phi, &
       & grad_psi, grad_alpha, grad_theta, &
       & grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, grad_alpha_dot_grad_theta, &
       & lapx, lapy ! laplacian of x and y
  real(p_), dimension(:), allocatable :: grad_alpha_r_left_bdry_minus_one, grad_alpha_z_left_bdry_minus_one

end  module magnetic_coordinates

module contour_mod
contains
  subroutine find_contour(psival, xc, zc)
    ! Given a value of the poloidal flux, find the magnetic surface
    use constants, only: p_
    use boundary, only: x_lcfs, z_lcfs, np_lcfs
    use radial_module, only: r_axis, z_axis
    implicit none
    real(p_), intent(in) :: psival
    real(p_), intent(out) :: xc(np_lcfs), zc(np_lcfs)
    real(p_), parameter :: xacc = 1.0d-6 !tolerance used in bi-section root-finder
    real(p_), parameter :: huge = 1d30
    real(p_) :: x1, x2, z1, z2, slope(np_lcfs), slope2(np_lcfs)
    integer :: i

    do i = 1, np_lcfs
       if(x_lcfs(i)-r_axis .ne. 0._p_) then 
          slope(i) = (z_lcfs(i)-z_axis)/(x_lcfs(i)-r_axis) !the slope for function Z=Z(X)
       else
          slope(i) = huge !I use compiler option that catches all erroneous arithmetic operation, I need to avoid dividing by zero
       endif
       if(z_lcfs(i)-z_axis .ne. 0._p_) then
          slope2(i) = (x_lcfs(i)-r_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
       else
          slope2(i) = huge
       endif
    enddo

    do i = 1, np_lcfs-1 
       if(abs(slope(i)).le.1.0_p_) then ! to use Z=Z(X) function, to aviod large slope.
          x1 = r_axis
          x2 = x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
          xc(i) = rtbis(one_dim_psi_func,x1,x2,xacc,r_axis,z_axis,slope(i),psival)
          zc(i) = zfunc(r_axis,z_axis,slope(i),xc(i))
       else ! switch to using X=X(Z) function
          z1 = z_axis
          z2 = z_lcfs(i)
          zc(i) = rtbis(one_dim_psi_func2,z1,z2,xacc,r_axis,z_axis,slope2(i),psival)
          xc(i) = xfunc(r_axis,z_axis,slope2(i),zc(i)) 
       endif
    enddo

    xc(np_lcfs) = xc(1) !i=1 and i=np_lcfs are identical
    zc(np_lcfs) = zc(1) 

  end subroutine find_contour


  function one_dim_psi_func(r_axis,z_axis,slope,psival,x) 
    !poloidal flux as a function of x on a straight line with slope "slope" in poloidal plane
    use constants, only: p_
    use magnetic_field, only: psi_func
    implicit none
    real(p_) :: one_dim_psi_func,x,r_axis,z_axis,slope,psival
    one_dim_psi_func = psi_func(x,zfunc(r_axis,z_axis,slope,x))-psival
  end function one_dim_psi_func


  function zfunc(r_axis,z_axis,slope,x) !straight line Z=Z(x) with slope "slope" in poloidal plane starting from the location of magnetic axis
    use constants,only:p_
    implicit none
    real(p_):: zfunc,x,r_axis,z_axis,slope
    zfunc=z_axis+slope*(x-r_axis)
  end function zfunc


  function one_dim_psi_func2(r_axis,z_axis,slope,psival,z) result(fun_val)
    !poloidal flux as a function of z on a straight line with slope "slope" in poloidal plane
    use constants,only:p_
    use magnetic_field, only : psi_func
    implicit none
    real(p_):: fun_val,z
    real(p_):: r_axis,z_axis,slope,psival
    fun_val = psi_func(xfunc(r_axis,z_axis,slope,z),z)-psival
  end function one_dim_psi_func2


  function xfunc(r_axis,z_axis,slope,z) !straight line X=X(Z) with slope "slope" in poloidal plane starting from the location of magnetic axis
    use constants,only:p_
    implicit none
    real(p_):: xfunc, z
    real(p_)::r_axis,z_axis,slope
    xfunc=r_axis+slope*(z-z_axis)
  end function xfunc


  FUNCTION rtbis(func,x1,x2,xacc,xmaxis,zmaxis,slope,psival)
    !find a root of func by using the bisection method
    use constants,only: p_
    implicit none
    INTEGER, parameter :: JMAX=40
    REAL(p_) rtbis, x1, x2, func, xacc, xmaxis, zmaxis, slope, psival
    external :: func
    INTEGER j
    REAL(p_) dx,f,fmid,xmid
    fmid=func(xmaxis,zmaxis,slope,psival,x2)
    f=   func(xmaxis,zmaxis,slope,psival,x1)
    !      write(*,*) 'f1=', f, 'f2=',fmid
    if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis'

    if(f.lt.0.)then
       rtbis=x1
       dx=x2-x1
    else
       rtbis=x2
       dx=x1-x2
    endif
    do  j=1,JMAX
       dx=dx*.5
       xmid=rtbis+dx
       fmid=func(xmaxis,zmaxis,slope,psival,xmid)
       if(fmid.le.0.)rtbis=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) return
    enddo
    stop 'too many bisections in rtbis'
  end FUNCTION rtbis

end module contour_mod


subroutine construct_magnetic_coordinates()
  use constants, only : p_, zero, twopi,pi
  use control_parameters, only: diagnosis
  use boundary, only: x_lcfs, z_lcfs, np_lcfs
  use radial_module,only: r_axis,z_axis, minor_a
  use magnetic_coordinates, only: mpol, mtor, nsegment, nrad, pfn, tfn, GSpsi_array, &
       & r_mc,z_mc, & !as output
       & zgrid, dtheta, & !as output
       & minor_r_array,minor_r_prime_array, itheta0, dl_mc, & !output
       & toroidal_range, dtor, ygrid, vol !output
  use domain_decomposition, only: myid
  use contour_mod, only: find_contour
  use splines, only: spline3ders
  use math, only:  one_dimensional_derivative, arc_length
  implicit none
  real(p_), allocatable :: r_mag_surf0(:,:), z_mag_surf0(:,:)
  integer :: i, j

  call choose_radial_grids()

  allocate(r_mag_surf0(np_lcfs, nrad))
  allocate(z_mag_surf0(np_lcfs, nrad))
  
  do j = 1, nrad
     call find_contour(GSpsi_array(j), r_mag_surf0(:,j), z_mag_surf0(:,j))
  enddo

  if(myid.eq.0) call diagnostic1()
  minor_a = (maxval(x_lcfs) - minval(x_lcfs))/2
  if(myid==0) write(*,*) 'minor_a=', minor_a
  allocate(minor_r_array(nrad))
  allocate(minor_r_prime_array(nrad))
  do j=1,nrad
     !minor_r_array(j) = r_axis - r_mag_surf0(1,j) !defined on the high-field-side midplane
     minor_r_array(j) = (maxval(r_mag_surf0(:,j))-minval(r_mag_surf0(:,j)))/2
  enddo
  call one_dimensional_derivative(nrad, pfn, minor_r_array, minor_r_prime_array) !to be used as interpolating table
  !call spline3ders(pfn(:), minor_r_array(:), pfn(:), dynew=minor_r_prime_array(:))

  dtheta = twopi/(mpol-1) !poloidal grid spacing for equilibrium
  allocate(zgrid(mpol))
  zgrid = [ (-pi+dtheta*(i-1), i = 1, mpol) ] !equilibrium theta grid
  itheta0 = (mpol+1)/2 !poloidal index corresponding to theta=0 (mpol is assumed odd)

  allocate(r_mc(mpol, nrad))
  allocate(z_mc(mpol, nrad))

  do j = 1, nrad
     call construct_poloidal_coordinate(r_mag_surf0(:,j), z_mag_surf0(:,j), np_lcfs, &
          & mpol, zgrid, r_mc(:,j), z_mc(:,j))
  enddo

  allocate(dl_mc(mpol-1, nrad))
  do j = 1, nrad
     call arc_length(r_mc(:,j), z_mc(:,j), mpol, dl_mc(:,j))
  enddo

  call calculate_metric() 

  toroidal_range = twopi/nsegment
  dtor = toroidal_range/mtor
  allocate(ygrid(mtor+1)) 
  ygrid = [ (zero + dtor*(i-1), i = 1, mtor+1) ]

  call plasma_volume_of_computational_region(vol)
  if ((myid==0) ) call diagnostic2()
  if((myid==0) .and. (diagnosis.eqv..true.)) call diagnostic3()
  deallocate(r_mag_surf0, z_mag_surf0)

contains

  subroutine diagnostic1()
    integer:: u
    open(newunit=u,file='mag_surf_shape0.txt')
    do j=1,nrad
       do i=1,np_lcfs
          write(u,*) r_mag_surf0(i,j), z_mag_surf0(i,j)
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
       do j=1,nrad
          write(u,*) r_mc(i,j),z_mc(i,j),zgrid(i)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)

    open(newunit=u,file='mag_surf_shape.txt')
    do j=1,nrad
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
    do j=1,nrad
       write(u,*) pfn(j), minor_r_array(j),  minor_r_prime_array(j)
    enddo
    close(u)
  end subroutine diagnostic3

end subroutine construct_magnetic_coordinates



subroutine choose_radial_grids()
  use constants, only: p_
  use domain_decomposition, only: myid
  use radial_module, only: psi_axis, psi_lcfs, j_fixed, radcor_fixed
  use magnetic_field, only: radcor_as_func_of_pfn, tfn_func_pfn
  use magnetic_coordinates, only: nrad, pfn_inner, pfn_bdry, &
       & pfn, tfn, GSpsi_array,GSpsi_prime, xgrid,dradcor, & !output
       & xlow, xupp, & !as output
       & radial_width !as output
  implicit none
  real(p_) :: dpfn
  integer :: j, u

  allocate(pfn(nrad))
  allocate(tfn(nrad))  
  allocate(GSpsi_array(nrad))
  allocate(xgrid(nrad))
  dpfn=(pfn_bdry-pfn_inner)/real(nrad-1)
  do j = 1, nrad !select some flux surfaces (labeld by GSpsi_array)
     pfn(j) = pfn_inner +dpfn*(j-1)
     GSpsi_array(j) = psi_axis + pfn(j)*(psi_lcfs-psi_axis) !GSpsi=Aphi*R, this array is used in finding magnetic surfaces
     xgrid(j) = radcor_as_func_of_pfn(pfn(j))
  enddo

  xlow = xgrid(1)
  xupp = xgrid(nrad)

  radial_width = xgrid(nrad) - xgrid(1)
  dradcor = xgrid(2) - xgrid(1) !radial grid interval
  GSpsi_prime = psi_lcfs - psi_axis !dGSpsi/dx, x is the normalized poloidal magnetic flux

  j_fixed = nrad/2
  radcor_fixed = xgrid(j_fixed) !the radcor of the center of computational region

  do j = 1, nrad
     tfn(j) = tfn_func_pfn(pfn(j))
  enddo

  if(myid==0) then
     open(newunit=u, file='xgrid.txt')
     do j = 1, nrad
        write(u, *) xgrid(j), tfn(j)
     enddo
     close(u)
  endif

end subroutine choose_radial_grids


subroutine construct_poloidal_coordinate(r_old, z_old, mpol_old, mpol, theta_new, r_new, z_new) !on a magnetic surface
  use constants,only: p_, two,pi,twopi
  use control_parameters,only: poloidal_angle_type
  use magnetic_field, only : psi_gradient_func, b
  use math, only: arc_length
  use interpolate_module, only: linear_1d_interpolate_nonuniform
  use splines, only: spline3
  implicit none
  integer, intent(in) :: mpol_old, mpol
  real(p_), intent(in) :: r_old(mpol_old), z_old(mpol_old)
  real(p_), intent(in) :: theta_new(mpol)
  real(p_), intent(out) :: r_new(mpol), z_new(mpol)
  real(p_) :: theta_old(mpol_old), dl(mpol_old-1)
  real(p_) :: rmid, zmid, y2(mpol_old)
  integer :: i

  call arc_length(r_old, z_old, mpol_old, dl)
  
  theta_old(1) = 0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i = 2, mpol_old
        theta_old(i) = theta_old(i-1) + dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i = 2, mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'straight-field-line') then
     do i = 2, mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)/(rmid*psi_gradient_func(rmid,zmid)) !straight-field-line poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'Boozer') then
     do i = 2, mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)*b(rmid,zmid)**2*rmid/psi_gradient_func(rmid,zmid)
     enddo
  else
     stop 'please choose poloidal angle type among equal-arc/equal-volume/Boozer/straight-field-line'
  endif

  !normalize and shift to the range [-pi:pi]
  theta_old = theta_old*twopi/theta_old(mpol_old) - pi 

  !interpolate R  to theta_new gridpoints
!!$  call spline(theta_old, r_old, mpol_old, 2.d30, 2.d30, y2) !prepare the second order derivative 
!!$  do i = 2, mpol-1
!!$     call splint(theta_old, r_old, y2, mpol_old, theta_new(i), r_new(i))
!!$  enddo
  do i = 2, mpol-1
     call linear_1d_interpolate_nonuniform(mpol_old, theta_old, r_old, theta_new(i), r_new(i))
  enddo
  
  !interpolate Z  to theta_new gridpoints
!!$  call spline(theta_old, z_old, mpol_old, 2.d30, 2.d30, y2) !prepare the second order derivative 
!!$  do i = 2, mpol-1
!!$     call splint(theta_old, z_old, y2, mpol_old, theta_new(i), z_new(i))
!!$  enddo
  do i = 2, mpol-1
     call linear_1d_interpolate_nonuniform(mpol_old, theta_old, z_old, theta_new(i), z_new(i))
  enddo

!!$    r_new(2:mpol-1) = spline3(theta_old, r_old(:), theta_new(2:mpol-1))
!!$    z_new(2:mpol-1) = spline3(theta_old, z_old(:), theta_new(2:mpol-1))

  r_new(1) = r_old(1) !ending points are not included in the above interpolation
  z_new(1) = z_old(1)
  r_new(mpol) = r_new(1)
  z_new(mpol) = z_new(1)
end subroutine construct_poloidal_coordinate


subroutine plasma_volume_of_computational_region(vol)
  use constants,only: p_, two, twopi
  use domain_decomposition,only: myid
  use magnetic_coordinates,only: mpol,nrad, dradcor, dtheta, toroidal_range, jacobian
  implicit none
  real(p_), intent(out) :: vol
  real(p_) :: dv(mpol, nrad)
  integer :: i, j

  do i=1,mpol
     do j=1,nrad
        dv(i,j)=abs(jacobian(i,j))*dradcor*dtheta*toroidal_range
     enddo
  enddo

  vol=0._p_
  do i=1,mpol-1
     do j=1,nrad
        vol = vol + dv(i,j)
     enddo
  enddo

  if(myid==0) write(*,*) 'Volume of the toroidal wedge (m^3):',  vol
  if(myid==0) write(*,*) 'Volume of full torus of the computational domain (m^3):', vol*twopi/toroidal_range

end subroutine plasma_volume_of_computational_region


module misc

contains

subroutine magnetic_coordinates_to_cylindrical_coordinates(theta,radcor,r,z) !given (theta,radcor), return (R,Z)
  use constants, only: p_
  use magnetic_coordinates, only: mpol,nrad,r_mc,z_mc,zgrid,xgrid
  use interpolate_module, only: linear_2d_interpolate
  implicit none
  real(p_), intent(in) :: theta, radcor
  real(p_), intent(out) :: r,z
  call linear_2d_interpolate(mpol,nrad, zgrid,xgrid, r_mc, theta,radcor,R)  !uniform 1darray is assumed
  call linear_2d_interpolate(mpol,nrad, zgrid,xgrid, z_mc, theta,radcor,Z)  !uniform 1darray is assumed
end subroutine magnetic_coordinates_to_cylindrical_coordinates


function abs_jacobian_func(theta, radcor) result (z) !!used in generating non-uniformly distributed random numbers that satisfied a probability density function proportional to the |jacobian|
  use constants, only: p_
  use magnetic_coordinates, only: mpol,nrad,zgrid,xgrid, abs_jacobian !as input
  use interpolate_module, only: linear_2d_interpolate
  implicit none
  real(p_) :: radcor,theta,z
  !real(p_) :: abs_jacobian_min,abs_jacobian_max

  call linear_2d_interpolate(mpol,nrad, zgrid, xgrid, abs_jacobian,theta,radcor,z)  !uniform 1d array is assumed

!!$  do i=1,mpol
!!$     do j=1,nrad
!!$        jshift=j
!!$        jac1(i,j)=jacobian0(i,jshift)
!!$     enddo
!!$  enddo
!!$  abs_jacobian_max=maxval(jac1)
  !abs_jacobian_max=maxval(jacobian0(:,1:nrad))
  !  abs_jacobian_min=minval(jac1)

  !  z=(z-abs_jacobian_min)/(abs_jacobian_max-abs_jacobian_min) !turns out to be wrong! a subtle bug which I spent much time in finding, the shift is wrong. This bug is due to my misunderstanding about the rejection method, not due to programming mistakes.
  !z=z/abs_jacobian_max
end function abs_jacobian_func


subroutine calculate_dvol(m, dvol) !one dtheta2 cell can include multiple equilibrium dtheta
  use constants, only: p_
  use domain_decomposition, only: ipol_eq, dtheta2
  use magnetic_coordinates, only: nrad, mpol, nrad, abs_jacobian, dradcor, dtor
  implicit none
  integer, intent(in) :: m ! dtheta2 = m*dtheta
  real(p_), intent(out), allocatable :: dvol(:)
  real(p_) :: jac(-m:mpol+m, nrad), jac0
  integer :: i, i1, i2, s, j, jeq

  allocate(dvol(nrad))

  jac(1:mpol, :) = abs_jacobian(:,:)
  do i = 0, -m, -1
     jac(i,:)= abs_jacobian(mpol-1+i,:)
  enddo
  do i = mpol+1, mpol+m
     jac(i,:)= abs_jacobian(i+1-mpol,:)
  enddo

  if(mod(m,2) == 0) then
     s = m/2
     i1 = ipol_eq - s
     i2 = ipol_eq + s
     do j = 1, nrad
        jeq = j
        jac0 = (sum(jac(i1+1:i2-1, jeq)) + 0.5*(jac(i1, jeq) + jac(i2, jeq)))/m
        dvol(j) = jac0*dradcor*dtheta2*dtor
     enddo
  else
     s = (m-1)/2
     i1 = ipol_eq - s
     i2 = ipol_eq + s
     do j = 1, nrad
        jeq = j
        jac0 = sum(jac(i1:i2, jeq))/m
        dvol(j) = jac0*dradcor*dtheta2*dtor
     enddo
  endif
end subroutine calculate_dvol

end module misc


module table_in_mc !used in computing guiding-center drift
  use constants, only: p_
  implicit none
  save
  real(p_), dimension(:,:), allocatable :: br_mc, bz_mc, bphi_mc, &
       & bp_mc, b_mc, bdgxcgy, bdgxcgz, bdgycgz, &
       & w1,w2,w3,w4,w5,w5p,w6,w7,w8,w8p,w9,w10,w12,w13
   !here bdgxcgy is B0_dot_grad_x_cross_grad_y. Similar names for other cooordinate combination.
contains

  subroutine prepare_table_in_mc()
    use constants, only: one
    use domain_decomposition, only: myid
    use magnetic_coordinates, only: mpol, nrad, r_mc,z_mc, jacobian, xgrid, &
         &  GSpsi_prime, grad_psi, grad_alpha, &
         & grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, grad_alpha_dot_grad_theta, &
         & grad_psi_r, grad_psi_z, grad_theta_r, grad_theta_z, &
         & grad_alpha_r, grad_alpha_z
    use magnetic_field, only: br, bz, bphi, b, &
         & b_r, b_z, unitbr_z, unitbz_r, unitbphi_r,unitbphi_z
    implicit none

    real(p_) :: brval,bzval,bphival,bval
    real(p_) :: b_rval,b_zval
    real(p_) :: dradial_dr_func,dradial_dz_func !function names
    real(p_) :: dradial_dr_func2,dradial_dz_func2 !function names
    real(p_) :: dtheta_dr_func,dtheta_dz_func !function names
    real(p_) :: ddelta_dr_func,ddelta_dz_func,ddelta_dr_lsf_midplane_twopi,ddelta_dz_lsf_midplane_twopi !function names
    real(p_) :: dradial_dr_val,dradial_dz_val
    real(p_) :: dtheta_dr_val,dtheta_dz_val
    real(p_) :: ddelta_dr_val,ddelta_dz_val
    real(p_) :: unitbr,unitbz,unitbphi
    real(p_) :: curl_unitb_rcomp,curl_unitb_zcomp,curl_unitb_phicomp
    real(p_) :: unitb_dot_curl_unitb
    real(p_) :: r,z,radcor
    real(p_) :: grad_psi_val,grad_alpha_val, grad_psi_dot_grad_alpha_val
    real(p_) :: grad_psi_dot_grad_theta_val,grad_alpha_dot_grad_theta_val
    real(p_) :: dalpha_dr_val,dalpha_dz_val,dalpha_dphi_val
    integer :: i, j

    allocate(w1(mpol,nrad))
    allocate(w2(mpol,nrad))
    allocate(w3(mpol,nrad))
    allocate(w4(mpol,nrad))
    allocate(w5(mpol,nrad))
    allocate(w5p(mpol,nrad))
    allocate(w6(mpol,nrad))
    allocate(w7(mpol,nrad))
    allocate(w8(mpol,nrad))
    allocate(w8p(mpol,nrad))
    allocate(w9(mpol,nrad))
    allocate(w10(mpol,nrad))
    allocate(w12(mpol,nrad))
    allocate(w13(mpol,nrad))
    allocate(bdgxcgz(mpol,nrad))
    allocate(bdgycgz(mpol,nrad))
    allocate(bdgxcgy(mpol,nrad))
    allocate(b_mc(mpol,nrad))
    allocate(br_mc(mpol,nrad))
    allocate(bz_mc(mpol,nrad))
    allocate(bphi_mc(mpol,nrad))
    allocate(bp_mc(mpol,nrad))


    do j = 1, nrad
       radcor = xgrid(j)
       do i=1,mpol
          r = r_mc(i,j)
          z = z_mc(i,j)
          brval=br(r,z)
          bzval=bz(r,z)
          bphival=bphi(r,z)
          bval=sqrt(brval**2+bzval**2+bphival**2)
          b_zval=b_z(r,z)
          b_rval=b_r(r,z)
          unitbr=brval/bval
          unitbz=bzval/bval
          unitbphi=bphival/bval
          curl_unitb_rcomp = -unitbphi_z(r,z)
          curl_unitb_phicomp = unitbr_z(r,z)-unitbz_r(r,z)
          curl_unitb_zcomp = unitbphi_r(r,z)+unitbphi/r

          unitb_dot_curl_unitb=unitbr*curl_unitb_rcomp +unitbphi*curl_unitb_phicomp &
               & +unitbz*curl_unitb_zcomp

          dradial_dr_val = grad_psi_r(i,j)
          dradial_dz_val = grad_psi_z(i,j) 
          dtheta_dr_val = grad_theta_r(i,j)
          dtheta_dz_val = grad_theta_z(i,j)
          ddelta_dr_val = -grad_alpha_r(i,j)
          ddelta_dz_val = -grad_alpha_z(i,j)

          dalpha_dr_val = grad_alpha_r(i,j)
          dalpha_dz_val = grad_alpha_z(i,j)
          dalpha_dphi_val=one

          !        if(i.eq.mpol)    ddelta_dr_val=ddelta_dr_lsf_midplane_twopi(r) !at theta=twopi cut
          !        if(i.eq.mpol)    ddelta_dz_val=ddelta_dz_lsf_midplane_twopi(r) !at theta=twopi cut

!!$!write(*,*) 'i,j=',dradial_dr_val, dradial_dr_func2(r,z),dradial_dz_val, dradial_dz_func2(r,z)
          w1(i,j) = unitb_dot_curl_unitb/bval
          !w1(i,j) = 0 ! approximation, tested, no difference from the above line
          w2(i,j)=unitbr*dtheta_dr_val+unitbz*dtheta_dz_val
          !      w2(i,j)=-(psi_lcfs-psi_axis)/(bval*jacobian(i,j))
          w3(i,j)=curl_unitb_rcomp/bval*dradial_dr_val+curl_unitb_zcomp/bval*dradial_dz_val
          w4(i,j)=curl_unitb_rcomp/bval*dtheta_dr_val+curl_unitb_zcomp/bval*dtheta_dz_val
          w5(i,j)=curl_unitb_phicomp/(bval*r)-curl_unitb_rcomp/bval*ddelta_dr_val-curl_unitb_zcomp/bval*ddelta_dz_val
          w5p(i,j)=curl_unitb_phicomp/(bval*r)
          w6(i,j)=(bphival*b_zval)*dradial_dr_val+(-bphival*b_rval)*dradial_dz_val
          w6(i,j)=w6(i,j)/bval**2
          w7(i,j)=(bphival*b_zval)*dtheta_dr_val+(-bphival*b_rval)*dtheta_dz_val
          w7(i,j)=w7(i,j)/bval**2
          w8(i,j)=(bzval*b_rval-brval*b_zval)/r-(bphival*b_zval)*ddelta_dr_val-(-bphival*b_rval)*ddelta_dz_val
          w8(i,j)=w8(i,j)/bval**2
          w8p(i,j)=(bzval*b_rval-brval*b_zval)/(r*bval**2)
          w9(i,j)=unitbr*b_rval+unitbz*b_zval
          w10(i,j)=(b_rval*curl_unitb_rcomp+b_zval*curl_unitb_zcomp)/bval

!!$        grad_psi(i,j)=sqrt(dradial_dr_val**2+dradial_dz_val**2) !now calculated in another subroutine
!!$        grad_alpha(i,j)=sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)
!!$        grad_psi_dot_grad_alpha(i,j)=dradial_dr_val*ddelta_dr_val+dradial_dz_val*ddelta_dz_val

          b_mc(i,j)=bval
          br_mc(i,j)=brval
          bz_mc(i,j)=bzval
          bphi_mc(i,j)=bphival
          bp_mc(i,j)=sqrt(brval**2+bzval**2)

          grad_psi_val=grad_psi(i,j)
          grad_alpha_val=grad_alpha(i,j)
          grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha(i,j)
          grad_psi_dot_grad_theta_val=grad_psi_dot_grad_theta(i,j)
          grad_alpha_dot_grad_theta_val=grad_alpha_dot_grad_theta(i,j)

          w12(i,j)=GSpsi_prime/bval*(grad_psi_dot_grad_theta_val*grad_psi_dot_grad_alpha_val &
               & -grad_psi_val**2*grad_alpha_dot_grad_theta_val)
          w13(i,j)=GSpsi_prime/bval*(grad_psi_dot_grad_theta_val*grad_alpha_val**2 &
               & -grad_psi_dot_grad_alpha_val*grad_alpha_dot_grad_theta_val)

          bdgxcgz(i,j) = bphival*(dradial_dz_val*dtheta_dr_val-dradial_dr_val*dtheta_dz_val)
          bdgycgz(i,j) = brval*dtheta_dz_val/r &
               & + bphival*(dalpha_dz_val*dtheta_dr_val - dalpha_dr_val*dtheta_dz_val) &
               & + bzval*(-dtheta_dr_val)/r
          bdgxcgy(i,j) = bval**2/GSpsi_prime
       enddo

    enddo

    if(myid.eq.0) call diagnostic()
  contains
    subroutine diagnostic()
      character(8)::filename
      integer:: file_unit,u
      open(newunit=u,file='grad_theta_alpha.txt')
      do j=1,nrad
         radcor=xgrid(j)
         do i=1,mpol
!!$          write(u,*) r_mc(i,j), z_mc(i,j), dtheta_dr_val,dtheta_dz_val,&
!!$               & sqrt(dtheta_dz_val**2+dtheta_dr_val**2), ddelta_dr_val,  ddelta_dz_val,&
!!$               sqrt(ddelta_dr_val**2+ddelta_dz_val**2), sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)
!!$          write(u,*) r_mc(i,j), z_mc(i,j),jacobian(i,j),&
!!$               & -w2(i,j)

            !write(u,*) r_mc(i,j), z_mc(i,j),w12(i,j),w13(i,j),w5(i,j)
            write(u,*) r_mc(i,j), z_mc(i,j), w5(i,j), w8(i,j)

            !               & -b_mc(i,j)/GSpsi_prime*w2(i,j)


         enddo
         write(u,*) 
      enddo
      close(u)
!!$  write(filename,'(a4,i4.4)') 'metr',myid
!!$  file_unit=myid+311
!!$  open(file_unit,file=filename)
!!$  do i=1,mpol
!!$     do j=1,nrad
!!$        write(file_unit,'(2i8.4,3(1pe14.5))')  i,j,grad_alpha(i,j),grad_psi_dot_grad_alpha(i,j), grad_psi(i,j)
!!$     enddo
!!$     write(file_unit,*)
!!$  enddo
!!$  close(file_unit)

      ! write(*,*) 'min max w1-4=', maxval(w1),minval(w1),maxval(w2),minval(w2),maxval(w3),minval(w3),maxval(w4),minval(w4)
      !write(*,*) 'min max w5-8=', maxval(w5),minval(w5),maxval(w6),minval(w6),maxval(w7),minval(w7),maxval(w8),minval(w8)
      !write(*,*) 'min max w9-10=', maxval(w9),minval(w9),maxval(w10),minval(w10)

    end subroutine diagnostic

  end subroutine prepare_table_in_mc

end module table_in_mc



module func_in_mc
contains
  function b_mc_func(theta, radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: b_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,b_mc,theta,radcor,f)
  end function b_mc_func

  function br_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: br_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,br_mc,theta,radcor,f)
  end function br_mc_func

  function bz_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: bz_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bz_mc,theta,radcor,f)
  end function bz_mc_func


  function bphi_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: bphi_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bphi_mc,theta,radcor,f)
  end function bphi_mc_func

pure real(p_)  function tor_shift_func(theta, radcor) result(z)
    use constants, only: p_
    use magnetic_coordinates, only: mpol, nrad, zgrid, xgrid, tor_shift_mc
    use interpolate_module, only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: theta, radcor
    call linear_2d_interpolate(mpol, nrad, zgrid, xgrid, tor_shift_mc, theta, radcor, z) !interpolate in magnetic coordinates
  end function tor_shift_func


pure real(p_)   function grad_psi_func(theta,radcor) result(f)
    use constants, only: p_
    use magnetic_coordinates, only: mpol, nrad, zgrid, xgrid, grad_psi
    use interpolate_module, only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: theta, radcor
    call linear_2d_interpolate(mpol, nrad, zgrid, xgrid, grad_psi, theta, radcor, f)
  end function grad_psi_func

  pure real(p_) function minor_r_radcor(radcor) result (z)
    use constants,only: p_, two,twopi
    use magnetic_coordinates,only: nrad, xgrid, minor_r_array
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_),intent(in):: radcor

    call linear_1d_interpolate(nrad, xgrid, minor_r_array, radcor, z)  
  end function minor_r_radcor

  pure real(p_) function radcor_minor_r(minor_r) result (z)
    use constants,only: p_
    use magnetic_coordinates,only: nrad,xgrid,minor_r_array
    use interpolate_module,only: linear_1d_interpolate_nonuniform
    implicit none
    real(p_), intent(in):: minor_r

    call linear_1d_interpolate_nonuniform(nrad, minor_r_array, xgrid, minor_r, z)

  end function radcor_minor_r


  pure real(p_) function minor_r_prime(radcor) result (z) !derivative of minor_r with respect to the radial coordinate
    use constants, only: p_, two
    use magnetic_coordinates, only: nrad,xgrid,minor_r_prime_array
    use interpolate_module, only: linear_1d_interpolate
    implicit none
    real(p_), intent(in) :: radcor

    call linear_1d_interpolate(nrad, xgrid, minor_r_prime_array, radcor, z)  
  end function minor_r_prime

  
end module func_in_mc

module calculate_toroidal_shift_module
contains
  subroutine calculate_toroidal_shift(psival,x_contour,z_contour,np_lcfs,end_i,r,z,tor_shift) 
    use constants, only: p_, one_half
    use magnetic_coordinates, only: sign_of_jacobian,sign_of_GSpsi_prime, itheta0
    use radial_module, only: psi_axis,psi_lcfs
    use domain_decomposition, only: myid
    use magnetic_field, only : psi_z_func,psi_r_func, g_func !toroidal field function
    implicit none
    real(p_), intent(in) :: psival
    integer, intent(in) :: np_lcfs,end_i
    real(p_) :: x_contour(np_lcfs),z_contour(np_lcfs)
    real(p_), intent(in) :: r,z
    real(p_), intent(out) :: tor_shift
    real(p_) :: x_mid,z_mid,gx0,dl, g_value
    real(p_) :: pfn,mr,costh,sinth,r0
    integer :: i

    g_value = g_func(psival)

    tor_shift = 0._p_
    if (end_i .lt. itheta0) then 
       do i = itheta0-1, end_i+1, -1
          x_mid=(x_contour(i)+x_contour(i-1))*one_half
          z_mid=(z_contour(i)+z_contour(i-1))*one_half
          gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
          dl=-sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
          tor_shift=tor_shift+g_value/(x_mid*gx0)*dl
       enddo
       x_mid=(x_contour(end_i+1)+r)*one_half
       z_mid=(z_contour(end_i+1)+z)*one_half
       gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
       dl=-sqrt((x_contour(end_i+1)-r)**2+(z_contour(end_i+1)-z)**2)
       tor_shift=tor_shift+g_value/(x_mid*gx0)*dl

    elseif(end_i.ge.itheta0) then
       x_contour(end_i+1)=r !replace No. end_i+1 point by the given point (r,z)
       z_contour(end_i+1)=z
       do i=itheta0, end_i
          x_mid=(x_contour(i)+x_contour(i+1))*one_half
          z_mid=(z_contour(i)+z_contour(i+1))*one_half
          gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
          dl=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
          tor_shift=tor_shift+g_value/(x_mid*gx0)*dl
       enddo
    endif
    !--for testing----analytical formula, for concentric circular configuration
!!$  r0=1.32_p_
!!$  mr=sqrt((x_contour(end_i+1)-r0)**2+(z_contour(end_i+1)-0._p_)**2)
!!$  costh=(x_contour(end_i+1)-r0)/mr
!!$  sinth=z_contour(end_i+1)/mr
!!$  pfn=(psival-psi_axis)/(psi_lcfs-psi_axis)
!!$  tor_shift=2*qfunc(pfn)*atan((r0-mr)/sqrt(r0**2-mr**2)*sinth/(costh+1))
    !----for tesing----------------------

!!$  do i=2,end_i+1
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     tor_shift=tor_shift+g_value/(x_mid*gx0(i-1))*dl(i-1)
!!$  enddo
    tor_shift=tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
  end subroutine calculate_toroidal_shift


  pure subroutine calculate_toroidal_shift2(j, g_value, x_contour,z_contour, jacob, m, &
       & tor_shift, tor_shift_left_bdry_minus_one)
    !cf. calculate_toroidal_shift, the differences are (1) g_value instead of psival is provided as input;
    !(2) compute a series of tor_shift on a magnetic surface (rather than a single value). The theta range is [-pi:pi]
    use constants, only: p_, one_half, two
    use magnetic_coordinates, only: dtheta, sign_of_jacobian,sign_of_GSpsi_prime,itheta0, GSpsi_array, GSpsi_prime, qhat
    use domain_decomposition, only: dtheta2, multi_eq_cells
    use magnetic_field, only: psi_z_func, psi_r_func
    implicit none
    integer,intent(in) :: j, m
    real(p_),intent(in) :: g_value
    real(p_),intent(in) :: x_contour(m), z_contour(m), jacob(m)
    real(p_),intent(out) :: tor_shift(m), tor_shift_left_bdry_minus_one
    real(p_) :: gx0, dl, sum, x_mid, z_mid, jac_mid, qhat_mid
    integer :: i

!!$  tor_shift(itheta0) = 0._p_
!!$  do i=itheta0+1, m !for theta in (0:pi]
!!$     x_mid=(x_contour(i-1)+x_contour(i))/two
!!$     z_mid=(z_contour(i-1)+z_contour(i))/two
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
!!$     tor_shift(i) = tor_shift(i-1) + g_value/(x_mid*gx0)*dl
!!$  enddo
!!$
!!$  do i=itheta0-1, 1, -1 !for theta in (0:-pi]
!!$     x_mid=(x_contour(i+1)+x_contour(i))/two
!!$     z_mid=(z_contour(i+1)+z_contour(i))/two
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2) 
!!$     dl=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
!!$     tor_shift(i) = tor_shift(i+1) - g_value/(x_mid*gx0)*dl
!!$  enddo
!!$  tor_shift = tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign

    !-------------another method------------:
    tor_shift(itheta0) = 0._p_
    do i = itheta0 + 1, m !for theta in (0:pi]
       x_mid = (x_contour(i-1) + x_contour(i))/two
       jac_mid = (jacob(i-1) + jacob(i))/two
       qhat_mid = (qhat(i-1,j) + qhat(i,j))/two
       !tor_shift(i) = tor_shift(i-1) - g_value/(x_mid**2)*jac_mid/GSpsi_prime*dtheta
       tor_shift(i) = tor_shift(i-1) + qhat_mid*dtheta
    enddo

    do i = itheta0 - 1, 1, -1 !for theta in (0:-pi]
       x_mid = (x_contour(i+1)+x_contour(i))/two
       jac_mid = (jacob(i+1)+jacob(i))/two
       qhat_mid = (qhat(i+1,j) + qhat(i,j))/two
       !tor_shift(i) = tor_shift(i+1) + g_value/(x_mid**2)*jac_mid/GSpsi_prime*dtheta
       tor_shift(i) = tor_shift(i+1) - qhat_mid*dtheta
    enddo
    !---------------------------------------
    !another way of calculating tor_shift
!!$ tor_shift(1)=0._p_
!!$  do i=1,m-1 
!!$     x_mid=(x_contour(i+1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i+1)+z_contour(i))*one_half
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i+1)-x_contour(i))**2+(z_contour(i+1)-z_contour(i))**2)
!!$     tor_shift(i+1)=tor_shift(i)+g_value/(x_mid*gx0)*dl
!!$  enddo
!!$  tor_shift=tor_shift-tor_shift(itheta0)

!!$  tor_shift(1)=0._p_
!!$  do i=2,m
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i-1)+z_contour(i))*one_half
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
!!$     tor_shift(i)=tor_shift(i-1)+g_value/(x_mid*gx0)*dl
!!$  enddo
    !tor_shift = tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign

    tor_shift_left_bdry_minus_one = tor_shift(1) - (tor_shift(m) - tor_shift(m-multi_eq_cells))

  end subroutine calculate_toroidal_shift2


  subroutine calculate_toroidal_shift_at_theta_cut(psival,x_contour,z_contour,np,tor_shift_a,tor_shift_b)
    use constants,only: p_, one_half
    use magnetic_coordinates,only: sign_of_jacobian, sign_of_GSpsi_prime, itheta0
    use domain_decomposition,only:myid
    use magnetic_field, only : psi_z_func,psi_r_func, g_func
    implicit none
    real(p_),intent(in):: psival
    integer,intent(in):: np
    real(p_),intent(in):: x_contour(np),z_contour(np)
    real(p_),intent(out):: tor_shift_a, tor_shift_b
    real(p_):: x_mid,z_mid,gx0,dl, g_value
    integer:: i

    g_value=g_func(psival)

    tor_shift_b=0._p_
    do i=itheta0+1,np ! for location near the cut above the midplane
       x_mid=(x_contour(i)+x_contour(i-1))*one_half
       z_mid=(z_contour(i)+z_contour(i-1))*one_half
       gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
       dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
       tor_shift_b=tor_shift_b+g_value/(x_mid*gx0)*dl
    enddo

    tor_shift_a=0._p_
    do i=itheta0,2,-1 ! for location near the cut below the midplane
       x_mid=(x_contour(i)+x_contour(i-1))*one_half
       z_mid=(z_contour(i)+z_contour(i-1))*one_half
       gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
       dl=-sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
       tor_shift_a=tor_shift_a+g_value/(x_mid*gx0)*dl
    enddo

    tor_shift_a=tor_shift_a*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
    tor_shift_b=tor_shift_b*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
  end subroutine calculate_toroidal_shift_at_theta_cut
end module calculate_toroidal_shift_module


subroutine calculate_metric()
  use constants, only: p_, twopi
  use magnetic_coordinates,only: mpol,nrad,r_mc,z_mc,dtheta,dradcor, &
       & tor_shift_mc, tor_shift_left_bdry_minus_one, & !output
       & qhat, jacobian, & !as output
       & pfn_inner,pfn_bdry,GSpsi_array, &
       & zgrid, xgrid, &
       & sign_of_jacobian,sign_of_GSpsi_prime, GSpsi_prime , grad_psi 

  use radial_module, only: psi_lcfs,psi_axis, npsi, sign_bphi, qpsi
  use radial_module, only: q_with_sign, qrad !as output
  use magnetic_field, only : g_func, qfunc
  use calculate_toroidal_shift_module, only : calculate_toroidal_shift2
  use domain_decomposition, only: myid
  use control_parameters, only: diagnosis
  implicit none
  integer :: i, j, ierr
  real(p_) :: g

  call calculate_gradients_of_psi_and_theta() !computing gradients at mc gridpoints

  if(myid==0) write(*,*) 'max_jacobian=', maxval(jacobian), 'min_jacobian=', minval(jacobian)

  sign_of_jacobian=sign(1._p_,jacobian(mpol/3,nrad/3))
  sign_of_GSpsi_prime=sign(1._p_,GSpsi_prime) !radial coordinator is assumed to be pfn

  allocate(q_with_sign(npsi))
  q_with_sign = abs(qpsi)*(-sign_bphi)*sign_of_jacobian*sign_of_GSpsi_prime !include the correct sign
  ! after this call, function "qfunc" is ready to be used
  allocate(qrad(nrad))
  do i = 1, nrad
     qrad(i) = qfunc(xgrid(i))
  enddo
  if(myid==0) write(*,*) 'q_inner=', qrad(1), 'q_bdry=', qrad(nrad)

  allocate(tor_shift_mc(mpol, nrad))
  allocate(tor_shift_left_bdry_minus_one(nrad)) 

  ! local safety factor
  allocate(qhat(mpol,nrad))
  do j = 1, nrad
     g = g_func(GSpsi_array(j))
     do i = 1, mpol
        qhat(i,j) = -g/r_mc(i,j)**2*jacobian(i,j)/GSpsi_prime
     enddo
  enddo

  do j = 1, nrad
     g = g_func(GSpsi_array(j))
     call calculate_toroidal_shift2(j, g, r_mc(:,j), z_mc(:,j), jacobian(:,j), mpol, &
          & tor_shift_mc(:,j), tor_shift_left_bdry_minus_one(j))

  enddo

!!$ if(myid.eq.0) call draw_grids_on_theta_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc)
!!$ if(myid.eq.0) call draw_alpha_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc)
!!$ if(myid.eq.0) call draw_alpha_contours_on_a_magnetic_surface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !i.e., magnetic field lines

  call calculate_gradients_of_generalized_toroidal_angle()

  if(myid==0) then
     call diagnostic1()
     call diagnostic2()
  endif

contains
  subroutine diagnostic1()
    integer:: i,u
    write(*,*) 'maximum of tor_shift_mc=', maxval(tor_shift_mc)
    write(*,*) 'minimum of tor_shift_mc=', minval(tor_shift_mc)
    open(newunit=u,file='qhat.txt')
    do j = 1, nrad
       do i = 1, mpol
          write(u,*) r_mc(i,j), z_mc(i,j), qhat(i,j), tor_shift_mc(i,j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine diagnostic1

  subroutine diagnostic2()
    use magnetic_coordinates, only : pfn, dl_mc
    use func_in_mc, only: minor_r_radcor
    integer :: i, u
    real(p_) :: s, q2
!!$     do j=1,nrad
!!$        write(*,*) j,xgrid(j),minor_r_radcor(xgrid(j)),r_mag_surf0(1,j)-r_axis,&
!!$             & minor_r_prime(xgrid(j))
!!$     enddo
    open(newunit=u,file="q2.txt")
    do j = 1, nrad
       s = 0
       do i = 1, mpol-1
          s = s + g_func(GSpsi_array(j))/(grad_psi(i,j)*(psi_lcfs-psi_axis)*r_mc(i,j))*dl_mc(i,j)
       enddo
       q2 = s/twopi
       write(u,'(20(ES16.4E2))') pfn(j), minor_r_radcor(pfn(j)), q2, qrad(j)
    enddo
    close(u)
  end subroutine diagnostic2

end subroutine calculate_metric


subroutine calculate_gradients_of_psi_and_theta()
  use constants, only: p_, zero, one, two
  use magnetic_coordinates, only: m=>mpol, n=>nrad, r=>r_mc, z=>z_mc, dtheta, dradcor, xgrid, zgrid, &
       & jacobian, abs_jacobian, jacobian_av, & !as output
       & rth, zth, rpsi, zpsi, grad_psi, grad_theta, grad_psi_dot_grad_theta, & !as output
       & grad_psi_r, grad_psi_z, grad_theta_r, grad_theta_z, lapx !as output
  use domain_decomposition, only: myid
  use splines, only: spline3ders
  use math, only: one_dimensional_derivative
  !use magnetic_field, only: minor_r_prime, minor_r_radcor
  implicit none
  integer :: i, j, u
  real(p_) ::  minor_r_prime_val
  real(p_), allocatable :: a(:,:), b(:,:), apsi(:,:), bth(:,:)

  allocate(jacobian(m,n))
  allocate(abs_jacobian(m,n))
  allocate(jacobian_av(n))
  allocate(grad_psi(m,n))
  allocate(grad_theta(m,n))
  allocate(grad_psi_dot_grad_theta(m,n))
  allocate(grad_theta_r(m,n))
  allocate(grad_theta_z(m,n))
  allocate(grad_psi_r(m,n))
  allocate(grad_psi_z(m,n))
  allocate(rpsi(m,n))
  allocate(zpsi(m,n))
  allocate(rth(m,n))
  allocate(zth(m,n))
  allocate(lapx(m,n))

  call partial_derivative_in_mc(m,n,r,z,dtheta,dradcor,rpsi,rth,zpsi,zth,jacobian)
  !call partial_derivative_in_mc2(m,n,r,z,dtheta,dradcor,rpsi,rth,zpsi,zth,jacobian) 

  abs_jacobian = abs(jacobian)

  grad_psi = r/abs_jacobian*sqrt(zth**2+rth**2)
  grad_theta = r/abs_jacobian*sqrt(zpsi**2+rpsi**2)
  grad_psi_dot_grad_theta = -(r/jacobian)**2*(zth*zpsi+rth*rpsi)

  grad_psi_r = -r/jacobian*zth
  grad_psi_z = r/jacobian*rth
  grad_theta_r = r/jacobian*zpsi
  grad_theta_z = -r/jacobian*rpsi

  do j = 1, n
     jacobian_av(j) = sum(abs_jacobian(1:m-1,j))/(m-1) !uniform theta grid is assumed
  enddo

  allocate(a(m,n), apsi(m,n))
  allocate(b(m,n), bth(m,n))
  a = jacobian * grad_psi**2
  b = jacobian * grad_psi_dot_grad_theta
  do i = 1, m
     !call spline3ders(xgrid, a(i,:), xgrid, dynew=apsi(i,:))
     call one_dimensional_derivative(n, xgrid, a(i,:), apsi(i,:))
  enddo

  do j = 1, n
     !call spline3ders(zgrid, b(:,j), zgrid, dynew=bth(:,j))
     call one_dimensional_derivative(m, zgrid, b(:,j), bth(:,j))
  enddo

  lapx = (apsi + bth)/jacobian
  
  if(myid==0) call diagnostic()
  !if(myid==0) call plot_psi_r_z_theta_r_z_mc()
contains

  subroutine diagnostic()
    open(newunit=u,file="gradxz.txt")
    do i = 1, m
       do j = 1, n
          write(u,'(99(1pe14.5))')  r(i,j), z(i,j),  grad_psi(i,j), &
               & grad_psi_dot_grad_theta(i,j), jacobian(i,j), lapx(i,j), zgrid(i), xgrid(j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine diagnostic

!!$  subroutine plot_psi_r_z_theta_r_z_mc()
!!$    use magnetic_coordinates,only: zgrid,xgrid
!!$    integer:: u
!!$    real(p_):: minor_r_val,minor_r_prime_val
!!$
!!$    open(newunit=u, file='mc_derivatives_mc1.txt')
!!$    do i=1,m
!!$       do j=1,n
!!$          write(u,*) r(i,j),z(i,j), jacobian(i,j)
!!$
!!$          !minor_r_prime_val=minor_r_prime(xgrid(j))
!!$          !write(u,*) r(i,j),z(i,j),jacobian(i,j)/minor_r_prime_val !for equal arc-length poloidal angle, jacobian(i,j)/minor_r_prime_val is equal -R*minor_r.
!!$          !write(u,*) r(i,j),z(i,j),grad_psi_r(i,j)*minor_r_prime_val, grad_psi_z(i,j)*minor_r_prime_val,&
!!$          !write(u,*) r(i,j),z(i,j),grad_theta_r(i,j),grad_theta_z(i,j),grad_theta(i,j),jacobian(i,j)
!!$       enddo
!!$       write(u,*)
!!$    enddo
!!$    close(u)
!!$
!!$    open(newunit=u, file='psi_r_z_theta_r_z_mc_analytic.txt') !for concentric-circular magnetic field (equal-arc-length theta)
!!$    do i=1,m
!!$       do j=1,n
!!$          !minor_r_val=minor_r_radcor(xgrid(j))
!!$          write(u,*) r(i,j),z(i,j),cos(zgrid(i)),sin(zgrid(i)),&
!!$               &-sin(zgrid(i))/minor_r_val,cos(zgrid(i))/minor_r_val,-R(i,j)*minor_r_val
!!$       enddo
!!$       write(u,*)
!!$    enddo
!!$    close(u)
!!$  end subroutine plot_psi_r_z_theta_r_z_mc
end subroutine calculate_gradients_of_psi_and_theta


subroutine calculate_gradients_of_generalized_toroidal_angle() 
  use constants,only: p_, zero,one,two
  use control_parameters, only : diagnosis
  use magnetic_coordinates,only: m=>mpol, n=>nrad, r=>r_mc, z=>z_mc, xgrid, zgrid, &
       & tor_shift_mc,jacobian,dtheta,dradcor, &
       & rpsi,zpsi,rth,zth, &
       & grad_alpha_r, grad_alpha_z, grad_alpha_phi, & !as output 
       & grad_alpha, grad_psi_dot_grad_alpha, grad_alpha_dot_grad_theta, lapy, & !as output
       & grad_alpha_r_left_bdry_minus_one, grad_alpha_z_left_bdry_minus_one !as output
  use splines, only: spline3ders
  use math, only: one_dimensional_derivative
  use domain_decomposition,only: myid, dtheta2
  implicit none
  real(p_) :: tor_shift_psi(m,n), tor_shift_th(m,n)
  real(p_) :: tor_shift_psi_left_bdry_minus_one(n)
  real(p_), allocatable :: a(:,:), b(:,:), apsi(:,:), bth(:,:)
  integer :: i, j, u

  allocate(grad_alpha_r(m,n))
  allocate(grad_alpha_z(m,n))
  allocate(grad_alpha_phi(m,n))
  allocate(grad_alpha(m,n))
  allocate(grad_psi_dot_grad_alpha(m,n))
  allocate(grad_alpha_dot_grad_theta(m,n))
  allocate(grad_alpha_r_left_bdry_minus_one(n))
  allocate(grad_alpha_z_left_bdry_minus_one(n))

  call partial_derivative_of_tor_shift_in_mc(m,n,tor_shift_mc,dtheta,dradcor,tor_shift_psi,tor_shift_th, &
       & tor_shift_psi_left_bdry_minus_one) 
  !call partial_derivative_of_tor_shift_in_mc2(m,n,tor_shift_mc,dtheta,dradcor,tor_shift_psi,tor_shift_th, &
  !     & tor_shift_psi_left_bdry_minus_one) 
  grad_alpha_r = tor_shift_psi*r/jacobian*zth - tor_shift_th*r/jacobian*zpsi
  grad_alpha_z = tor_shift_th*r/jacobian*rpsi - tor_shift_psi*r/jacobian*rth
  grad_alpha_phi = one/r
  grad_alpha=sqrt(grad_alpha_phi**2+grad_alpha_r**2+grad_alpha_z**2)
  grad_psi_dot_grad_alpha=-r/jacobian*zth*grad_alpha_r+r/jacobian*rth*grad_alpha_z
  grad_alpha_dot_grad_theta=grad_alpha_r*r/jacobian*zpsi-grad_alpha_z*r/jacobian*rpsi

  !i=m-1 !wrong
  i = m - NINT(dtheta2/dtheta)
  !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
  do j = 1, n
     grad_alpha_r_left_bdry_minus_one(j)=tor_shift_psi_left_bdry_minus_one(j)*r(i,j)/jacobian(i,j)*zth(i,j)&
          & -tor_shift_th(i,j)*r(i,j)/jacobian(i,j)*zpsi(i,j)
     grad_alpha_z_left_bdry_minus_one(j)=tor_shift_th(i,j)*r(i,j)/jacobian(i,j)*rpsi(i,j)&
          & -tor_shift_psi_left_bdry_minus_one(j)*r(i,j)/jacobian(i,j)*rth(i,j)
  enddo

  allocate(a(m,n), apsi(m,n))
  allocate(b(m,n), bth(m,n))

  a = jacobian * grad_psi_dot_grad_alpha
  b = jacobian * grad_alpha_dot_grad_theta

  do i = 1, m
     !call spline3ders(xgrid, a(i,:), xgrid, dynew=apsi(i,:))
     call one_dimensional_derivative(n, xgrid, a(i,:),  apsi(i,:))
  enddo

  do j = 1, n
     !call spline3ders(zgrid, b(:,j), zgrid, dynew=bth(:,j))
     call one_dimensional_derivative(m, zgrid, b(:,j), bth(:,j))
  enddo
  allocate(lapy(m,n))
  lapy = (apsi + bth)/jacobian

  !if(myid==0 .and. (diagnosis .eqv. .true.)) call plot_alpha_r_z_mc()
  if(myid==0 .and. (diagnosis .eqv. .true.)) call verification2()
  if(myid==0) call diagnostic()
contains

  subroutine diagnostic()
    open(newunit=u,file="grady.txt")
    do i = 1, m
       do j = 1, n
          write(u,'(99(1pe14.5))')  r(i,j), z(i,j), tor_shift_mc(i,j), &
               & grad_alpha(i,j), grad_psi_dot_grad_alpha(i,j), grad_alpha_dot_grad_theta(i,j), &
               & lapy(i,j), zgrid(i), xgrid(j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine diagnostic
  
  ! subroutine plot_alpha_r_z_mc()
  !   use magnetic_coordinates,only: xgrid,zgrid
  !   use magnetic_field,only: minor_r_radcor, qfunc

  !   integer :: u
  !   real(p_) :: theta,minor_r,major_r0,major_r,q,local_q,dq_dr,factor
  !   real(p_) :: ddelta_minor_r,alpha_r,alpha_z
  !   real(p_) :: tmp,tmp2,dpsi_dr
  !   real(p_), parameter :: g0=1.32*1.91,psi_axis=0.,  psi_lcfs=  0.146018378_p_
  !   open(newunit=u,file='mc_derivatives_mc2.txt')
  !   do i=1,m
  !      do j=1,n
  !         !write(u,*) r(i,j),z(i,j),grad_alpha_r(i,j),grad_alpha_z(i,j), grad_alpha(i,j)
  !         write(u,*) r(i,j),z(i,j),tor_shift_th(i,j),tor_shift_psi(i,j)
  !         !write(u,*) r(i,j),z(i,j),tor_shift_mc(i,j)
  !      enddo
  !      write(u,*)
  !      write(u,*)
  !   enddo

  !   i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
  !   do j=1,n
  !      write(u,*) r(i,j),z(i,j), tor_shift_th(i,j),tor_shift_psi_left_bdry_minus_one(j)
  !   enddo
  !   close(u)
  !   !write(*,*) 'maxval(grad_alpha)=',maxval(grad_alpha),'minval(grad_alpha)=',minval(grad_alpha)
  !   open(newunit=u,file='alpha_r_z_mc_analytic.txt') !for concentric-circular magnetic field (equal-arc-length theta and geometric minor radius as the flux surface label)
  !   do i=1,m
  !      theta=zgrid(i)
  !      do j=1,n
  !         major_r0=1.32_p_!for DIII-D cyclone base case parameter
  !         minor_r=minor_r_radcor(xgrid(j))
  !         major_r=major_r0+minor_r*cos(theta)
  !         q=qfunc(xgrid(j))
  !         local_q=q*sqrt(major_r0**2-minor_r**2)/major_r
  !         dq_dr=0.78*1.4/0.24_p_ !for DIII-D cyclone base case parameter, a linear q profile is assumed, the same as what Ben did
  !         factor=(major_r0-minor_r)/sqrt(major_r0**2-minor_r**2)*tan(theta/two)
  !         ddelta_minor_r=two*dq_dr*atan(factor)+two*q/(1+factor**2)*tan(theta/2._p_)*&
  !              & (-major_r0)/((major_r0+minor_r)*sqrt(major_r0**2-minor_r**2))
  !         alpha_r=-ddelta_minor_r*cos(theta)+local_q*sin(theta)/minor_r
  !         alpha_z=-local_q*cos(theta)/minor_r-ddelta_minor_r*sin(theta)
  !         !write(u,*) r(i,j),z(i,j), alpha_r,alpha_z
  !         tmp=(major_r0-minor_r)/sqrt(major_r0**2-minor_r**2)
  !         tmp2=tan(theta/2)
  !         !dpsi_dr=one/minor_r_prime(xgrid(j)) 
  !         dpsi_dr=g0*minor_r/(q*sqrt(major_r0**2-minor_r**2)*(psi_lcfs-psi_axis)) !another way to calculate dpsi_dr
  !         !write(u,*) r(i,j),z(i,j), local_q, 2*q*tmp*(tan(theta/2)**2/2+0.5)/(tmp**2*tan(theta/2)**2+1)
  !         write(u,*) r(i,j),z(i,j), local_q, ddelta_minor_r/dpsi_dr !local_q is equal to ddelta_dtheta
  !         !write(u,*) r(i,j),z(i,j), (2*atan(tmp2*(major_r0 - minor_r)/sqrt(major_r0**2 - minor_r**2))*dq_dr&
  !         !     &  + 2*(tmp2*minor_r*(major_r0 - minor_r)/(major_r0**2 - minor_r**2)**(3./2) &
  !         !     & - tmp2/sqrt(major_r0**2 - minor_r**2))*q/(tmp2**2*(major_r0 - minor_r)**2&
  !         !     & /(major_r0**2 - minor_r**2) + 1))/dpsi_dr !ddelta/dradcor=ddelta/dr*(dr/dradcor), where the formula for computing ddelta/dr is obtained by using sympy

  !         !write(u,*) r(i,j),z(i,j), 2*q*atan(tmp*tan(theta/2))
  !      enddo
  !      write(u,*)
  !      write(u,*)
  !   enddo
  !   close(u)
  ! end subroutine plot_alpha_r_z_mc

  subroutine verification2() !verify grad_psi_cross_grad_alpha=B0/GSpsi_prime
    use magnetic_coordinates,only: xgrid, grad_psi_r,grad_psi_z, GSpsi_prime
    use math,only:  cross_product_in_cartesian
    use magnetic_field,only : br,bz,bphi
    real(p_):: ax,ay,az, cx,cy,cz,dx,dy,dz
    real(p_):: brval,bzval,bphival
    integer:: u,i,j
    open(newunit=u,file='cross_product_comparision.txt')

    !    i=1
    !    i=m

    i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
    do j=1,n
       ax=grad_psi_r(i,j)
       ay=0._p_
       az=grad_psi_z(i,j)

!!$       cx=grad_alpha_r(i,j)
!!$       cy=grad_alpha_phi(i,j)
!!$       cz=grad_alpha_z(i,j)

       cx=grad_alpha_r_left_bdry_minus_one(j)
       cy=grad_alpha_phi(i,j)
       cz=grad_alpha_z_left_bdry_minus_one(j)

       call cross_product_in_cartesian(ax,ay,az,cx,cy,cz,dx,dy,dz) !grad_psi_cross_grad_alpha

       brval=br(r(i,j),z(i,j))
       bzval=bz(r(i,j),z(i,j))
       bphival=bphi(r(i,j),z(i,j))
       write(u,*) dx,dy,dz, brval/GSpsi_prime,bphival/GSpsi_prime,bzval/GSpsi_prime
    enddo

    close(u)
  end subroutine verification2
end subroutine calculate_gradients_of_generalized_toroidal_angle


subroutine partial_derivative_in_mc(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)
  !calculate the partial derivative of R and Z with respect to the magnetic cooordinates (psi,theta)
  !jacob is also calculated in this subroutine
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: rpsi(m,n),rth(m,n),zpsi(m,n),zth(m,n),jacob(m,n)
  real(p_):: tmp0
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        rpsi(i,j)=(r(i,j+1)-r(i,j-1))/(two*dpsi)
        zpsi(i,j)=(z(i,j+1)-z(i,j-1))/(two*dpsi)
     enddo

     !use linear interpolation to get the value  j=n
     tmp0=(r(i,n)-r(i,n-1))/dpsi
     rpsi(i,n)=two*tmp0-rpsi(i,n-1)
     tmp0=(z(i,n)-z(i,n-1))/dpsi
     zpsi(i,n)=two*tmp0-zpsi(i,n-1)

     !use linear interpolation to get the value j=1
     tmp0=(r(i,2)-r(i,1))/dpsi
     rpsi(i,1)=two*tmp0-rpsi(i,2)

     tmp0=(z(i,2)-z(i,1))/dpsi
     zpsi(i,1)=two*tmp0-zpsi(i,2)

  enddo

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        rth(i,j)= (r(i+1,j)-r(i-1,j))/(two*dtheta)
        zth(i,j)=(z(i+1,j)-z(i-1,j))/(two*dtheta)
     enddo

     !use peroidic property of r and z to calculate the partial derivative for boundary points at theta=0 and 2pi
     rth(1,j)=(r(2,j)-r(m-1,j))/(two*dtheta)
     zth(1,j)=(z(2,j)-z(m-1,j))/(two*dtheta)
     rth(m,j)=rth(1,j)
     zth(m,j)=zth(1,j)
  enddo

  !calculate the Jacobian:
  do i=1,m
     !          do j=2,n !the jacobian at the magnetic axis is zero
     do j=1,n
        jacob(i,j)=r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,fai)
     enddo
     !jacob(i,n)=two*jacob(i,n-1)-jacob(i,n-2)
     !use linear interpolation to get the value of jacobian at the magnetic surface near the magnetic axis
  enddo

  !write(*,*) jacob(1:m,5)
end subroutine partial_derivative_in_mc

subroutine partial_derivative_in_mc2(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)
  !calculate the partial derivative of R and Z with respect to the magnetic cooordinates (psi,theta)
  !jacob is also calculated in this subroutine
  use constants,only: p_
  use magnetic_coordinates, only : radcor=>xgrid, theta=>zgrid
  use splines, only: spline3ders
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: r(m,n), z(m,n)
  real(p_), intent(in) :: dtheta, dpsi
  real(p_), intent(out) :: rpsi(m,n), rth(m,n), zpsi(m,n), zth(m,n), jacob(m,n)
  integer :: i, j

  do i = 1, m
     call spline3ders(radcor, r(i,:), radcor, dynew=rpsi(i,:))
     call spline3ders(radcor, z(i,:), radcor, dynew=zpsi(i,:))
  enddo

  do j=1,n
     call spline3ders(theta, r(:,j), theta, dynew=rth(:,j))
     call spline3ders(theta, z(:,j), theta, dynew=zth(:,j))
  enddo

  do i=1,m
     do j=1,n
        jacob(i,j) = r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,phi)
     enddo
     !jacob(i,1) = 2*jacob(i,2) - jacob(i,3)
  enddo

end subroutine partial_derivative_in_mc2


subroutine partial_derivative_of_tor_shift_in_mc(m,n,tor_shift,dtheta,dpsi,tor_shift_psi,tor_shift_th,&
     & tor_shift_psi_left_bdry_minus_one)
  !calculate the partial derivative with respect to the magnetic cooordinates (psi,theta)
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  use magnetic_coordinates,only: tor_shift_left_bdry_minus_one
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: tor_shift(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: tor_shift_psi(m,n),tor_shift_th(m,n)
  real(p_),intent(out):: tor_shift_psi_left_bdry_minus_one(n)
  real(p_):: tmp0,twopi_q
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        tor_shift_psi(i,j)=(tor_shift(i,j+1)-tor_shift(i,j-1))/(two*dpsi)
     enddo
     !use linear interpolation to get the value at j=n
     tmp0=(tor_shift(i,n)-tor_shift(i,n-1))/dpsi
     tor_shift_psi(i,n)=two*tmp0-tor_shift_psi(i,n-1)
     !use linear interpolation to get the value j=1
     tmp0=(tor_shift(i,2)-tor_shift(i,1))/dpsi
     tor_shift_psi(i,1)=two*tmp0-tor_shift_psi(i,2)
  enddo

  do j=2,n-1
     tor_shift_psi_left_bdry_minus_one(j)=&
          & (tor_shift_left_bdry_minus_one(j+1)-tor_shift_left_bdry_minus_one(j-1))/(two*dpsi)
  enddo
  !use linear interpolation to get the value at j=n
  tmp0=(tor_shift_left_bdry_minus_one(n)-tor_shift_left_bdry_minus_one(n-1))/dpsi
  tor_shift_psi_left_bdry_minus_one(n)=two*tmp0-tor_shift_psi_left_bdry_minus_one(n-1)
  !use linear interpolation to get the value j=1
  tmp0=(tor_shift_left_bdry_minus_one(2)-tor_shift_left_bdry_minus_one(1))/dpsi
  tor_shift_psi_left_bdry_minus_one(1)=two*tmp0-tor_shift_psi_left_bdry_minus_one(2)

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        tor_shift_th(i,j)= (tor_shift(i+1,j)-tor_shift(i-1,j))/(two*dtheta)
     enddo
     !for boundary points at theta cut
!!$     twopi_q=tor_shift(m,j)
!!$     tor_shift_th(1,j)=(tor_shift(2,j)+twopi_q-tor_shift(m-1,j))/(two*dtheta)
!!$     tor_shift_th(m,j)=tor_shift_th(1,j)
     !tor_shift_th(1,j)=(tor_shift(2,j)-tor_shift_left_bdry_minus_one(j))/(two*dtheta), wrong! since left_bdry_minus_one is defined on the grids for perturbations, which are different from equilibrium grids.
     tor_shift_th(1,j)=two*tor_shift_th(2,j)-tor_shift_th(3,j)   !use linear interpolation to get the value at i=1 and i=m
     tor_shift_th(m,j)=tor_shift_th(1,j) !tor_shift_th is periodic

  enddo
end subroutine partial_derivative_of_tor_shift_in_mc


subroutine partial_derivative_of_tor_shift_in_mc2(m,n,tor_shift,dtheta,dpsi,tor_shift_psi,tor_shift_th,&
     & tor_shift_psi_left_bdry_minus_one)
  !calculate the partial derivative with respect to the magnetic cooordinates (psi,theta)
  use constants,only:p_,zero,one,two,twopi,one_half
  use magnetic_coordinates,only: tor_shift_left_bdry_minus_one, xgrid, zgrid
  use splines, only: spline3ders
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: tor_shift(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: tor_shift_psi(m,n),tor_shift_th(m,n)
  real(p_),intent(out):: tor_shift_psi_left_bdry_minus_one(n)
  integer:: i,j
  real(p_) :: tmp(m,n), tmp0(n)

  do i = 1, m
     call spline3ders(xgrid, tor_shift(i,:), xgrid, tmp(i,:), tor_shift_psi(i,:), tmp(i,:))
  enddo

  do j = 1, n
     call spline3ders(zgrid, tor_shift(:,j), zgrid, tmp(:,j), tor_shift_th(:,j), tmp(:,j))
  enddo

  call spline3ders(xgrid, tor_shift_left_bdry_minus_one(:), xgrid, &
       & tmp0(:), tor_shift_psi_left_bdry_minus_one(:), tmp0(:))


end subroutine partial_derivative_of_tor_shift_in_mc2



function jacobian_func(theta,radcor) result (z) 
  use constants,only: p_
  use magnetic_coordinates,only: mpol,nrad,zgrid,xgrid,jacobian !as input
  use interpolate_module,only: linear_2d_interpolate

  implicit none
  real(p_)::radcor,theta,z
  call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,theta,radcor,z)  !uniform 1darray is assumed
end function jacobian_func
