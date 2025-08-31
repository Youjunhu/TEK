module contour_mod
contains
  subroutine contour(psival, x_contour, z_contour)
    ! Given a value of the poloidal flux, psival,
    ! this subroutine find the magnetic surface corresponding to psival
    use constants, only: p_
    use boundary, only: x_lcfs, z_lcfs, np_lcfs
    use radial_module, only: r_axis, z_axis
    implicit none
    real(p_), intent(in) :: psival
    real(p_), intent(out) :: x_contour(np_lcfs), z_contour(np_lcfs)
    real(p_), parameter:: xacc = 1.0d-6 !tolerance used in bi-section root-finder
    real(p_), parameter:: huge = 1d30
    real(p_) :: x1, x2, z1, z2, slope(np_lcfs), slope2(np_lcfs)
    integer :: i

    do i=1,np_lcfs
       if(x_lcfs(i)-r_axis .ne. 0._p_) then 
          slope(i)= (z_lcfs(i)-z_axis)/(x_lcfs(i)-r_axis) !the slope for function Z=Z(X)
       else
          slope(i) = huge !I use compiler option that catches all erroneous arithmetic operation, I need to avoid dividing by zero
       endif
       if(z_lcfs(i)-z_axis .ne. 0._p_) then
          slope2(i)=(x_lcfs(i)-r_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
       else
          slope2(i) = huge
       endif
    enddo

    do i=1,np_lcfs -1 
       if(abs(slope(i)).le.1.0_p_) then !use Z=Z(X) function, the reason that I switch between using function X=X(Z) and Z=Z(X) is to aviod large slope.
          x1=r_axis
          x2=x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
          x_contour(i)=rtbis(one_dim_psi_func,x1,x2,xacc,r_axis,z_axis,slope(i),psival)
          z_contour(i)=zfunc(r_axis,z_axis,slope(i),x_contour(i))
       else !switch to using X=X(Z) function
          z1=z_axis
          z2=z_lcfs(i)
          z_contour(i)=rtbis(one_dim_psi_func2,z1,z2,xacc,r_axis,z_axis,slope2(i),psival)
          x_contour(i)=xfunc(r_axis,z_axis,slope2(i),z_contour(i)) 
       endif
    enddo

    x_contour(np_lcfs) = x_contour(1) !i=1 and i=np_lcfs are identical
    z_contour(np_lcfs) = z_contour(1) 

  end subroutine contour


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
  use magnetic_coordinates, only: mpol, mtor, nsegment, nrad, pfn, GSpsi_array, &
       & r_mc,z_mc, & !as output
       & zgrid, dtheta, & !as output
       & minor_r_array,minor_r_prime_array, i_theta_zero, dl_mc, & !output
       & toroidal_range, dtor, ygrid, vol !output
  use domain_decomposition, only: myid
  use contour_mod, only: contour
  use splines, only: spline3ders
  use math, only:  one_dimensional_derivative, arc_length
  implicit none
  real(p_), allocatable :: r_mag_surf0(:,:), z_mag_surf0(:,:)
  integer :: i, j

  call choose_radial_grids()
  
  allocate(r_mag_surf0(np_lcfs,nrad))
  allocate(z_mag_surf0(np_lcfs,nrad))
  
  do j=1,nrad
     call contour(GSpsi_array(j),r_mag_surf0(:,j),z_mag_surf0(:,j))
  enddo

  if((myid.eq.0) .and. (diagnosis.eqv..true.)) call diagnostic1()
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
  zgrid = [ (-pi+dtheta*(i-1), i=1,mpol) ]
  i_theta_zero = (mpol+1)/2 !poloidal index corresponding to theta=0 (mpol is assumed odd)

  allocate(r_mc(mpol,nrad))
  allocate(z_mc(mpol,nrad))

  do j=1,nrad
     call construct_poloidal_coordinate(r_mag_surf0(:,j),z_mag_surf0(:,j),np_lcfs, &
          & mpol,zgrid, r_mc(:,j), z_mc(:,j))
  enddo

  allocate(dl_mc(mpol-1,nrad))
  do j=1,nrad
     call arc_length(r_mc(:,j), z_mc(:,j), mpol, dl_mc(:,j))
  enddo

  call calculate_metric() 

  toroidal_range = twopi/nsegment
  dtor = toroidal_range/mtor
  allocate(ygrid(mtor+1)) 
  ygrid = [ (zero+dtor*(i-1), i=1,mtor+1) ]

  call plasma_volume_of_computational_region(vol)
  if ((myid==0) .and. (diagnosis .eqv. .true.)) call diagnostic2()
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
  use radial_module, only: psi_axis, psi_lcfs, j_fixed, radcor_fixed
  use magnetic_field, only: radcor_as_func_of_pfn
  use magnetic_coordinates, only: nrad, pfn_inner, pfn_bdry, &
       & nrad, pfn,GSpsi_array,GSpsi_prime, xgrid,dradcor, & !output
       & xlow, xupp, & !as output
       & radial_width !as output
  implicit none
  real(p_) :: dpfn
  integer :: j

  allocate(pfn(nrad))  
  allocate(GSpsi_array(nrad))
  allocate(xgrid(nrad))
  dpfn=(pfn_bdry-pfn_inner)/real(nrad-1)
  do j = 1, nrad !select some flux surfaces (labeld by GSpsi_array)
     pfn(j) = pfn_inner +dpfn*(j-1)
     GSpsi_array(j)=psi_axis+pfn(j)*(psi_lcfs-psi_axis) !GSpsi=Aphi*R, this array is used in finding magnetic surfaces
     xgrid(j)=radcor_as_func_of_pfn(pfn(j))
  enddo

  xlow = xgrid(1)
  xupp = xgrid(nrad)
  
  radial_width = xgrid(nrad) - xgrid(1)
  dradcor = xgrid(2) - xgrid(1) !radial grid interval
  GSpsi_prime = psi_lcfs - psi_axis !dGSpsi/dx, x is the normalized poloidal magnetic flux

  j_fixed = nrad/2
  radcor_fixed = xgrid(j_fixed) !the radcor of the center of computational region, used in flux tube model

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

  theta_old = theta_old*twopi/theta_old(mpol_old) - pi !normalized to the range [-pi:pi]

  !interpolate R  to theta_new gridpoints
  call spline(theta_old, r_old, mpol_old, 2.d30, 2.d30, y2) !prepare the second order derivative needed in the cubic spline interpolation
  do i=2,mpol-1
     call splint(theta_old,r_old,y2,mpol_old,theta_new(i), r_new(i))
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

  if(myid.eq.0) write(*,*) 'volume of the toroidal wedge=', vol
  if(myid.eq.0) write(*,*) 'volume of full torus of the computational domain=', vol*twopi/toroidal_range

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


function abs_jacobian_func(theta,radcor) result (z) !!used in generating non-uniformly distributed random numbers that satisfied a probability density function proportional to the |jacobian|
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
  integer, intent(in) :: m
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
