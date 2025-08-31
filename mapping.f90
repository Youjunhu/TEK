subroutine mapping_cylindrical_to_magnetic_coordinates() !this to prepare numerical tables radcor(R,Z), theta(R,Z), and tor_shift(R,Z). Boris pusher works in cylindrical coordinates, we need these pre-mapping data to interpolate the Boris orbit into magnetic coordinates. Later I decide that tor_shift(R,Z) is not used in the rest of the program. The value of tor_shift and the various gradients of psi, theta, and alpha will be computed by interpolating the numerical table in magnetic coordinates (numerical talbes are created in "calculate_metric" routine)
  use constants, only: p_, pi
  use boundary,only: np_lcfs !,x_lcfs,z_lcfs
  use mapping_module,only:nx_mapping,nz_mapping
  use mapping_module,only:r_cyl,z_cyl,dr,dz,i0,j0 !as output
  use mapping_module,only: radcor,theta_a,theta_b,tor_shift_a,tor_shift_b !as output
  use domain_decomposition,only:myid
  use math, only : pnpoly
  use magnetic_field, only : psi_func, qfunc0
  implicit none
  logical:: within_region(nx_mapping,nz_mapping)
  integer:: i,j,inout1(nx_mapping,nz_mapping),inout2(nx_mapping,nz_mapping)
  real(p_):: theta(nx_mapping,nz_mapping),tor_shift(nx_mapping,nz_mapping)
  real(p_):: r_inner_surf(np_lcfs),z_inner_surf(np_lcfs),r_outer_surf(np_lcfs),z_outer_surf(np_lcfs)

  call choose_boundary_magnetic_surfaces_for_the_mapping(r_inner_surf,z_inner_surf,r_outer_surf,z_outer_surf)
  call create_cylindrical_grids(r_outer_surf,z_outer_surf,np_lcfs,nx_mapping,nz_mapping,r_cyl,z_cyl,dr,dz,i0,j0)

  do i=1,nx_mapping
     do j=1,nz_mapping
!!$        call PNPOLY(r_cyl(i),z_cyl(j),r_mc(:,nrad),z_mc(:,nrad),mpol,INOUT1(i,j))
!!$        call PNPOLY(r_cyl(i),z_cyl(j),r_mc(:,2),z_mc(:,2),mpol,INOUT2(i,j))
        call PNPOLY(r_cyl(i),z_cyl(j),r_outer_surf,z_outer_surf,np_lcfs,INOUT1(i,j))
        call PNPOLY(r_cyl(i),z_cyl(j),r_inner_surf,z_inner_surf,np_lcfs,INOUT2(i,j))
     enddo
  enddo

  within_region=.true.
  do i=1,nx_mapping !check whether a point is within the specifed region or not.
     do j=1,nz_mapping
        if((inout1(i,j).eq. -1) .or. (inout2(i,j).eq.1)) within_region(i,j)=.false.
     enddo
  enddo
  ! if(myid.eq.0) write(*,*) 'q_edge=',qfunc0(psi_func(r_mc(1,nrad),z_mc(1,nrad))) !value of q at boundary
!!$ if(myid.eq.0) call diagnostic1()
  ! call arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis)

  radcor=0._p_ !initialized
  theta=0._p_
  tor_shift=0._p_

  do i=1,nx_mapping
     do j=1,nz_mapping
        if(i<i0 .and. j.eq.j0) then !theta cut, special treatment is needed here
           !done at the end of this subroutine
        else
           if(within_region(i,j).eqv..true.)  call mapping(r_cyl(i),z_cyl(j),radcor(i,j),theta(i,j),tor_shift(i,j))  !calculate (radcor, theta, tor_shift) of the point (r_cyl(i),z_cyl(j)).
        endif
     enddo
  enddo

  theta_a=theta
  theta_b=theta
  tor_shift_a=tor_shift
  tor_shift_b=tor_shift

  do i=1,nx_mapping
     do j=1,nz_mapping
        if(i<i0 .and. j.eq.j0) then !theta cut, special treatment is needed here
           theta_a(i,j)=-pi
           theta_b(i,j)=pi
           call mapping2(r_cyl(i),z_cyl(j),radcor(i,j),tor_shift_a(i,j),tor_shift_b(i,j))
        endif
     enddo
  enddo

!!$  tor_shift_b=tor_shift
!!$  do i=i0,nx_mapping
!!$     j=j0 !re-calculate tor_shift_angle at low-field-side midplane (2pi*q, instead of zero) and store them in tor_shift_b array
!!$     if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$  enddo
!!$  tor_shift_b=tor_shift
!!$  do i=1,nx_mapping
!!$     do j=1,nz_mapping
!!$        if(i<i0 .and. j.eq.j0) then !the value of theta and tor_shift at the high-field-side midplane
!!$           if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$           if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$        endif
!!$     enddo
!!$  enddo
  if(myid.eq.0) call diagnostic2()
contains
  subroutine diagnostic1()
    integer:: u1,u2
    open(newunit=u1,file='mapping_in.txt')
    open(newunit=u2,file='mapping_out.txt')
    do i=1,nx_mapping
       do j=1,nz_mapping
          if(within_region(i,j) .eqv. .false.) then 
             write(u2,*) r_cyl(i),z_cyl(j) !out
          else
             write(u1,*) r_cyl(i),z_cyl(j) !in
          endif
       enddo
    enddo
    close(u1)
    close(u2)
  end subroutine diagnostic1
  subroutine diagnostic2()
    integer:: u
    open(newunit=u,file='mapping_table.txt')
    do i=1,nx_mapping
       do j=1,nz_mapping 
          !do j=j0,j0
          if(within_region(i,j).eqv..true.) then
             write(u,*) r_cyl(i),z_cyl(j),theta_a(i,j),theta_b(i,j),tor_shift_a(i,j),tor_shift_b(i,j)
             !write(u,*) r_cyl(i),z_cyl(j),psi_func(r_cyl(i),z_cyl(j))
          else
             write(u,*) r_cyl(i),z_cyl(j), 'NaN',' NaN',' NaN',' NaN'
          endif
       enddo
       write(u,*) 
    enddo
    close(u)
    !  write(*,*) 'maximum of tor_shift=',maxval(tor_shift_a),'minimum of tor_shift=',minval(tor_shift_a)
    !  write(*,*) 'maximum of tor_shift=',maxval(tor_shift_b),'minimum of tor_shift=',minval(tor_shift_b)
  end subroutine diagnostic2
end subroutine mapping_cylindrical_to_magnetic_coordinates


subroutine mapping(r,z,radcor,theta,tor_shift)
  !given (R,Z), this subroutine finds the magnetic surface that passes through the point and calculates its poloidal angle and toroidal shift
  use constants,only : p_, zero,one,two,twopi,pi
  use math, only : arc_length
  use boundary, only: x_lcfs,z_lcfs,np_lcfs
  use radial_module,only:r_axis,z_axis,psi_axis,psi_lcfs
  use magnetic_field,only:radcor_as_func_of_pfn
  use control_parameters,only: poloidal_angle_type
  use domain_decomposition,only:myid
  use calculate_toroidal_shift_module
  use contour_mod,only : contour
  use magnetic_field, only : psi_func,psi_gradient_func, b
  implicit none
  real(p_),intent(in):: r,z
  real(p_),intent(out):: radcor,theta,tor_shift

  real(p_)::psival
  real(p_):: x_contour(np_lcfs+1),z_contour(np_lcfs+1)
  real(p_)::dl(np_lcfs), sum

  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
  real(p_):: x1,x2,z1,z2

  real(p_):: slope(np_lcfs),slope2(np_lcfs)
  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
  integer:: i,end_i
  real(p_):: value1,value2,value3, rmid,zmid,normalization,tmpx,tmpz
  real(p_),parameter:: large_number=1d30

  psival=psi_func(r,z)
  radcor=radcor_as_func_of_pfn((psival-psi_axis)/(psi_lcfs-psi_axis))

  call contour(psival,x_contour,z_contour)

  call arc_length(x_contour,z_contour,np_lcfs,dl)
!!$     sum=0.
!!$     do i=1,np_lcfs-1
!!$        sum=sum+dl(i)
!!$     enddo
!!$     circumference=sum

  normalization=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
        normalization=normalization+dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        normalization=normalization+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'straight-field-line') then
     do i=2,np_lcfs
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        normalization=normalization+dl(i-1)/(rmid*psi_gradient_func(rmid,zmid))
     enddo

  elseif(poloidal_angle_type .eq. 'Boozer') then
     do i=2,np_lcfs
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        normalization=normalization+dl(i-1)*b(rmid,zmid)**2*rmid/psi_gradient_func(rmid,zmid)
     enddo

  else
     stop 'please choose poloidal angle type between equal-arc/equal-volume/straight-field-line, in mapping'
  endif

  call locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
  ! if(myid.eq.0) call diagnostic1()
  tmpx=x_contour(end_i+1) !backup the original value
  tmpz=z_contour(end_i+1)
  x_contour(end_i+1)=r !replace No. end_i+1 point by the given point (r,z)
  z_contour(end_i+1)=z
  dl(end_i)=sqrt((x_contour(end_i+1)-x_contour(end_i))**2+(z_contour(end_i+1)-z_contour(end_i))**2)

  !calculate poloidal angle 
  theta=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,end_i+1
        theta=theta+dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,end_i+1
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        theta=theta+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'straight-field-line') then
     do i=2,end_i+1
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        theta=theta+dl(i-1)/(rmid*psi_gradient_func(rmid,zmid))
     enddo
  elseif(poloidal_angle_type .eq. 'Boozer') then
     do i=2,end_i+1
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        theta=theta+dl(i-1)*b(rmid,zmid)**2*rmid/psi_gradient_func(rmid,zmid)
     enddo

  else
     stop 'please choose poloidal angle type between equal-arc and equal-volume'
  endif
  theta=theta*twopi/normalization-pi !normalized to the range [-pi:pi]

  x_contour(end_i+1)=tmpx !restore to the original value
  z_contour(end_i+1)=tmpz
  call calculate_toroidal_shift(psival,x_contour,z_contour,np_lcfs,end_i,r,z,tor_shift) !calculate toroidal shift (\int_0_theta{q_hat dtheta}) which is needed in the definition of the generalized toroidal angle
contains
  subroutine diagnostic1()
    integer:: u
    open(newunit=u,file='tmp_contour')
    do i=1,end_i
       write(u,*) x_contour(i),z_contour(i)
    enddo
    write(u,*) 
    write(u,*) 
  end subroutine diagnostic1
end subroutine mapping

subroutine mapping2(r,z,radcor,tor_shift_a,tor_shift_b)
  !given (R,Z), this subroutine find the magnetic surface that passes throught the point and calculates its toroidal angle shift
  !the same as mapping(), the difference is that this only takes care of the total_tor_shift at the theta cut, see the comments where this subroutine is called
  use constants,only:p_
  use constants,only:zero,one,two,twopi,pi
  use boundary, only: x_lcfs,z_lcfs,np_lcfs
  use radial_module,only:r_axis,z_axis,psi_axis,psi_lcfs
  use calculate_toroidal_shift_module
  use contour_mod,only : contour
  use magnetic_field, only : psi_func, radcor_as_func_of_pfn
  implicit none
  real(p_),intent(in):: r,z
  real(p_),intent(out):: radcor,tor_shift_a,tor_shift_b
  real(p_)::psival
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  real(p_)::dl(np_lcfs-1)
  integer:: i,end_i

  psival=psi_func(r,z)
  radcor=radcor_as_func_of_pfn((psival-psi_axis)/(psi_lcfs-psi_axis))
  call contour(psival,x_contour,z_contour)

  !  call calculate_toroidal_shift_total(psival,x_contour,z_contour,np_lcfs,dl,tor_shift_a) 
  call calculate_toroidal_shift_at_theta_cut(psival,x_contour,z_contour,np_lcfs,tor_shift_a,tor_shift_b) !calculate toroidal shift which is needed in the definition of the generalized toroidal angle
end subroutine mapping2


subroutine locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
  use constants,only:p_
  use constants,only:twopi,pi
  use radial_module,only:r_axis,z_axis
  implicit none
  real(p_),intent(in):: r,z
  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  integer,intent(out):: end_i !!poloidal index of point (r,z) is between end_i and end_i+1
  real(p_):: xn,zn,xn0,zn0,angle0,angle1,angle2
  integer:: i

  xn0=r-r_axis
  zn0=z-z_axis
  angle0=acos(xn0/sqrt(xn0*xn0+zn0*zn0))
  if(zn0.le.0.0) angle0=-angle0 !theta in [-pi:pi]

  do i=1,np_lcfs-1
     xn=x_lcfs(i)-r_axis
     zn=z_lcfs(i)-z_axis
     angle1=acos(xn/sqrt(xn**2+zn**2)) 
     if(zn<0.0) angle1=-angle1
     if(i.eq.1) angle1=-pi

     xn=x_lcfs(i+1)-r_axis
     zn=z_lcfs(i+1)-z_axis
     angle2=acos(xn/sqrt(xn**2+zn**2))
     if(zn<0.0) angle2=-angle2 !theta in [-pi:pi]
     if(i+1.eq.np_lcfs) angle2=pi

     if((angle0-angle1)*(angle0-angle2).le.0._p_) exit

  enddo


  end_i=i

end subroutine locate_poloidal_index

subroutine choose_boundary_magnetic_surfaces_for_the_mapping(r_inner_surf,z_inner_surf,r_outer_surf,z_outer_surf)
  use constants,only:p_
  use constants,only: two
  use magnetic_coordinates,only: pfn_inner,pfn_bdry !boundary for the region in which magnetic_coordinates are constructed
  use radial_module,only: r_axis,z_axis,psi_axis,psi_lcfs
  use boundary,only: x_lcfs,z_lcfs,np_lcfs
  use contour_mod,only : contour
  implicit none
  real(p_),intent(out):: r_inner_surf(np_lcfs),z_inner_surf(np_lcfs),r_outer_surf(np_lcfs),z_outer_surf(np_lcfs)
  real(p_):: psi_val,psi_bdry,psi_inner

  psi_inner=psi_axis+pfn_inner*(psi_lcfs-psi_axis) 
  psi_val=(psi_axis+psi_inner)/two !the inner boundary for doing the mapping
  call contour(psi_val,r_inner_surf,z_inner_surf)

  psi_bdry=psi_axis+pfn_bdry*(psi_lcfs-psi_axis)
  psi_val=(psi_lcfs+psi_bdry)/two !the outer boundary for doing the mapping, which is chosen to be between lcfs and the outer boundary of the region in which magnetic_coordinates are availabe
  call contour(psi_val,r_outer_surf,z_outer_surf)

end subroutine choose_boundary_magnetic_surfaces_for_the_mapping

subroutine create_cylindrical_grids(r_outer_surf,z_outer_surf,np_lcfs,nx,nz,r,z,dr,dz,i0,j0)
  !create rectangular box (with cylindrical grids in it) in the poloidal plane with the boundary flux surface within the box and (r_axis,z_axis) is exactly on a grid point
  use constants,only:p_
  !  use magnetic_coordinates,only:r_mc, z_mc,mpol,nrad
  use radial_module,only: r_axis,z_axis
  implicit none
  integer,intent(in):: np_lcfs,nx,nz
  real(p_),intent(out):: r(nx),z(nz),dr,dz
  integer,intent(out):: i0,j0 !index of the point at magnetic axis
  real(p_):: r_min, r_max,z_min,z_max,r_width,z_width
  !  real(p_):: r_mag_surf1(mpol), z_mag_surf1(mpol)
  real(p_):: r_outer_surf(np_lcfs),z_outer_surf(np_lcfs)
  integer:: i,j,nxp,nzp
  real(p_):: rp(nx-1),zp(nz-1)


  nxp=nx-1 !using a reduced number, so that I can append the array with an additional element
  nzp=nz-1 !using a reduced number, so that I can append the array with an additional element

!!$  do i=1,mpol !select the boundary magnetic surface
!!$     r_mag_surf1(i)=r_mc(i,nrad)
!!$     z_mag_surf1(i)=z_mc(i,nrad)
!!$  enddo

  r_min=minval(r_outer_surf)
  r_max=maxval(r_outer_surf)
  r_width=r_max-r_min

  z_min=minval(z_outer_surf)
  z_max=maxval(z_outer_surf)
  z_width=z_max-z_min

  dr=r_width/(nxp-1)
  dz=z_width/(nzp-1)

  !i0=nx/2 !i index at magnetic axis
  i0=floor((r_axis-r_min)/dr)+1 !i index near the magnetic axis
  j0=floor((z_axis-z_min)/dz)+1 !j index near the magnetic axis
  rp(i0)=r_axis !set (i0,j0) to be on the magnetic axis, since there is a floor operation when getting (i0,j0), this shifts the box to the low-field-side (lfs) and upward
  zp(j0)=z_axis 

  !I do the above steps because I want the magnetic axis to lie exactly on a grid so that I can get grids that exactly represent the midplane.
  do i=i0-1,1,-1 !the grids at the high-field-side (hfs) of the magnetic axis
     rp(i)=rp(i+1)-dr
  enddo
  do i=i0+1,nxp,1 !the girds at the low-field-side (lfs) of the magnetic axis
     rp(i)=rp(i-1)+dr
  enddo

  do j=j0-1,1,-1 !the grids below the midplane
     zp(j)=zp(j+1)-dz
  enddo
  do j=j0+1,nzp,1 !the grids above the midplane
     zp(j)=zp(j-1)+dz
  enddo

  !add an additional element. This is needed because setting the magnetic axis to be on a grid usually shifts the box to the lfs, which will make some points on the hfs not included in the box
  r(1)=rp(1)-dr
  do i=1,nxp
     r(i+1)=rp(i)
  enddo

  z(1)=zp(1)-dz !similar reason as mentioned above
  do j=1,nzp
     z(j+1)=zp(j)
  enddo

  i0=i0+1 !the i index of the magnetic axis at the new array r
  j0=j0+1 !the j index of the magnetic axis at the new array z
  block
    use domain_decomposition, only : myid
    integer :: u
    if(myid==0) then
       open(newunit=u,file='rzgrid.txt')
       do i =1,nx
          do j =1,nz
             write(u,*) r(i), z(j)
          enddo
       enddo
       close(u)
    end if
  endblock
end subroutine create_cylindrical_grids


module map_to_mc
contains
  subroutine interpolate_from_cylindrical_to_magnetic_coordinates(r0,z0,theta0, tor_shift0)
    !calculate the magnetic coordinats (theta0,tor_shift0) of (R0,Z0) by using interpolation,
    !radcor0 is not calculated in this subroutine because I can use another reliable way to calculate radcor0 from (r0,z0), i.e., radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))
    !to use this iterpolation, first to make sure that (r0,z0) is within the specfied region
    use constants,only:p_,two,twopi
    use mapping_module,only: r_cyl,z_cyl,dr,dz,radcor,theta_a,theta_b,tor_shift_a,i0,j0,tor_shift_b
    use interpolate_module,only: linear_2d_interpolate_kernel
    implicit none
    real(p_),intent(in):: r0,z0
    real(p_),intent(out):: theta0,tor_shift0
    real(p_):: radcor_tmp(2,2),theta_tmp(2,2),tor_shift_tmp(2,2),qval
    integer::i,j,ii,jj
    real(p_):: radcor_as_func_of_pfn,pfn_func

    i=floor((r0-r_cyl(1))/dr)+1
    j=floor((z0-z_cyl(1))/dz)+1

    !  2D interpolations to get radcor0, theta0, and tor_shift0
!!$  do ii=1,2
!!$     do jj=1,2
!!$        radcor_tmp(ii,jj)=radcor(i+ii-1,j+jj-1)
!!$     enddo
!!$  enddo
!!$
!!$  call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),radcor_tmp,r0,z0,radcor0)

    !  radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))


    do ii=1,2
       do jj=1,2
          theta_tmp(ii,jj)=theta_a(i+ii-1,j+jj-1)
       enddo
    enddo

!!$ if(i>i0 .and. j+1.eq.j0) then !handle the boundary at the low-field-side midplane
!!$theta_tmp(1,2)=twopi 
!!$theta_tmp(2,2)=twopi 
!!$endif
    if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side (theta cut)
       theta_tmp(1,1)=theta_b(i,j0)
       theta_tmp(2,1)=theta_b(i+1,j0)
    endif
    call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),theta_tmp,r0,z0,theta0)

    do ii=1,2
       do jj=1,2
          tor_shift_tmp(ii,jj)=tor_shift_a(i+ii-1,j+jj-1)
       enddo
    enddo

!!$if(i>i0 .and. j+1.eq.j0) then !handle the boundary at the low-field-side midplne
!!$!qval=2.3074808101594884*1.0004
!!$!tor_shift_tmp(1,2)=two*tor_shift(i,j)-tor_shift(i,j-1) 
!!$!tor_shift_tmp(1,2)=qval*twopi
!!$!tor_shift_tmp(2,2)=qval*twopi
!!$tor_shift_tmp(1,2)=tor_shift_b(i,j0)
!!$tor_shift_tmp(2,2)=tor_shift_b(i+1,j0)
!!$endif

    if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side midplane (theta cut)
       tor_shift_tmp(1,1)=tor_shift_b(i,j0)
       tor_shift_tmp(2,1)=tor_shift_b(i+1,j0)
    endif

    call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),tor_shift_tmp,r0,z0,tor_shift0)

    !alpha0=phi0-tor_shift0

  end subroutine interpolate_from_cylindrical_to_magnetic_coordinates

  pure  subroutine interpolate_from_cylindrical_to_magnetic_coordinates1(r0,z0,theta0)
    use constants,only:p_,two,twopi, pi
    use mapping_module,only: r_cyl,z_cyl,dr,dz,radcor,theta_a,theta_b,i0,j0
    use interpolate_module

    implicit none
    real(p_),intent(in):: r0,z0
    real(p_),intent(out):: theta0
    real(p_):: radcor_tmp(2,2),theta_tmp(2,2),tor_shift_tmp(2,2),qval
    integer::i,j,ii,jj
    real(p_):: radcor_as_func_of_pfn,pfn_func
    !locate
    i=floor((r0-r_cyl(1))/dr)+1
    j=floor((z0-z_cyl(1))/dz)+1

    !  2D interpolations to get radcor0, theta0, and tor_shift0
!!$  do ii=1,2
!!$     do jj=1,2
!!$        radcor_tmp(ii,jj)=radcor(i+ii-1,j+jj-1)
!!$     enddo
!!$  enddo
!!$
!!$  call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),radcor_tmp,r0,z0,radcor0)

    !  radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))

    do ii=1,2
       do jj=1,2
          theta_tmp(ii,jj)=theta_a(i+ii-1,j+jj-1)
       enddo
    enddo

    if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side (theta cut)
       theta_tmp(1,1)=theta_b(i,j0)
       theta_tmp(2,1)=theta_b(i+1,j0)
    endif
    call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),theta_tmp,r0,z0,theta0)


!!$  do ii=1,2
!!$     do jj=1,2
!!$        tor_shift_tmp(ii,jj)=tor_shift_a(i+ii-1,j+jj-1)
!!$     enddo
!!$  enddo
!!$
!!$ if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side midplane (theta cut)
!!$     tor_shift_tmp(1,1)=tor_shift_b(i,j0)
!!$     tor_shift_tmp(2,1)=tor_shift_b(i+1,j0)
!!$  endif
!!$
!!$  call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),tor_shift_tmp,r0,z0,tor_shift0)

  end subroutine interpolate_from_cylindrical_to_magnetic_coordinates1
end module map_to_mc
