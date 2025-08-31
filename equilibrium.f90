subroutine read_and_process_equilibrium()
  use constants,only: p_, zero,one,two,three,four,five,twopi 
  use poloidal_flux_2d,only:xarray,zarray,nx,nz,psi,psi_gradient,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx
  use radial_module,only:psi_axis,psi_lcfs,npsi,psi_1d,fpsi,qpsi,pfn_npsi,tfn_npsi,baxis,r_axis,z_axis 
  use domain_decomposition,only:myid,numprocs
  implicit none
  integer :: i, u

  call read_gfile()
  
  allocate(psi_x(nx,nz))
  allocate(psi_z(nx,nz))
  allocate(psi_xx(nx,nz))
  allocate(psi_zz(nx,nz))
  allocate(psi_xz(nx,nz))
  allocate(psi_zx(nx,nz))
  allocate(psi_gradient(nx,nz))
  
  !call calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  call calculate_poloidal_flux_partial_derivatives2(nx,nz,xarray,zarray,psi,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)

  npsi=nx !nx in g-file is also used to define the number of radial grid points.
  allocate(psi_1d(npsi))
  do i=1,npsi
     psi_1d(i)=psi_axis+(psi_lcfs-psi_axis)/(npsi-1)*(i-1) 
  enddo

  allocate(tfn_npsi(npsi))
  allocate(pfn_npsi(npsi))

  call calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
  pfn_npsi=(psi_1d-psi_1d(1))/(psi_1d(npsi)-psi_1d(1))

  if(myid==0) then
     open(newunit=u,file='q1.txt')
     do i=1,npsi
        write(u,*) psi_1d(i), pfn_npsi(i), tfn_npsi(i), qpsi(i)
     enddo
     close(u)
  endif
  
end subroutine read_and_process_equilibrium


subroutine read_gfile()
  use constants, only : p_, two
  use poloidal_flux_2d, only: xarray, zarray, nx, nz, psi !as output
  use radial_module,only: psi_axis,psi_lcfs ,fpsi,ffprime,fprime,qpsi,baxis,r_axis,z_axis, pressure, pprime, sign_bphi
  use radial_module, only : 
  use boundary,only: nlim,rlim,zlim,np_lcfs, x_lcfs,z_lcfs !as output
  use math, only : arrange_lcfs
  use domain_decomposition, only : myid
  implicit none
  character(100) :: gfile_name
  logical :: reverse_tf,reverse_ip
  character(len=100) :: format1, format2, format3
  character(len=8) :: ntitle(5), vid
  integer:: neq,ipestg
  real(p_)::xdim, zdim, rleft, zmid
  real(p_)::  btorus, current, major_R
  real(p_):: dumaraya4(4),dumarayb5(5)
  integer:: i,j, u

  namelist /magnetic_configuration/gfile_name, reverse_tf, reverse_ip
  open(newunit=u, file='input.nmlt')
  read(u, magnetic_configuration)
  close(u)
  if(myid==0)  write(*,magnetic_configuration)

  format1='(5e16.9)'
  format2='(6a8, 3i4)'
  format3='(2i5)'
  !open and read in eqdsk file (refer to G EQDSK.pdf (or weqdsku.f in onetwo) for the gfile format)
  open(newunit=neq,file=gfile_name,status='old')
  ipestg = 4
  read (neq, format2) (ntitle(i), i=1,5), vid, ipestg, nx, nz
  allocate(psi(nx,nz))
  allocate(qpsi(nx))
  allocate(fpsi(nx))
  allocate(ffprime(nx))
  allocate(pressure(nx))
  allocate(pprime(nx))
  read (neq, format1) xdim, zdim, major_R, rleft, zmid
  read (neq, format1) r_axis, z_axis, psi_axis, psi_lcfs, btorus
  read (neq, format1) current, dumaraya4
  read (neq, format1) dumarayb5
  read (neq ,format1) (fpsi(i), i=1,nx)
  read (neq ,format1) (pressure(i), i=1,nx)
  read (neq ,format1) (ffprime(i), i=1,nx)
  read (neq ,format1) (pprime(i), i=1,nx)
  read (neq ,format1) ((psi(i,j), i=1,nx), j=1,nz)
  read (neq ,format1) (qpsi(i), i=1,nx)
  read (neq ,format3) np_lcfs, nlim
  allocate(x_lcfs(np_lcfs))
  allocate(z_lcfs(np_lcfs))
  !    nlim=56
  allocate(rlim(nlim))
  allocate(zlim(nlim))
  read (neq ,format1) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
  read (neq ,format1) (rlim(i), zlim(i), i=1,nlim)
  close(neq)


!!$    open(newunit=u,file='hl/first_wall.txt')
!!$    do i=1,nlim
!!$       read(u,*) rlim(i), zlim(i)
!!$    enddo
!!$    close(u)

  !to verify that I have read the eqdsk file correctly, I write the data read to a new file called 'tmp.gfile'.
  !After the program finished, I compare the file 'tmp.gfile' with the original file using diff command
  !the output of diff command indicates that the two files are idential, which shows I have read the eqdsk file correctly
  if(myid==0) then
!!$         open(newunit=neq,file='tmp.gfile')
!!$         write (neq, format2) (ntitle(i), i=1,5), vid, ipestg, nx, nz
!!$         write (neq, format1) xdim, zdim, major_R, rleft, zmid
!!$         write (neq, format1) r_axis, z_axis, psi_axis, psi_lcfs, btorus
!!$         write (neq, format1) current, dumaraya4
!!$         write (neq, format1) dumarayb5
!!$         write (neq ,format1) (fpsi(i), i=1,nx)
!!$         write (neq ,format1) (pressure(i), i=1,nx)
!!$         write (neq ,format1) (ffprime(i), i=1,nx)
!!$         write (neq ,format1) (pprime(i), i=1,nx)
!!$         write (neq ,format1) ((psi(i,j), i=1,nx), j=1,nz)
!!$         write (neq ,format1) (qpsi(i), i=1,nx)
!!$         write (neq ,format3) np_lcfs, nlim
!!$         write (neq ,format1) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
!!$         write (neq ,format1) (rlim(i), zlim(i), i=1,nlim)
!!$         close(neq)
  endif


  !Somtimes, I alter some quantities (e.g. reverse the toroidal magnetic field, or increase the pressure by a constant),in this case, the out gfile is different from the original one
!!$  do i=1,nx
!!$     pressure(i)=pressure(i)+0.5*pressure(1) !increase the presure
!!$  enddo

  if(reverse_tf.eqv..true.)  then
     fpsi=-fpsi !revert the toroidal magnetic field
  endif
  baxis=fpsi(1)/r_axis

  if(reverse_ip.eqv..true.) then
     psi=-psi !revert direction of the torodial current
     psi_axis=-psi_axis
     psi_lcfs=-psi_lcfs
     ffprime=-ffprime
     pprime=-pprime
  endif

  if(fpsi(1)>0) then 
     sign_bphi=1
     if(myid==0) write(*,*) 'bphi>0'
  else
     sign_bphi=-1
     if(myid==0) write(*,*) 'bphi<0'
  endif

  if(myid == 0) then
     open(newunit=u,file='lcfs.txt')
     do i=1,np_lcfs
        write(u,*) x_lcfs(i), z_lcfs(i)
     enddo
     close(u)
  endif

  allocate(xarray(nx))
  allocate(zarray(nz))

  do i=1,nx !construct the X array
     xarray(i)=rleft+xdim/(nx-1)*(i-1)
  enddo

  do j=1,nz !construct the Z array
     zarray(j)=(zmid-zdim/two)+zdim/(nz-1)*(j-1)
  enddo

  allocate(fprime(nx))
  fprime=ffprime/fpsi
  call arrange_lcfs(x_lcfs, z_lcfs, np_lcfs, r_axis, z_axis)     


  if(myid==0) then
     write(*,*) 'nR, nZ =', nx, nz
     write(*,*) 'rleft=',rleft, 'rright=', rleft+xdim,'zlow=',zmid-zdim/2._p_,'zupp=',zmid+zdim/2._p_
     write(*,*) 'magnetic location (meter) (r,z)=', r_axis, z_axis
     write(*,*) 'baxis (Tesla)=', baxis
     write(*,*) 'rcenter=',major_R, 'vacuum magnetic field at rcenter=',btorus
     !write(*,*)  'total current in all TF coils (MegaAmpere)=', btorus*major_R/(2d-7)/10**6 !Ampere's circuital law
     write(*,*) 'psi_axis=',psi_axis,'psi_lcfs=',psi_lcfs
     write(*,*) 'np_lcfs=',np_lcfs
     !write(*,*) 'x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)=', x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)
     write(*,*) 'Cyclotron angular frequency of Deuterium ion at magnetic axis (10^6 rad/s)=', &
          & fpsi(1)/r_axis*1.6022d-19/3.3452d-27/1.d6

     if(psi_lcfs>psi_axis) then
        write(*,*) 'Iphi<0'
     else
        write(*,*) 'Iphi>0'
     endif
     block !find out the direction of the troidal current
       use math, only : laplace_cylindrical2d
       real(p_), allocatable:: jphi(:,:)
       real(p_) :: dx, dz, s

       allocate(jphi(nx,nz))
       call laplace_cylindrical2d(psi(1:nx,1:nz), xarray, zarray, nx, nz, jphi)
       if (sum(jphi(nx/5:nx/2,nz/5:nz/2))>0) then
          write(*,*) 'Jphi>0'
       else
          write(*,*) 'Jphi<0'
       endif
       dx = xarray(2) -xarray(1)
       dz = zarray(2) -zarray(1)
       if(myid==0) then
          open(newunit=u,file='psi_and_jphi.txt')
          s = 0
          do i=1,nx
             do j =1,nz
                write(u,'(4e16.5)') xarray(i), zarray(j), psi(i,j), jphi(i,j)
                s = s + jphi(i,j)*dx*dz
             enddo
             write(u,*)
          enddo
          write(*,*) 'total plasma current calculated from the poloidal magnetic flux (kA)=', s/1000.
          write(*,*) "total plasma current given in g-file (kA)= ", current/1000.
          close(u)
       endif

     end block
     open(newunit=u,file='limiter.txt')
     do i=1,nlim
        write(u,*) rlim(i), zlim(i)
     enddo
     close(u)
     write(*,*) 'wall r_min=',minval(rlim), 'wall r_max=',maxval(rlim), 'wall z_min=',minval(zlim), 'wall z_max=',maxval(zlim)
  endif

end subroutine read_gfile


subroutine calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,&
     & psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  !subroutine calculate_poloidal_flux_gradient(nx,nz,xarray,zarray,psi,psi_gradient)
  use constants,only:p_
  use constants,only:one,two,twopi
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: psi(nx,nz)
  real(p_),intent(in):: xarray(nx),zarray(nz)
  real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz),psi_xx(nx,nz),psi_zz(nx,nz),psi_xz(nx,nz),psi_zx(nx,nz)
  real(p_),intent(out):: psi_gradient(nx,nz)

  integer:: i,j,i1,i2,j1,j2

  !first-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_x(i,j)=(psi(i2,j)-psi(i1,j))/(xarray(i2)-xarray(i1))
        psi_z(i,j)=(psi(i,j2)-psi(i,j1))/(zarray(j2)-zarray(j1))
     enddo
  enddo

  !second-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_xx(i,j)=(psi_x(i2,j)-psi_x(i1,j))/(xarray(i2)-xarray(i1))
        psi_zz(i,j)=(psi_z(i,j2)-psi_z(i,j1))/(zarray(j2)-zarray(j1))
        psi_xz(i,j)=(psi_x(i,j2)-psi_x(i,j1))/(zarray(j2)-zarray(j1))
        psi_zx(i,j)=(psi_z(i2,j)-psi_z(i1,j))/(xarray(i2)-xarray(i1))
     enddo
  enddo


  psi_gradient=sqrt(psi_x**2+psi_z**2)


end subroutine calculate_poloidal_flux_partial_derivatives


subroutine calculate_poloidal_flux_partial_derivatives2(nx,nz,xarray,zarray,psi,psi_x,psi_z,&
     & psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  use constants,only:p_, one,two,twopi
  use splines, only: spline3ders
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: xarray(nx),zarray(nz), psi(nx,nz)
  real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz),psi_xx(nx,nz),psi_zz(nx,nz),psi_xz(nx,nz),psi_zx(nx,nz)
  real(p_),intent(out):: psi_gradient(nx,nz)
  real(p_) :: tmp(nx,nz)
  integer :: i, j
  
  do i = 1, nx
     call spline3ders(zarray, psi(i,:), zarray, dynew=psi_z(i,:), d2ynew=psi_zz(i,:))
  enddo

  do j = 1, nz
     call spline3ders(xarray, psi(:,j), xarray, dynew=psi_x(:,j), d2ynew=psi_xx(:,j))
  enddo

  do i = 1, nx
     call spline3ders(zarray, psi_x(i,:), zarray, dynew=psi_xz(i,:))
  enddo

!!$  psi_zx = psi_xz
  do j = 1, nz
     call spline3ders(xarray, psi_z(:,j), xarray, dynew=psi_zx(:,j))
  enddo

  
 psi_gradient=sqrt(psi_x**2+psi_z**2)

end subroutine calculate_poloidal_flux_partial_derivatives2

subroutine calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
  !calculate the toroidal magnetic flux
  use constants,only:p_
  use constants,only: two,twopi
  use domain_decomposition,only:myid
  implicit none
  integer,intent(in):: npsi
  real(p_),intent(in):: psi_1d(npsi),qpsi(npsi)
  real(p_),intent(out):: tfn_npsi(npsi)
  real(p_):: dpsi,tf_npsi(npsi),total_tf
  integer:: j

  dpsi=(psi_1d(npsi)-psi_1d(1))/(npsi-1)
  tf_npsi(1)=0._p_
  do j=2,npsi 
     tf_npsi(j)=tf_npsi(j-1)+(qpsi(j)+qpsi(j-1))/two*twopi*dpsi  !using the formula dtf=q*dpf=q*d(pf_gs)*twopi
  enddo

  if(myid.eq.0) then
     open(213,file='pf_tf.txt')
     do j=1,npsi
        write(213,*) psi_1d(j),tf_npsi(j)
     enddo
     close(213)
  endif
  total_tf=tf_npsi(npsi)
  tfn_npsi= tf_npsi/total_tf !normalized toroidal magnetic flux
end subroutine calculate_tfn


subroutine draw_rect_region(nx,nz,r_1d,z_1d)
  use constants,only:p_
  implicit none
  integer,intent(in):: nx,nz
  real(p_),intent(in):: r_1d(nx),z_1d(nz)
  integer:: i,j

  open(12,file='rectangular')
  do j=1,nz
     write(12,*) r_1d(1),z_1d(j)
  enddo
  do i=1,nx
     write(12,*) r_1d(i),z_1d(nz)
  enddo
  do j=1,nz
     write(12,*) r_1d(nx),z_1d(nz-j+1)
  enddo

  do i=1,nx
     write(12,*) r_1d(nx-i+1),z_1d(1)
  enddo
  close(12)
end subroutine draw_rect_region


module magnetic_field
contains
  !all independent and dependent variables are in SI units
  pure real(p_) function psi_func(x,z) result(f)!SI units, poloidal magnetic flux function
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi,x,z,f)  
  end function psi_func

  pure real(p_) function pfn_func(x,z) result(f)!normalized poloidal magnetic flux
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi
    use radial_module,only:psi_axis,psi_lcfs
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z

    call linear_2d_interpolate(nx,nz,xarray,zarray,psi,x,z,f)  
    f=(f-psi_axis)/(psi_lcfs-psi_axis)
  end function pfn_func

  pure real(p_) function psi_r_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_x
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_x,x,z,psi_r_func)  
  end function psi_r_func


  pure real(p_) function psi_z_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_z !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in):: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_z,x,z,psi_z_func)  
  end function psi_z_func


  pure real(p_) function psi_rr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xx !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_xx,x,z,psi_rr_func)  
  end function psi_rr_func

  pure real(p_) function psi_zz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zz !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_zz,x,z,psi_zz_func)  
  end function psi_zz_func


  pure real(p_) function psi_rz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xz 
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_xz,x,z,psi_rz_func)  
  end function psi_rz_func


  pure real(p_) function psi_zr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zx 
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_zx,x,z,psi_zr_func)  
  end function psi_zr_func

  pure real(p_) function psi_gradient_func(xval,zval) result (f)
    use constants,only:p_
    use poloidal_flux_2d,only: xarray,zarray,psi_gradient,y2a_gradient,nx,nz
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: xval,zval
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_gradient,xval,zval,f)
  end function psi_gradient_func


  pure real(p_) function br(r,z) !R component of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in) :: r,z
    br=-psi_z_func(r,z)/r
  end function br

  pure real(p_) function bz(r,z) !Z component of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in) :: r,z
    bz=psi_r_func(r,z)/r
  end function bz

  pure real(p_) function bphi(r,z) !phi component of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    bphi=g_func(psi_func(r,z))/r
  end function bphi


  pure real(p_) function b(r,z) !strength of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in) :: r,z
    b=sqrt(br(r,z)**2+bz(r,z)**2+bphi(r,z)**2)
  end function b

  pure real(p_) function br_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    br_r=psi_z_func(r,z)/r**2-psi_rz_func(r,z)/r
  end function br_r


  pure real(p_) function br_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    br_z=-psi_zz_func(r,z)/r
  end function br_z

  pure real(p_) function bz_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_) , intent(in) :: r,z
    bz_r=-psi_r_func(r,z)/r**2+psi_rr_func(r,z)/r
  end function bz_r


  pure real(p_) function bz_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    bz_z=psi_rz_func(r,z)/r
  end function bz_z


  pure real(p_) function bphi_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    bphi_r=gprime(psi_func(r,z))*psi_r_func(r,z)/r-g_func(psi_func(r,z))/r**2
  end function bphi_r

  pure real(p_) function bphi_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one,zero
    implicit none
    real(p_), intent(in):: r,z
    bphi_z=gprime(psi_func(r,z))*psi_z_func(r,z)/r
  end function bphi_z

  pure real(p_) function b_r(r,z) !partial derivative of  magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    b_r=-b(r,z)/r+one/(b(r,z)*r*r)*(psi_r_func(r,z)*psi_rr_func(r,z) &
         & +psi_z_func(r,z)*psi_rz_func(r,z)+g_func(psi_func(r,z))*g_r(r,z))
  end function b_r

  pure real(p_) function b_z(r,z) !partial derivative of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in):: r,z
    b_z=one/(b(r,z)*r*r)*(psi_r_func(r,z)*psi_rz_func(r,z)+&
         & psi_z_func(r,z)*psi_zz_func(r,z)+g_func(psi_func(r,z))*g_z(r,z))
  end function b_z

  pure real(p_) function b_phi(r,z) result(f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    f=0._p_
  end function b_phi

  pure real(p_) function unitbr_r(r,z)  result(f)!partial derivative of the component of magnetic unit vector
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    bval=b(r,z)
    f=(br_r(r,z)*bval-b_r(r,z)*br(r,z))/bval**2
  end function unitbr_r

  pure real(p_) function unitbr_z(r,z) result(f)
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval

    bval=b(r,z)
    f=(br_z(r,z)*bval-b_z(r,z)*br(r,z))/bval**2
  end function unitbr_z


  pure real(p_) function unitbr_phi(r,z) result(f)
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval

    !bval=b(r,z)
    !f=(0*bval-b_phi(r,z)*br(r,z))/bval**2
    f=0.
  end function unitbr_phi


  pure real(p_) function unitbz_r(r,z) result(f) 
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_)::bval
    bval=b(r,z)
    f=(bz_r(r,z)*bval-b_r(r,z)*bz(r,z))/bval**2
  end function unitbz_r

  pure real(p_) function unitbz_z(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval

    bval=b(r,z)
    f=(bz_z(r,z)*bval-b_z(r,z)*bz(r,z))/bval**2
  end function unitbz_z

  pure real(p_) function unitbz_phi(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    !bval=b(r,z)
    !f=(0*bval-b_phi(r,z)*bz(r,z))/bval**2
    f=0.
  end function unitbz_phi


  pure real(p_) function unitbphi_r(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    bval=b(r,z)
    f=(bphi_r(r,z)*bval-b_r(r,z)*bphi(r,z))/bval**2
  end function unitbphi_r


  pure real(p_) function unitbphi_z(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    bval=b(r,z)
    f=(bphi_z(r,z)*bval-b_z(r,z)*bphi(r,z))/bval**2
  end function unitbphi_z


  pure real(p_) function g_r(r,z)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: psival
    psival=psi_func(r,z)
    g_r=gprime(psival)*psi_r_func(r,z)
  end function g_r

  pure real(p_) function g_z(r,z)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: psival
    psival=psi_func(r,z)
    g_z=gprime(psival)*psi_z_func(r,z)
  end function g_z

  pure real(p_) function g_func(psival) result(f) !all quantities are in S.I units.
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fpsi
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_), intent(in) :: psival
    
    call linear_1d_interpolate(npsi,psi_1d,fpsi,psival,f)  
  end function g_func

  pure real(p_) function gprime(psival) result(f)
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fprime
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_), intent(in):: psival
    call linear_1d_interpolate(npsi,psi_1d,fprime,psival,f)  
  end function gprime

  pure  real(p_) function qfunc(pfn) result(f) !safety factor, with correct sign
    use constants,only:p_
    use radial_module,only:npsi,pfn_npsi,q_with_sign
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_),intent(in) :: pfn
    
    call linear_1d_interpolate(npsi,pfn_npsi,q_with_sign,pfn,f)  
  end function qfunc

  pure real(p_) function qfunc0(psival) result(f) !safety factor
    !not used in computing the orbits, only used to do analytical estimation of some quantities, such as bounce frequency
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,qpsi
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_), intent(in):: psival
    call linear_1d_interpolate(npsi,psi_1d,qpsi,psival,f)  
  end function qfunc0


  pure real(p_) function radcor_as_func_of_pfn(pfn) result (z)
    use constants,only: p_
    implicit none
    real(p_), intent(in):: pfn
    z=pfn
  end function radcor_as_func_of_pfn

  pure real(p_) function tfn_func_pfn(pfn0) result(z)
    use constants,only:p_
    use radial_module, only : npsi, pfn_npsi, tfn_npsi
    use interpolate_module, only: linear_1d_interpolate
    implicit none
    real(p_),intent(in):: pfn0

    call linear_1d_interpolate(npsi, pfn_npsi, tfn_npsi, pfn0, z)  
  end function tfn_func_pfn


end module magnetic_field
