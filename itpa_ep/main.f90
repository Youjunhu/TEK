module  precision
  implicit none
  !integer,parameter:: p_=kind(1.0)
  integer,parameter:: p_=kind(1.0d0)
end module precision


module constants
  use precision,only:p_
  implicit none

  real(p_),parameter:: pi=3.141592653589793_p_
  real(p_),parameter:: twopi=pi*2.0_p_
  real(p_),parameter:: fourpi=pi*4.0_p_
  real(p_),parameter:: mu0=fourpi*1.0d-7 !permeability in SI unit
  real(p_),parameter:: zero=0.0_p_
  real(p_),parameter:: one=1.0_p_
  real(p_),parameter:: two=2.0_p_
  real(p_),parameter:: three=3.0_p_
  real(p_),parameter:: one_half=0.5_p_
  real(p_),parameter:: one_third=one/three
  real(p_),parameter:: one_fifth=0.2_p_
  real(p_),parameter:: three_halfs=1.5_p_
  real(p_),parameter:: four=4._p_
  real(p_),parameter :: kev=1.6022d-16   !unit J
end module constants


module magnetic_parameters 
  use precision,only:p_
  implicit none
  !real(p_),parameter:: a=0.60_p_ !minor radius in meter
  real(p_),parameter:: a = 1.0_p_ !minor radius in meter
  !real(p_),parameter::r0=1.67_p_ !major radius in unit of meter
  real(p_),parameter:: R0 = 10.0_p_
  !real(p_),parameter::b0 = 2.0_p_ !in Tesla
  real(p_),parameter:: b0 = 3.0_p_ !in Tesla
  real(p_),parameter:: g0 = r0*b0
  real(p_),parameter:: minor_r0 = a*0.5_p_

end module magnetic_parameters


module equilibrium
  contains
subroutine calculate_equilibrium()
  use precision,only:p_
  use constants,only:two,twopi, pi, kev
  use magnetic_parameters,only:r0,a,g0, minor_r0
  implicit none
  integer,parameter:: nx=257, nz=257
  integer,parameter:: np_lcfs = 257, n=nx
  real(p_):: psi_2d(nx,nz)
  real(p_):: r_1d(nx),z_1d(nz)
  real(p_):: r_lcfs(np_lcfs),z_lcfs(np_lcfs)
  real(p_):: f(n), ffp(n), q(n), p(n), pp(n), pfn(n) !on uniform poloidal flux grid
  integer:: i,j,u
  real(p_)::  rleft,  zmid,  rdim,  zdim !specify the rectangular region for which the value of psi_2d is specified
  real(p_):: psi_axis,psi_bdry,xmaxis,zmaxis,psival
  real(p_):: minor_r ,theta !magnetic surface coordinates

  real(p_) :: ti(n), ni(n), kappa_ti, kappa_ni
  real(p_),parameter:: epsilon=1.0d-6, xacc=1.0d-6


  do i=1,np_lcfs !the shape of the LCFS
     theta = -pi+(i-1)*twopi/(np_lcfs-1)
     r_lcfs(i) = r0 + a*cos(theta)
     z_lcfs(i) = a*sin(theta)
  enddo

  open(11,file='lcfs.txt') !store the location of LCFS
  do i=1,np_lcfs
     write(11,*) r_lcfs(i),z_lcfs(i)
  enddo
  close(11)

  !set a rectangular region
  rleft = R0 - 1.1_p_*a
  rdim = 2.2*a
  zmid = 0.0_p_
  zdim = 2.2*a

  do i=1,nx
     r_1d(i) = rleft + rdim/(nx-1)*(i-1)
  enddo
  do j=1,nz
     z_1d(j) = (zmid-zdim/two) + zdim/(nz-1)*(j-1)
  enddo

  call draw_rect_region(nx,nz,r_1d,z_1d)

  open(112,file='psi_2d.txt')
  do i=1,nx !generate 2d poloidal magnetic flux data within the rectangular region on (r,z) plane
     do j=1,nz
        minor_r = sqrt((r_1d(i)-r0)**2+(z_1d(j)-0._p_)**2)
        psi_2d(i,j) = psi_func(minor_r)
        write(112,*) r_1d(i), z_1d(j), psi_2d(i,j)
     enddo
     write(112,*) ! for gnuplot grids data format
  enddo
  close(112)

  psi_bdry = psi_func(a)
  psi_axis = psi_func(epsilon)
  write(*,*) 'psi_bdry=',psi_bdry,'psi_axi=',psi_axis

  kappa_ti = 6.96/R0
  kappa_ni = 2.23/R0
  do j=1,n
     !minor_r=epsilon+(a-epsilon)/(n-1)*(j-1)
     psival = psi_axis + (psi_bdry-psi_axis)/(n-1)*(j-1)
     pfn(j) = (psival - psi_axis)/(psi_bdry -psi_axis)
     minor_r = rtbis(func0, epsilon*0.01_p_, a*1.01_p_, psival, xacc) !find the minor_r where poloida_flux=psival
     q(j) = q_func(minor_r) !q array corresponding to uniform poloidal flux array
     Ti(j)= 1.0d0
     ni(j)= 2.0d19

  enddo

  f = g0 !Bphi*R is a constant
  ffp = 0._p_ !Bphi*R is a constant
  !p = 0._p_ !zero pressure
   p = 2*Ti*kev*ni
  pp = 0._p_ 
  xmaxis = r0
  zmaxis = 0.

  call  write_gfile(psi_2d, nx, nz, rleft, zmid, rdim, zdim, psi_axis, psi_bdry, xmaxis, zmaxis, &
       & r_lcfs, z_lcfs, np_lcfs, f, q, p, pp, ffp)

  open(newunit=u,file='ti.txt')
  do j=1,n
     write(u,*)  pfn(j), Ti(j)
  enddo
  close(u)

  open(newunit=u,file='ni.txt')
  do j=1,n
     write(u,*)  pfn(j), ni(j)
  enddo
  close(u)

  
call tor_shift()

end subroutine calculate_equilibrium

subroutine write_gfile(psi,nx,nz,rleft,zmid,xdim,zdim,psimax,psi_lcfs,xmaxis,zmaxis,&
     & r_lcfs,z_lcfs,np_lcfs,fpsi,qpsi,presspsi,pprime,ffprime)
  use precision,only:p_
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: psi(nx,nz)
  character(len=8):: ntitle(6)
  integer:: neq,ipestg
  real(p_),intent(in)::xdim, zdim,  rleft, zmid
  real(p_),intent(in)::xmaxis, zmaxis, psimax, psi_lcfs

  real(p_):: dumaraya5(5),dumarayb5(5)
  real(p_),intent(in):: fpsi(nx),qpsi(nx),presspsi(nx),ffprime(nx),pprime(nx)
  integer,intent(in):: np_lcfs
  real(p_),intent(in):: r_lcfs(np_lcfs),z_lcfs(np_lcfs)
  integer,parameter:: nlim_eqd=100
  real(p_):: rlim_eqd(nlim_eqd), zlim_eqd(nlim_eqd) !no data, just a household array
  real(p_):: rcenter,btorus !no data, just a household array
  integer:: i,j
  rlim_eqd=0.
  zlim_eqd=0.
  rcenter=1.0 !arbitrary data
  btorus=1.0 !arbitrary data
  ntitle=''
  ntitle(1)='circular'
  ipestg=22 !arbitrary number
  dumaraya5=0;dumarayb5=0 !arbitrary number
  neq=111
  open(neq,file='gfile_itpa_ep257x257')
  write (neq, '(6a8, 3i4)') (ntitle(i), i=1,6), ipestg, nx, nz
  write (neq,300) xdim, zdim, rcenter, rleft, zmid
  write (neq,300) xmaxis, zmaxis, psimax, psi_lcfs, btorus
  write (neq,300) dumaraya5
  write (neq,300) dumarayb5
  write (neq ,300) (fpsi(i), i=1,nx)
  write (neq ,300) (presspsi(i), i=1,nx)
  write (neq ,300) (ffprime(i), i=1,nx)
  write (neq ,300) (pprime(i), i=1,nx)
  write (neq ,300) ((psi(i,j), i=1,nx), j=1,nz)
  write (neq ,300) (qpsi(i), i=1,nx)
  write (neq ,'(2i5)') np_lcfs, nlim_eqd
  write (neq ,300) (r_lcfs(i), z_lcfs(i), i=1,np_lcfs)
  write (neq ,300) (rlim_eqd(i), zlim_eqd(i), i=1,nlim_eqd)
  close(neq)

300 format (5e16.9)

  write(*,*) 'xdim=',xdim, 'zdim=',zdim, 'rleft=',rleft,'zmid=',zmid
  write(*,*) 'xmaxis=',xmaxis, 'zmaxis=',zmaxis,'psimax=',psimax,'psi_lcfs=',psi_lcfs
  write(*,*) 'np_lcfs=',np_lcfs
  !write(*,*) 'r_lcfs(1),z_lcfs(1),r_lcfs(np_lcfs),z_lcfs(np_lcfs)=', r_lcfs(1),z_lcfs(1),r_lcfs(np_lcfs),z_lcfs(np_lcfs)
!  EFITD    04/03/2009    # 13606  7104ms           3 129 129
end subroutine write_gfile


subroutine tor_shift() !using analytic formula to compute the toroidal shift used in the defintion of generalized toroidal angle
  use precision,only:p_
  use constants,only:two,twopi,pi
  use magnetic_parameters,only:r0,a
  implicit none
  real(p_),parameter:: epsilon=1.0d-6,xacc=1.0d-6
  integer,parameter:: m=50,n=21
  integer::i,j
  real(p_):: minor_r,r,theta

  real(p_)::pfn_inner=0.4d0,pfn_bdry=0.7d0
  real(p_):: psi_inner,psi_outer
  real(p_):: psi_axis,psi_bdry,psival

  psi_bdry=psi_func(a)
  psi_axis=psi_func(epsilon)

  psi_inner=psi_axis+(psi_bdry-psi_axis)*pfn_inner
  psi_outer=psi_axis+(psi_bdry-psi_axis)*pfn_bdry

  open(113,file='tor_shift.txt')
  do j=1,n
     psival=psi_inner+(psi_outer-psi_inner)/(n-1)*(j-1)
     minor_r=rtbis(func0,epsilon*0.01_p_,a*1.01_p_,psival,xacc) !find the minor_r where poloida_flux=psival
     r=minor_r
     do i=1,m
        theta=-pi+twopi/(m-1)*(i-1)
        write(113,*) r,theta,2*q_func(r)*atan(tan(theta/two)*(r0-r)/sqrt(r0**2-r**2))
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)
end subroutine tor_shift

function func0(minor_r,psival) result(z)
  use precision,only:p_
  implicit none
  real(p_):: z,psival,minor_r

z=psi_func(minor_r)-psival

end function func0

function psi_func(minor_r) result(z)
  use precision,only:p_
  use constants,only: twopi,one_half
  use magnetic_parameters,only:r0,g0
  implicit none
  real(p_):: minor_r, z
  integer, parameter:: n = 400
  integer:: j
  real(p_):: sum, xmid, low_limit, upp_limit, step

  low_limit = 1.0d-6
  upp_limit = minor_r
  step = (upp_limit-low_limit)/(n-1)

  sum = 0._p_
  do j = 1, n-1
     xmid = low_limit + (j-1)*step + one_half*step
     sum = sum + twopi*g0*xmid/(q_func(xmid)*sqrt(r0**2-xmid**2))*step
  enddo
 z = sum/twopi
end function psi_func

function q_func(minor_r) result(z)
  use precision,only:p_
  use magnetic_parameters,only: minor_r0, a
  implicit none
  real(p_):: minor_r, z, s
  !real(p_),parameter:: q0=1.41_p_
    !real(p_),parameter:: shear0=0.78_p_
  !real(p_),parameter:: shear0=0.84_p_

  !real(p_):: qprime0
!z=q0*(minor_r/minor_r0)**shear0 !q profile with a constant shear
!!$qprime0=shear0*q0/minor_r0
!!$z=q0+(minor_r-minor_r0)*qprime0 !a linear q profile over r, with value of q and the shear at r=minor_r0 being q0 and shear0, respectively. This profile is used in Ben's toroidal ITG simulation
s = minor_r/a
!z  =0.86_p_-0.16_p_*s+2.52_p_*s*s !DIII-D cyclone case
z = 1.71 + 0.16*s**2 !ITPA-EP TAE benchmarking

!z=1.41+(s-r0a)*(0.35*1.41/r0a) 
end function q_func


subroutine draw_rect_region(nx,nz,r_1d,z_1d)
  use precision,only:p_
!  use constants,only:two,twopi !,four
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

FUNCTION rtbis(func,x1,x2,psival,xacc)
  !find a root of func by using the bisection method
  use precision,only: p_
  implicit none
  INTEGER JMAX
  REAL(p_) rtbis,x1,x2,psival,xacc,func
  EXTERNAL func
  PARAMETER (JMAX=40)
  INTEGER j
  REAL(p_) dx,f,fmid,xmid
  fmid=func(x2,psival)
  f=   func(x1,psival)
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
     fmid=func(xmid,psival)
     if(fmid.le.0.)rtbis=xmid
     if(abs(dx).lt.xacc) return
  enddo
  stop 'too many bisections in rtbis'
end FUNCTION rtbis
end module equilibrium



program main
  use equilibrium, only: calculate_equilibrium
  implicit none

  call calculate_equilibrium()

end  program main
