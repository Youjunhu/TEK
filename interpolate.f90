module interpolate_module
  use constants,only:p_
  use constants,only:one
  implicit none
  private
  public linear_2d_interpolation_kernel, linear_2d_interpolation, linear_2d_interpolation0, linear_1d_interpolation, &
       & linear_1d_interpolation_general_case, locate
contains


pure subroutine linear_2d_interpolation(nx,nz,xarray,zarray,psi,x,z,psival)  !uniform xarray and zarray are assumed
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in):: xarray(nx),zarray(nz),psi(nx,nz)
    real(p_),intent(in):: x,z
    real(p_),intent(out)::psival
    real(p_):: dx,dz,t1,t2,slope
    integer:: i,j ,ii,jj

    dx=xarray(2)-xarray(1)
    i=floor(one+(x-xarray(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval

    dz=zarray(2)-zarray(1)
    j=floor(one+(z-zarray(1))/dz)

    if(i.ge.nx) i=nx-1
    if(j.ge.nz) j=nz-1
    if(i.le.1) i=1 
    if(j.le.1) j=1 
    slope=(psi(i+1,j)-psi(i,j))/dx
    t1=psi(i,j)+slope*(x-xarray(i))
    slope=(psi(i+1,j+1)-psi(i,j+1))/dx
    t2=psi(i,j+1)+slope*(x-xarray(i))
    slope=(t2-t1)/dz
    psival=t1+slope*(z-zarray(j))
  end subroutine linear_2d_interpolation
  
  pure subroutine locate(m,x,dx,xval,i) !uniform grid is assumed
    use constants,only:p_
    implicit none
    integer, intent(in) :: m
    real(p_),intent(in) :: x(m),  dx,  xval
    integer,intent(out) :: i

    i = floor(one+(xval-x(1))/dx) 
    if(i==0) i=1
    if(i==m) i=m-1

  end subroutine locate


 pure subroutine linear_2d_interpolation0(nx,nz,xarray,zarray,dx,dz,psi,x,z,i,j,psival)  !uniform xarray and zarray are assumed
    use constants,only: p_, one
    implicit none
    integer,intent(in) :: nx,nz
    real(p_),intent(in) :: xarray(nx), zarray(nz), dx, dz, psi(nx,nz), x, z
    integer,intent(in) :: i,j
    real(p_),intent(out) :: psival
    real(p_) :: t1,t2,slope

    slope = (psi(i+1,j)-psi(i,j))/dx
    t1 = psi(i,j)+slope*(x-xarray(i))
    slope = (psi(i+1,j+1)-psi(i,j+1))/dx
    t2 = psi(i,j+1)+slope*(x-xarray(i))
    slope = (t2-t1)/dz
    psival = t1+slope*(z-zarray(j))
  end subroutine linear_2d_interpolation0
  
  pure  subroutine linear_2d_interpolation_kernel(x1a,x2a,ya,x1,x2,y)
    real(p_),intent(in) :: x1a(2), x2a(2), ya(2,2), x1, x2
    real(p_),intent(out) :: y
    real(p_) :: ytmp(2),slope
    integer :: j

    do j=1,2
       slope=(ya(2,j)-ya(1,j))/(x1a(2)-x1a(1))
       ytmp(j)=ya(1,j)+slope*(x1-x1a(1))
    enddo
    slope=(ytmp(2)-ytmp(1))/(x2a(2)-x2a(1))
    y=ytmp(1)+slope*(x2-x2a(1))

  end subroutine linear_2d_interpolation_kernel


  pure  subroutine linear_1d_interpolation(n,x,y,xval,yval)
    use constants,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval

    real(p_):: dx,slope
    integer:: i

    dx=x(2)-x(1)
    i=floor(one+(xval-x(1))/dx) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
    !if(i.le.0) write(*,*) 'warning****', xval, x(1),x(n),n
    if(i==0) i=1
    !call location(n,x,xval,i)
    if(i.ge.n) i=n-1
    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolation




  subroutine linear_1d_interpolation_general_case(n,x,y,xval,yval)  !non-uniform x array
    use constants,only:p_
    use constants,only:one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval
    real(p_):: slope
    integer:: i

    !dx=x(2)-x(1)
    !i=floor(one+(xval-x(1))/dx) !this for uniform x, otherwise we need to call location() subroutine to locate xval
    if(xval.ge.x(n) ) then
       i=n-1
    elseif(xval.le.x(1)) then
       i=1
    else
       call location(n,x,xval,i)
    endif

    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolation_general_case


  subroutine location(n,x,xval,k) !use bisection method to locate xval in an array
    !return k (xval is located between x(k) and x(k+1)
    use constants,only:p_
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),xval
    integer,intent(out)::k
    integer:: kl,ku,km

!!$    if(xval.gt.x(n) ) then
!!$       write(*,*) "***warning****, x provided is not in the range"
!!$       k=n-1
!!$       return
!!$    elseif(xval.lt.x(1)) then
!!$       write(*,*) "***warning****, x provided is not in the range"
!!$       k=1
!!$       return
!!$    endif
    kl=1
    ku=n
30  if(ku-kl .gt. 1) then  !use bisection method to search location of theta
       km=(ku+kl)/2
       if((x(n).ge.x(1)).eqv.(xval.ge.x(km))) then
          kl=km
       else
          ku=km
       endif
       goto 30
    endif
    k=kl
  end subroutine location



end module interpolate_module

