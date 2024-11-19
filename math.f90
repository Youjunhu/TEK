module math
contains

  function random_yj(seed) result (z) !return a random number uniform distribute in [0:1]
  !linear congruental method to genereate random number
  !This corresponds to Park and Miller's choice implemented using Schrange's algorithm to aviod integer overflow
  !refer to R. Fitzpatrick's book "Computational physics, an introduction course" for the details
  use constants,only:p_
  implicit none
  integer:: seed
  real(p_):: z 
!!$  integer, parameter:: a=106,c=1283,m=6075
!!$    z=mod(a*seed+c,m)
  integer, parameter:: a=16807,m=2147483647 !m=2^31-1,c=0, this choice is called Park and Miller method
  integer,parameter:: q=127773 !q=m/a
  integer,parameter:: r=2836 !r=mod(m,a)
  real(p_),parameter:: RANDMAX=2147483646._p_
  integer,save:: next=1
  !  write(*,*) m
  if (seed .ne.0) next=seed
  next=a*mod(next,q)-r*(next/q)
  if(next<0) next=next+m
  z=next/RANDMAX
end function random_yj


subroutine sub_random_yj(seed,next_seed,randnum)  !return a random number uniform distribute in [0:1] 
!modified form random_yj(), also return the seed for next generator which may be needed by another generator in another proc
  !linear congruental method to genereate random number
  !This corresponds to Park and Miller's choice implemented using Schrange's algorithm to aviod integer overflow
  !refer to R. Fitzpatrick's book "Computational physics, an introduction course" for the details
  use constants,only:p_
  implicit none
  integer,intent(in):: seed
  real(p_),intent(out):: randnum 
  integer,intent(out):: next_seed
  integer, parameter:: a=16807,m=2147483647 !m=2^31-1,c=0, this choice is called Park and Miller method
  integer,parameter:: q=127773 !q=m/a
  integer,parameter:: r=2836 !r=mod(m,a)
  real(p_),parameter:: RANDMAX=2147483646._p_
  integer,save:: next=1
  !  write(*,*) m
  if (seed .ne.0) next=seed
  next=a*mod(next,q)-r*(next/q)
  if(next<0) next=next+m
  randnum=next/RANDMAX
  next_seed=next
end subroutine sub_random_yj

  subroutine arc_length(x_contour,z_contour,np,dl)
    !calculate the poloidal arc length between neighbour points on a contour line.
    use constants,only: p_
    implicit none
    integer,intent(in) :: np
    real(p_),intent(in) :: x_contour(np), z_contour(np) 
    real(p_),intent(out) :: dl(np-1)
    integer :: i, i_plus_one, i_plus_two, i_minus_one
    real(p_) :: x(4), z(4)

    do i=1,np-1
       i_plus_one=i+1  !i_plus_one indicates the right point
       i_minus_one=i-1 !i_minus_one indicates the left point
       i_plus_two=i+2
       if (i == 1)  i_minus_one=np-1 !deal with boundary points
       if (i == np-1) i_plus_two=2 !deal with boundary points

       x(1)=x_contour(i_minus_one)
       x(2)=x_contour(i)
       x(3)=x_contour(i_plus_one)
       x(4)=x_contour(i_plus_two)

       z(1)=z_contour(i_minus_one)
       z(2)=z_contour(i)
       z(3)=z_contour(i_plus_one)
       z(4)=z_contour(i_plus_two)

       !call arc_between_two_points(x,z,dl(i))

       dl(i)=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2) !use simple formula to calculate the arc length
    enddo

  end subroutine arc_length

  subroutine   arc_between_two_points(x,z,dl)
    !calculate the arc length between point (x(2),z(2)) and point (x(3),z(3))
    use constants,only:p_
    use constants,only:one,two
    implicit none
    real(p_),intent(in):: x(4),z(4)
    real(p_),intent(out):: dl
    real(p_):: ds,a_length,b_length
    real(p_):: dot_a_and_ds,dot_b_and_ds, cos_tha,cos_thb,m1,m2

    !ds is the length of straight-line segment passing through (x(2),z(2)) and (x(3),z(3))
    ds=sqrt((x(3)-x(2))**2+(z(3)-z(2))**2) 

    a_length=sqrt((x(3)-x(1))**2+(z(3)-z(1))**2)
    b_length=sqrt((x(4)-x(2))**2+(z(4)-z(2))**2)
    dot_a_and_ds=(x(3)-x(1))*(x(3)-x(2)) &
         +(z(3)-z(1))*(z(3)-z(2))
    dot_b_and_ds=(x(4)-x(2))*(x(3)-x(2)) &
         +(z(4)-z(2))*(z(3)-z(2))
    cos_tha=dot_a_and_ds/(a_length*ds)
    cos_thb=dot_b_and_ds/(b_length*ds)

    m1=sqrt(one-cos_tha**2)/cos_tha 
    m2=sqrt(one-cos_thb**2)/cos_thb
    !the value of m1 and m2 should be positive for most cases
    dl=ds*(one+(two*m1**2+two*m2**2+m1*m2)/30._p_) !calculate arc length using Eq. (5.38) in  S. Jardin's book. Here I assume that the dot product of the slope is negative, so the -m1*m2 term is replaced with +abs(slope1)*abs(slope2), need checking the correctness for general case
    !dl=ds*(1._p_+0.) !use linear function to approximate the arc length
  end subroutine arc_between_two_points

  subroutine ZGETRS_wrapper(kn,matrix, ipiv, rhs_dft, solution_dft) !solve the field equation for a single toroidal harmonic
    use constants,only:p_
    implicit none
    integer,intent(in):: kn !index of the toroidal harmonic
    complex(p_), intent(in) :: matrix(:,:,0:) !LU factorization of the radial coefficient matrix
    integer, intent(in) :: ipiv(:,0:) 
    complex(p_),intent(in):: rhs_dft(:)
    complex(p_),intent(out):: solution_dft(:)
    complex(p_):: s(size(rhs_dft),1)
    integer:: info, n

    n = size(rhs_dft)
    s(:,1) = rhs_dft(:)
    call ZGETRS("N", n, 1, matrix(:,:,kn), n, IPIV(:,kn), s, n, INFO) !solve using the LU factorization computed by ZGETRF.
    !call svd_back_substitution(s(:,1),kn) !solve using the SVD computed by ZGESDD, the results agree with that given by ZGETRS.
    solution_dft=s(:,1)
  end subroutine ZGETRS_wrapper


  subroutine srcbes(biz,gam0,gam1)
    use constants,only:p_
    implicit none
    REAL(p_) :: t1,t2,biz,gam0,gam1
    !.....Calculates gamma nought and gamma 1. (Abramowitz and Stegun).
446 if (biz.gt.3.75d0) go to 148
    t1=(biz/3.75d0)**2
    t2=exp(-biz)
    gam0=t2*((((((.0045813d0*t1+.0360768d0)*t1+.2659732d0)*t1+ &
         1.2067492d0)*t1+3.0899424d0)*t1+3.5156229d0)*t1+1.d0)
    gam1=t2*biz*((((((.00032411d0*t1+.00301532d0)*t1+.02658733d0) &
         *t1+.15084934d0)*t1+.51498869d0)*t1+.87890594d0)*t1+.5d0)
    go to 149
148 t2=1.d0/sqrt(biz)
    t1=3.75d0/biz
    gam0=t2*((((((((.00392377d0*t1-.01647633d0)*t1+.02635537d0) &
         *t1-.02057706d0)*t1+.00916281d0)*t1-.00157565d0)*t1+ &
         .00225319d0)*t1+.01328592d0)*t1+.39894228d0)
    gam1=t2*((((((((-.00420059d0*t1+.01787654d0)*t1-.02895312d0) &
         *t1+.02282967d0)*t1-.01031555d0)*t1+.00163801d0)*t1- &
         .00362018d0)*t1-.03988024d0)*t1+.39894228d0)
149 continue
    return
  end subroutine srcbes

  !I have verified that gam0 obtained by calling srcbes(x,gam0,gam1) is equal to bessi0(x)*exp(-x)
  FUNCTION bessi0(x) ! the 0th modified Bessel function of the first kind.
    use constants,only:p_
    implicit none
    REAL(p_)::  bessi0,x
    REAL(p_):: ax
    real(p_):: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0, &
         & 1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
    DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,&
         & 0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1, &
         & -0.1647633d-1,0.392377d-2/
    if (abs(x).lt.3.75) then
       y=(x/3.75)**2
       bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
    else
       ax=abs(x)
       y=3.75/ax
       bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y* &
            &(q7+y*(q8+y*q9))))))))
    endif
  END FUNCTION bessi0
  !  (C) Copr. 1986-92 Numerical Recipes Software .

  subroutine cross_product_in_cartesian(ax,ay,az,bx,by,bz,cx,cy,cz)
    use constants,only:p_
    implicit none

    real(p_),intent(in):: ax,ay,az,bx,by,bz
    real(p_),intent(out)::cx,cy,cz

    cx=ay*bz-az*by
    cy=az*bx-ax*bz
    cz=ax*by-ay*bx
  end subroutine cross_product_in_cartesian


  subroutine shift_to_zero_twopi_range(a) !shift a into the range [0:twopi]
    use constants,only:p_
    use constants,only: twopi
    implicit none
    real(p_),intent(inout):: a
    integer:: ishift
!!$  a=a-int(a/twopi)*twopi !shift into the range [0:twopi]
!!$  if(a.lt.0) a=a+twopi !shift into the range [0:twopi]

    ishift=floor(a/twopi)
    a=a-ishift*twopi

  end subroutine shift_to_zero_twopi_range


  elemental subroutine shift_to_specified_toroidal_range(a) !shift "a" into the range [0:toroidal_range]
    use constants,only:p_
    !  use constants,only: twopi
    use magnetic_coordinates,only: toroidal_range !is positive
    implicit none
    real(p_),intent(inout):: a
    integer:: ishift

!!$  a=a-int(a/toroidal_range)*toroidal_range !shift into the range [0:toroidal_range]
!!$  if(a.lt.0) a=a+toroidal_range !shift into the range [0:toroidal_range]
!!$
!!$ ishift=floor(a/toroidal_range)
!!$ a=a-ishift*toroidal_range

    a = mod(a,toroidal_range)
    if(a<0) a = a + toroidal_range
  end subroutine shift_to_specified_toroidal_range


  subroutine shift_to_minus_pi_positive_pi_range(a) !shift "a" into the range [-pi:pi]
    use constants,only:p_
    use constants,only: twopi,pi
    implicit none
    real(p_),intent(inout):: a
    integer:: ishift

    ishift=floor(a/twopi)
    a=a-ishift*twopi
    if(a>pi) a=a-twopi
  end subroutine shift_to_minus_pi_positive_pi_range


  subroutine one_dimensional_derivative(n,x,y,dydx)
    use constants,only:p_
    use constants,only:zero,one,two,twopi,one_half
    implicit none

    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(out):: dydx(n)
    real(p_):: tmp0,dx
    integer:: j

    dx=x(2)-x(1) !uniform interval is assumed

    do j=2,n-1 !use center difference scheme for inner points
       dydx(j)=(y(j+1)-y(j-1))/(two*dx)
    enddo

    !use linear interpolation to get the value  j=n
    tmp0=(y(n)-y(n-1))/dx
    dydx(n)=two*tmp0-dydx(n-1)
    !use linear interpolation to get the value j=1
    tmp0=(y(2)-y(1))/dx
    dydx(1)=two*tmp0-dydx(2)
  end subroutine one_dimensional_derivative


  subroutine calculate_determinant_complex(a,n,d)
    use constants ,only: p_
    implicit none
    integer:: n
    complex(p_)::a(n,n),d  
    integer::info,j,ipiv(n)

    ! Compute LU Factorization
    if ( p_ == kind(1.0e1) ) then   !cgetrf and zgetrf are programs in Lapack
       call cgetrf(n,n,a,n,ipiv,info)
    else
       call zgetrf(n,n,a,n,ipiv,info)
    endif

    !  compute determinant
    if ( info .ge. 0 ) then
       d = (1.0,0.0)
       do j=1,n
          d = d*a(j,j)
          if ( ipiv(j) /= j ) d=-d
       enddo
    else
       print *," *** Error in computing determinant"
       print *," info =",info
    endif
  end subroutine calculate_determinant_complex

  !C PNP1                                                                
  !C       http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html                                                                
  !C     ..................................................................
  !C                                                                       
  !C        SUBROUTINE PNPOLY                                              
  !C                                                                       
  !C        PURPOSE                                                        
  !C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
  !C                                                                       
  !C        USAGE                                                          
  !C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
  !C                                                                       
  !C        DESCRIPTION OF THE PARAMETERS                                  
  !C           PX      - X-COORDINATE OF POINT IN QUESTION.                
  !C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
  !C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
  !C                     VERTICES OF POLYGON.                              
  !C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
  !C                     VERTICES OF POLYGON.                              
  !C           N       - NUMBER OF VERTICES IN THE POLYGON.                
  !C           INOUT   - THE SIGNAL RETURNED:                              
  !C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
  !C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
  !C                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
  !C                                                                       
  !C        REMARKS                                                        
  !C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
  !C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
  !C           OPTIONALLY BE INCREASED BY 1.                               
  !C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
  !C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
  !C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
  !C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
  !C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
  !C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
  !C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
  !C                                                                       
  !C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
  !C           NONE                                                        
  !C                                                                       
  !C        METHOD                                                         
  !C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
  !C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
  !C           POINT IS INSIDE OF THE POLYGON.                             
  !C                                                                       
  !C     ..................................................................
  !C                                                                       

  pure  SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)
    use constants,only:p_
    implicit none
    integer,intent(in):: n
    integer,intent(inout):: inout
    real(p_),intent(in):: px,py
    REAL(p_),intent(in):: XX(N),YY(N)                                    
    integer:: i,j,maxdim
    REAL(p_):: X(n),Y(n), tmp
    LOGICAL MX,MY,NX,NY                                               
    !      INTEGER O                                                         
    !      OUTPUT UNIT FOR PRINTED MESSAGES                                 
    !      DATA O/6/                                                         

    !yj    MAXDIM=200
    MAXDIM=n
    IF(N.LE.MAXDIM)GO TO 6                                            
    !      WRITE(O,7)                                                        
7   FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. RESULTS INVALID')     
    RETURN                                                            
6   DO  I=1,N                                                        
       X(I)=XX(I)-PX                                                     
       Y(I)=YY(I)-PY
    enddo
    INOUT=-1                                                          
    DO I=1,N                                                        
       J=1+MOD(I,N)                                                      
       MX=X(I).GE.0.0                                                    
       NX=X(J).GE.0.0                                                    
       MY=Y(I).GE.0.0                                                    
       NY=Y(J).GE.0.0                                                    
       IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) cycle
       IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
       INOUT=-INOUT                                                      
       cycle
!!$3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
!!$4     INOUT=0                                                           
!!$      RETURN                                                            
!!$5     INOUT=-INOUT                                                      
!!$2     CONTINUE                                                          

3      tmp=(Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))
       IF(tmp<0) then
          cycle
       elseif(tmp==0) then
          goto 4
       else
          goto 5
       endif
4      INOUT=0                                                           
       RETURN                                                            
5      INOUT=-INOUT                                                      
    enddo
  END SUBROUTINE PNPOLY


  subroutine laplace_cylindrical2d(psi, x, z, nx, nz, jphi)
    !  use math, only : partial_derivatives_2d
    use constants, only : p_, mu0
    implicit none
    integer, intent(in) :: nx, nz
    real(p_), intent(in) :: psi(nx,nz),  x(nx), z(nz)
    real(p_), intent(out) :: jphi(nx,nz)
    real(p_) :: psi_x(nx,nz), psi_z(nx,nz), psi_xx(nx,nz), psi_xz(nx,nz), psi_zx(nx,nz),psi_zz(nx,nz)
    real(p_) :: dx, dz
    integer :: i, j  

    dx = x(2) - x(1)
    dz = z(2) - z(1)
    call partial_derivatives_2d(nx,nz,x,z, psi,   psi_x,  psi_z)
    call partial_derivatives_2d(nx,nz,x,z, psi_x, psi_xx, psi_xz)
    call partial_derivatives_2d(nx,nz,x,z, psi_z, psi_zx, psi_zz)

    do i= 1, nx
       do j=1,nz
          jphi(i,j) = psi_zz(i,j) + psi_xx(i,j)  - 1/x(i)*psi_x(i,j)
          jphi(i,j) = -jphi(i,j)/(mu0*x(i)) !to toroidal current density
       enddo
    enddo

  end subroutine laplace_cylindrical2d

  subroutine partial_derivatives_2d(nx,nz,rarray,zarray,b,b_x,b_z)
    use constants,only:p_
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in):: rarray(nx),zarray(nz),b(nx,nz)
    real(p_),intent(out):: b_x(nx,nz),b_z(nx,nz)
    integer:: i,j,i1,j1,i2,j2

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
          b_x(i,j)=(b(i2,j)-b(i1,j))/(rarray(i2)-rarray(i1))
          b_z(i,j)=(b(i,j2)-b(i,j1))/(zarray(j2)-zarray(j1))
       enddo
    end do
  end subroutine partial_derivatives_2d


  subroutine arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,x_axis,z_axis,myid)
    !replacing one point near the high-field side of the midplane by a point that is exactly on the high-field side of the midplane
    !then arrange the arrays x_lcfs_new and z_lcfs_new so that the starting point (x_lcfs_new(1),z_lcfs_new(1)) is on the high-field-side of the midplane
    !midplane is defined as the z=z_axis plane
    use constants,only:p_
    implicit none
    integer,intent(in):: np_lcfs,myid
    real(p_),intent(inout):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
    real(p_),intent(in):: x_axis,z_axis
    real(p_):: x_lcfs_new(np_lcfs),z_lcfs_new(np_lcfs) 
    real(p_):: Rout
    real(p_):: r_major,r_minor,eps,direction
    integer:: kk(1),k,i

    !write(*,*) ' Z values of the uppermost point of LCFS: ', maxval(z_lcfs)
    !write(*,*) ' Z values of the lowest point of LCFS: ' ,minval(z_lcfs)

    !set the starting point of LCFS to be at the low-field-side of the midplane
    !maxval(x_lcfs)
    !kk=maxloc(x_lcfs) !return the index of the array for which R is the largest,in order to determine the low-field side of the midplane
    kk=minloc(x_lcfs) !return the index of the array for which R is the smallest,in order to determine the high-field side of the midplane
    k=kk(1)
    !k=10
    !write(*,*) 'index of the point on the lcfs that have the largest R, k=',k
    !if((z_lcfs(k+1)-z_axis)*(z_lcfs(k-1)-z_axis)>0) stop 'error in selecting the point on LCFS'

    call major_radius_on_midplane(np_lcfs,x_lcfs,z_lcfs,x_axis,z_axis,Rout)

    !  r_major=(maxval(x_lcfs)+minval(x_lcfs))/2._p_ !by definition
    !  r_minor=(maxval(x_lcfs)-minval(x_lcfs))/2._p_ !by definition

    ! write(*,*) 'inverse aspect ratio of LCFS is (definied at low-field side) ', (Rout-x_axis)/x_axis
    !  write(*,*) 'r_axis=',x_axis, 'r_major=', r_major, 'r_minor=',r_minor
    !  eps= r_minor/r_major !standard definition
    ! write(*,*) 'inverse aspect ratio of LCFS (i.e., r_minor/r_major) is ', eps
    ! write(*,*) 'ellipticity (elongation) of LCFS is ', (maxval(z_lcfs)-minval(z_lcfs))/2._p_/r_minor
    !write(*,*) 'upper triangularity of LCFS is ', (r_major-x_lcfs(maxloc(z_lcfs)))/r_minor, &
    !          & 'lower triangularity of LCFS is ', (r_major-x_lcfs(minloc(z_lcfs)))/r_minor
    !replace one point of LCFS with the new point
    x_lcfs(k)=Rout
    z_lcfs(k)=z_axis

    !arrange the arrays so that x_lcfs_new and z_lcfs_new start from the low/high-field-side of the midplane
    do i=1,np_lcfs
       if(k+i-1.le.np_lcfs) then
          x_lcfs_new(i)=x_lcfs(k+i-1)
          z_lcfs_new(i)=z_lcfs(k+i-1)
       else
          x_lcfs_new(i)=x_lcfs(k+i-np_lcfs)
          z_lcfs_new(i)=z_lcfs(k+i-np_lcfs)
       endif
    enddo

    !use x_lcfs and z_lcfs to store the new data
    x_lcfs=x_lcfs_new
    z_lcfs=z_lcfs_new

    !check wheter the direction of the sequecne (r(i),z(i)) with i increasing is clockwise or anticlockwise when viewed along grad_phi direction, if clockwise, switch it to anticlockwise
    !This is achieved by using the determination of the direction matrix (a well known method in graphic theory).
    !Because the contours of Psi considered here are always convex polygons (instead of concave polygons), we can select any vertex on the curve to calculate the direction matrix. (refer to wikipedia about the direction matrix)
    direction=(x_lcfs(2)-x_lcfs(1))*(z_lcfs(3)-z_lcfs(1))-(x_lcfs(3)-x_lcfs(1))*(z_lcfs(2)-z_lcfs(1))
    if(direction .lt. 0.) then
       if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) &
            & with i increasing is clockwise, switch it to anticlockwise'

       do i=1,np_lcfs !switch it to anticlockwise
          x_lcfs(i)=x_lcfs_new(np_lcfs+1-i)
          z_lcfs(i)=z_lcfs_new(np_lcfs+1-i)
       enddo
    else if (direction .gt. 0.) then
       if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) with i increasing is anticlockwise'
    else
       stop 'the three vertex (points) used in calculating the direction matrix is collinear'
    endif


    if (myid.eq.0) then
       open(123,file='lcfs2.txt')
       do i=1,np_lcfs
          write(123,*) x_lcfs(i),z_lcfs(i)
       enddo
       close(123)
    endif

  end subroutine arrange_lcfs



  subroutine major_radius_on_midplane(mpol,rs,zs,r_axis,z_axis,Rout)
    !given a flux surface, this subroutine determines the major radius of the point on the low/high-field side of the middle plane
    use constants,only:p_
    use interpolate_module,only: linear_1d_interpolation_general_case

    implicit none
    integer,intent(in):: mpol
    real(p_),intent(in):: rs(mpol), zs(mpol),r_axis,z_axis
    real(p_),intent(out):: Rout !the major radius of the point on the low/high-field-side of the mid-plane
    integer:: i,j,k1,k2,n
    real(p_):: r_select(mpol),z_select(mpol)
    !real(p_),dimension(:),allocatable:: x,z,tmp_y2
    real(p_):: tmp

    n=1
    do i=1,mpol-1 
       !     if(rs(i).gt.r_axis) then !select the low-field side (i.e., the part with r larger than r_axis) of a flux surface
       if(rs(i).lt.r_axis) then !select the high-field side (i.e., the part with r less than r_axis) of a flux surface
          r_select(n)=rs(i)
          z_select(n)=zs(i)
          n=n+1
       endif
    enddo
    !  write(*,*) 'n-1= ', n-1
    !order the array according to the value of z_select
    do i=1,n-1
       do j=i+1,n-1
          if(z_select(j).le.z_select(i)) then
             !exchange the z value 
             tmp=z_select(i)
             z_select(i)=z_select(j)
             z_select(j)=tmp
             !also exchange the r value (I forgot this step in the older version, which cause a serious error)
             tmp=r_select(i)
             r_select(i)=r_select(j)
             r_select(j)=tmp
          endif
       enddo
    end do

    call linear_1d_interpolation_general_case(n-1,z_select,r_select,z_axis,Rout)  

  end subroutine major_radius_on_midplane

end module math

