module fourn_module
contains
  subroutine twod_fourier_transform_nr(s,s_fft,m,n) !tested on Feb 27, 2018, note that the convension of FFT and inverse FFT definition is differnt between FFTW and Numerical-recipe book. I forgot about this when comparing the result given by this subroutine with that by FFTW subroutine, and it took about one hour for me to finally figure out why they are different.
    use constants,only:p_
!    use transform_module,only: my_fft
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(in):: s(m,n)
    complex(p_),intent(out):: s_fft(m,n)
    real(p_):: data(2*m*n)
    integer,parameter:: ndim=2
    complex(p_),parameter:: ii=(0._p_,1._p_)
    complex(p_):: in(m,n),out(m,n)
    integer:: i,j,k,nn(ndim),isign

    do j=1,n
       do i=1,m
          k=i+m*(j-1)
          data(2*k-1)=s(i,j)  !real part
          data(2*k)=0._p_ !imaginary part
       enddo
    enddo

!---for testing--------pass this testing, Feb. 27, 2018
!!$    in=s
!!$    do j=1,n
!!$       call my_fft(in(:,j),out(:,j),m)
!!$    enddo
!!$    !Fourier transformation along the second dimension
!!$    in=out
!!$    do i=1,m
!!$       call my_fft(in(i,:),out(i,:),n)
!!$    enddo
!!$    s_fft=out
!----for testing----
    nn(1)=m
    nn(2)=n
    isign=-1 !isign=-1 corresponds to inverse FFT in numerical-recipe-book, but corresponds to the forward FFT in FFTW
    call fourn(data,nn,ndim,isign)

    do j=1,n
       do i=1,m
          k=i+m*(j-1)
          s_fft(i,j)=data(2*k-1)+ii*data(2*k)
       enddo
    enddo
 end subroutine twod_fourier_transform_nr


subroutine twod_inverse_fourier_transform_nr(s_fft,s,m,n) 
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  complex(p_),intent(in):: s_fft(m,n)
  real(p_),intent(out):: s(m,n)
  real(p_):: data(2*m*n)
  integer,parameter:: ndim=2
  integer:: i,j,k,nn(ndim),isign

  do j=1,n
     do i=1,m
        k=i+m*(j-1)
        data(2*k-1)=real(s_fft(i,j))  !real part
        data(2*k)=imag(s_fft(i,j)) !imaginary part
     enddo
  enddo

  nn(1)=m
  nn(2)=n
  isign=1 !isign=1 corresponds to FFT in numerical-recipe-book, but corresponds the inverse FFT in FFTW
  call fourn(data,nn,ndim,isign)
  do i=1,m
     do j=1,n
        k=i+m*(j-1)
        s(i,j)=data(2*k-1)
     enddo
  enddo
s=s/(m*n)
end subroutine twod_inverse_fourier_transform_nr

SUBROUTINE fourn(data,nn,ndim,isign)
  use constants,only:p_
  implicit none
  INTEGER isign,ndim,nn(ndim)
  REAL(p_):: data(*)
  !Replaces data by its ndim -dimensional discrete Fourier transform, if isign is input as 1.
  !nn(1:ndim) is an integer array containing the lengths of each dimension (number of complex values), which MUST all be powers of 2. 
  !data is a real array of length twice the product of these lengths, in which the data are stored as in a multidimensional complex FORTRAN array
  !If isign is input as −1, data is replaced by its inverse transform times the product of the lengths of all dimensions. 
  !code from the numerical recipe book
  INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,&
       & ip3,k1,k2,n,nprev,nrem,ntot
  REAL(p_):: tempi,tempr
  real(p_):: theta,wi,wpi,wpr,wr,wtemp !Double precision for trigonometric re-currences.

  ntot=1
  do idim=1,ndim !Compute total number of complex values.
     ntot=ntot*nn(idim)
  enddo
  nprev=1
  do idim=1,ndim !Main loop over the dimensions.
     n=nn(idim)
     nrem=ntot/(n*nprev)
     ip1=2*nprev
     ip2=ip1*n
     ip3=ip2*nrem
     i2rev=1
     do  i2=1,ip2,ip1 !This is the bit-reversal section of the routine.
        if(i2.lt.i2rev)then
           do  i1=i2,i2+ip1-2,2
              do  i3=i1,ip3,ip2
                 i3rev=i2rev+i3-i2
                 tempr=data(i3)
                 tempi=data(i3+1)
                 data(i3)=data(i3rev)
                 data(i3+1)=data(i3rev+1)
                 data(i3rev)=tempr
                 data(i3rev+1)=tempi
              enddo
           enddo
        endif
        ibit=ip2/2
1       if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
           i2rev=i2rev-ibit
           ibit=ibit/2
           goto 1
        endif
        i2rev=i2rev+ibit
     enddo
     ifp1=ip1 !Here begins the Danielson-Lanczos section of the routine.
2    if(ifp1.lt.ip2)then
        ifp2=2*ifp1
        theta=isign*6.28318530717959d0/(ifp2/ip1) !Initialize for the trig. recur-rence.
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do  i3=1,ifp1,ip1
           do  i1=i3,i3+ip1-2,2
              do  i2=i1,ip3,ifp2
                 k1=i2 !Danielson-Lanczos formula:
                 k2=k1+ifp1
                 tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                 tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                 data(k2)=data(k1)-tempr
                 data(k2+1)=data(k1+1)-tempi
                 data(k1)=data(k1)+tempr
                 data(k1+1)=data(k1+1)+tempi
              enddo
           enddo
           wtemp=wr !Trigonometric recurrence.
           wr=wr*wpr-wi*wpi+wr
           wi=wi*wpr+wtemp*wpi+wi
        enddo
        ifp1=ifp2
        goto 2
     endif
     nprev=n*nprev
  enddo
  return
END SUBROUTINE fourn
end module fourn_module
