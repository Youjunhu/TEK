module gyro_average_mod
contains
  subroutine gyro_average0(flr, z, x_ring, y_ring,  touch_bdry,phix, phix_av)
    use constants, only :  p_
    use gk_module,only: gyro_npt
    implicit none
    real(p_),intent(in) :: z, x_ring(:), y_ring(:)
    logical,intent(in)  :: flr, touch_bdry
    real(p_),intent(in) :: phix(:,:,:)
    real(p_),intent(out) :: phix_av
    real(p_) :: phixp(gyro_npt)
    integer :: kr, npt

    if(touch_bdry .eqv. .true.) then
       phix_av=0._p_
       return
    endif
    if (flr .eqv. .false.) then
       npt=1
    else
       npt=gyro_npt
    endif

    do kr=1,npt
       call field_at_particle0(x_ring(kr), z, y_ring(kr), touch_bdry, phix, phixp(kr))
    enddo
    phix_av =sum(phixp(1:npt))/npt

  end subroutine gyro_average0

  subroutine field_at_particle0(radcor,theta,alpha,touch_bdry, phix, phixp)
    use constants,only:p_
    use constants,only: one,zero,twopi
    use magnetic_coordinates,only: mtor,nflux2,tor_1d_array,radcor_1d_array2
    use domain_decomposition,only: dtheta2,theta_start
    implicit none
    real(p_),intent(in) :: radcor,theta,alpha, phix(:,:,:)
    logical,intent(in) :: touch_bdry
    real(p_),intent(out) :: phixp
    real(p_) :: coeff1,coeff2,tmp1,tmp2

    if(touch_bdry.eqv..true.) then !force on particles outside the computational region is set to zero
       phixp=zero
    else
       coeff1=(theta-theta_start)/dtheta2
       coeff2=one-coeff1
       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phix(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phix(:,:,2) ,alpha,radcor) 
       phixp=tmp1*coeff2+tmp2*coeff1
    endif
  end subroutine field_at_particle0


  subroutine gyro_average(ns,nm, radcor_e, theta_e, alpha_e, x_ring, y_ring,  touch_bdry_e, &
       & phix_e, phiy_e, phiz_e, ax_e, ay_e, az_e, ahx_e, ahy_e, ahz_e, ah_e)
    !output are used in computing drift and pusing weight.
    use constants, only :  p_
    use gk_module,only: gk_flr, gyro_npt
    implicit none
    integer,intent(in)  :: ns,nm
    real(p_),intent(in) :: radcor_e(:), theta_e(:), alpha_e(:)
    real(p_),intent(in) :: x_ring(:,:), y_ring(:,:)
    logical,intent(in)  :: touch_bdry_e(:)
    real(p_),intent(out) :: phix_e(:), phiy_e(:), phiz_e(:)
    real(p_),intent(out) :: ax_e(:),ay_e(:),az_e(:), ahx_e(:),ahy_e(:),ahz_e(:), ah_e(:)
    real(p_),dimension(gyro_npt) :: phix0, phiy0, phiz0, ax0, ay0, az0, ahx0, ahy0, ahz0, ah0
    integer :: k, kr, npt

    if(gk_flr(ns) .eqv. .false.) then
       npt=1
    else
       npt=gyro_npt
    endif

    do k=1,nm
       if( touch_bdry_e(k).eqv..true.) cycle
       do kr=1,npt
          call field_at_particle(x_ring(kr,k),theta_e(k),y_ring(kr,k),touch_bdry_e(k), & 
               phix0(kr),phiy0(kr),phiz0(kr), ax0(kr),ay0(kr),az0(kr),ahx0(kr),ahy0(kr),ahz0(kr),ah0(kr))
       enddo
       phix_e(k) =sum(phix0(1:npt))/npt
       phiy_e(k) =sum(phiy0(1:npt))/npt
       phiz_e(k) =sum(phiz0(1:npt))/npt
       ax_e(k)   =sum(ax0(1:npt))/npt
       ay_e(k)   =sum(ay0(1:npt))/npt
       az_e(k)   =sum(az0(1:npt))/npt
       ahx_e(k)  =sum(ahx0(1:npt))/npt
       ahy_e(k)  =sum(ahy0(1:npt))/npt
       ahz_e(k)  =sum(ahz0(1:npt))/npt
       ah_e(k)   =sum(ah0(1:npt))/npt
    enddo

  end subroutine gyro_average


  subroutine field_at_particle(radcor,theta,alpha,touch_bdry, phix0,phiy0,phiz0,ax0,ay0,az0,ahx0,ahy0,ahz0,ah0)
    use constants,only:p_
    use constants,only: one,zero,twopi
    use perturbation_field_matrix,only: phix, phiy, phiz, ax,ay,az, ahx, ahy, ahz, apara_h
    use magnetic_coordinates,only: mtor,nflux2,tor_1d_array,radcor_1d_array2
    use domain_decomposition,only: dtheta2,theta_start
    implicit none
    real(p_),intent(in)::radcor,theta,alpha
    logical,intent(in):: touch_bdry
    real(p_),intent(out) :: phix0, phiy0, phiz0, ax0, ay0, az0, ahx0, ahy0, ahz0, ah0
    real(p_):: coeff1,coeff2,tmp1,tmp2

    if (touch_bdry .eqv. .true.) then !force on particles outside the computational region is set to zero
       phix0=zero
       phiy0=zero
       phiz0=zero
       ax0=zero
       ay0=zero
       az0=zero
       ahx=zero
       ahy=zero
       ahz=zero
       ah0=zero
    else
       coeff1=(theta-theta_start)/dtheta2
       coeff2=one-coeff1

       !  if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phix(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phix(:,:,2) ,alpha,radcor) 
       phix0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phiy(:,:,1) ,alpha,radcor)  
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phiy(:,:,2) ,alpha,radcor)  
       phiy0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phiz(:,:,1) ,alpha,radcor)  
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,phiz(:,:,2) ,alpha,radcor)  
       phiz0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ax(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ax(:,:,2) ,alpha,radcor) 
       ax0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ay(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ay(:,:,2) ,alpha,radcor) 
       ay0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,az(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,az(:,:,2) ,alpha,radcor) 
       az0=tmp1*coeff2+tmp2*coeff1


       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ahx(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ahx(:,:,2) ,alpha,radcor) 
       ahx0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ahy(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ahy(:,:,2) ,alpha,radcor) 
       ahy0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ahz(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,ahz(:,:,2) ,alpha,radcor) 
       ahz0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,apara_h(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nflux2,tor_1d_array,radcor_1d_array2,apara_h(:,:,2) ,alpha,radcor) 
       ah0=tmp1*coeff2+tmp2*coeff1

    endif
  end subroutine field_at_particle


  pure function interpolate2d(nx,nz,xarray,zarray,psi,x,z) result(fval)
    use constants,only : p_, one
    implicit none
    integer,intent(in) :: nx,nz
    real(p_),intent(in) :: xarray(nx), zarray(nz), psi(nx,nz)
    real(p_),intent(in) :: x,z
    real(p_) :: fval
    real(p_) :: dx,dz,t1,t2,slope
    integer :: i,j ,ii,jj, i_plus1

    dx=xarray(2)-xarray(1)
    i=floor(one+(x-xarray(1))/dx)
    i_plus1 = i+1
    if(i.eq.nx) i_plus1 = 1 !periodic condition

    dz=zarray(2)-zarray(1)
    j=floor(one+(z-zarray(1))/dz)
    if(j.eq.nz) j=nz-1

    slope=(psi(i_plus1,j)-psi(i,j))/dx
    t1=psi(i,j)+slope*(x-xarray(i))
    slope=(psi(i_plus1,j+1)-psi(i,j+1))/dx
    t2=psi(i,j+1)+slope*(x-xarray(i))
    slope=(t2-t1)/dz
    fval=t1+slope*(z-zarray(j))
  end function interpolate2d

end module gyro_average_mod
