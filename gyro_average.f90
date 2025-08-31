module gyro_average_mod
contains
pure  subroutine gyro_average0(flr, x_ring, y_ring, z_ring, touch_bdry, phix, phix_av)
    use constants, only:  p_
    use gk_module, only: gyro_npt
    implicit none
    real(p_), intent(in) :: x_ring(:), y_ring(:), z_ring(:)
    logical, intent(in)  :: flr, touch_bdry
    real(p_), intent(in) :: phix(:,:,:)
    real(p_), intent(out) :: phix_av
    real(p_) :: phixp(gyro_npt)
    integer :: kr, npt

    if(touch_bdry .eqv. .true.) then
       phix_av = 0._p_
       return
    endif
    
    if (flr .eqv. .false.) then
       npt = 1
    else
       npt = gyro_npt
    endif

    do kr = 1, npt
       call field_at_particle0(x_ring(kr), y_ring(kr), z_ring(kr), touch_bdry, phix, phixp(kr))
    enddo
    phix_av =sum(phixp(1:npt))/npt

  end subroutine gyro_average0

  pure subroutine field_at_particle0(radcor, alpha, theta, touch_bdry, phix, phixp)
    use constants, only: p_, one, zero, twopi
    use magnetic_coordinates, only: mtor, nrad, ygrid, xgrid
    use domain_decomposition, only: dtheta2, theta_start
    implicit none
    real(p_), intent(in) :: radcor, theta, alpha, phix(:,:,:)
    logical, intent(in) :: touch_bdry
    real(p_), intent(out) :: phixp
    real(p_) :: c1, c2, tmp1, tmp2

    if(touch_bdry.eqv..true.) then !force on particles outside the computational region is set to zero
       phixp=zero
    else
       c1 = (theta-theta_start)/dtheta2
       c2 = one-c1
       tmp1 = interpolate2d(mtor,nrad, ygrid, xgrid, phix(:,:,1), alpha,radcor) 
       tmp2 = interpolate2d(mtor,nrad, ygrid, xgrid, phix(:,:,2), alpha,radcor) 
       phixp = tmp1*c2 + tmp2*c1
    endif
  end subroutine field_at_particle0


 pure subroutine gyro_average(ns, nm, x_ring, y_ring, z_ring, touch_bdry_gc, &
       & phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga)
    !output are used in computing drift and pushing weight.
    use constants, only:  p_
    use gk_module, only: gk_flr, gyro_npt
    implicit none
    integer, intent(in)  :: ns, nm
    real(p_), intent(in) :: x_ring(:,:), y_ring(:,:), z_ring(:,:)
    logical, intent(in)  :: touch_bdry_gc(nm)
    real(p_), intent(out) :: phix_ga(nm), phiy_ga(nm), phiz_ga(nm)
    real(p_), intent(out) :: ax_ga(:),ay_ga(:),az_ga(:), ahx_ga(:),ahy_ga(:),ahz_ga(:), ah_ga(:)
    real(p_), dimension(gyro_npt) :: phix0, phiy0, phiz0, ax0, ay0, az0, ahx0, ahy0, ahz0, ah0
    integer :: k, kr, npt

    if(gk_flr(ns) .eqv. .false.) then
       npt=1
    else
       npt=gyro_npt
    endif

    do k = 1, nm
       if( touch_bdry_gc(k) .eqv. .true.) cycle
       do kr = 1, npt
          call field_at_particle(x_ring(kr,k), y_ring(kr,k), z_ring(kr,k), touch_bdry_gc(k), & 
                              & phix0(kr), phiy0(kr), phiz0(kr), ax0(kr), ay0(kr), az0(kr), &
                              & ahx0(kr), ahy0(kr), ahz0(kr), ah0(kr))
       enddo
       phix_ga(k) =sum(phix0(1:npt))/npt
       phiy_ga(k) =sum(phiy0(1:npt))/npt
       phiz_ga(k) =sum(phiz0(1:npt))/npt
       ax_ga(k)   =sum(ax0(1:npt))/npt
       ay_ga(k)   =sum(ay0(1:npt))/npt
       az_ga(k)   =sum(az0(1:npt))/npt
       ahx_ga(k)  =sum(ahx0(1:npt))/npt
       ahy_ga(k)  =sum(ahy0(1:npt))/npt
       ahz_ga(k)  =sum(ahz0(1:npt))/npt
       ah_ga(k)   =sum(ah0(1:npt))/npt
    enddo

  end subroutine gyro_average


  pure subroutine field_at_particle(radcor, alpha, theta, touch_bdry, phix0,phiy0,phiz0,ax0,ay0,az0,ahx0,ahy0,ahz0,ah0)
    use constants, only: p_, one, zero, twopi
    use perturbation_field,only: phix, phiy, phiz, ax, ay, az, ahx, ahy, ahz, apara_h
    use magnetic_coordinates,only: mtor,nrad, ygrid, xgrid
    use domain_decomposition,only: dtheta2,theta_start
    implicit none
    real(p_), intent(in) :: radcor, theta, alpha
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
       ahx0=zero
       ahy0=zero
       ahz0=zero
       ah0=zero
    else
       coeff1=(theta-theta_start)/dtheta2
       coeff2=one-coeff1

       !  if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,phix(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,phix(:,:,2) ,alpha,radcor) 
       phix0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,phiy(:,:,1) ,alpha,radcor)  
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,phiy(:,:,2) ,alpha,radcor)  
       phiy0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,phiz(:,:,1) ,alpha,radcor)  
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,phiz(:,:,2) ,alpha,radcor)  
       phiz0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ax(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ax(:,:,2) ,alpha,radcor) 
       ax0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ay(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ay(:,:,2) ,alpha,radcor) 
       ay0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,az(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,az(:,:,2) ,alpha,radcor) 
       az0=tmp1*coeff2+tmp2*coeff1


       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ahx(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ahx(:,:,2) ,alpha,radcor) 
       ahx0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ahy(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ahy(:,:,2) ,alpha,radcor) 
       ahy0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ahz(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ahz(:,:,2) ,alpha,radcor) 
       ahz0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,apara_h(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,apara_h(:,:,2) ,alpha,radcor) 
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
    if(j>=nz .or. j<=1) then
       fval =0
       return
    endif

    slope=(psi(i_plus1,j)-psi(i,j))/dx
    t1=psi(i,j)+slope*(x-xarray(i))
    slope=(psi(i_plus1,j+1)-psi(i,j+1))/dx
    t2=psi(i,j+1)+slope*(x-xarray(i))
    slope=(t2-t1)/dz
    fval=t1+slope*(z-zarray(j))
  end function interpolate2d

end module gyro_average_mod
