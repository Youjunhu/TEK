module deposit_gk_module
contains
  
  subroutine deposit_gk(ns, touch_bdry_gc, vpar_gk, w_gk, x_ring, y_ring, z_ring, &
       & density_left, density_right, jpar_left, jpar_right)
    !with markers' (v,x,w) given, do the deposition to get density and jpar on grids
    use constants, only: p_, one
    use magnetic_coordinates, only: n=>nrad,m=>mtor,xarray=>xgrid,yarray=>ygrid,dtor,dradcor
    use domain_decomposition, only: dtheta2, theta_start
    use gk_module, only: nm_gk, gk_flr, gyro_npt, charge_gk, vn_gk
    use normalizing, only: vu,  qu
    implicit none
    integer,intent(in)  :: ns
    logical,intent(in)  :: touch_bdry_gc(:)
    real(p_),intent(in) :: vpar_gk(:), w_gk(:), x_ring(:,:), y_ring(:,:), z_ring(:,:)
    real(p_), intent(inout) :: density_left(:,:), density_right(:,:)
    real(p_), intent(inout) :: jpar_left(:,:), jpar_right(:,:)
    real(p_) :: cz1, cz2, cy1, cy2, cx1, cx2, kernel, kernel2
    integer  :: i,j,k, iz, i_plus1, j_plus1, kr, npt0, nm

    nm= nm_gk(ns)
    if(gk_flr(ns) .eqv. .false.) then
       npt0 = 1
    else
       npt0 = gyro_npt
    endif

    do k=1,nm
       if (touch_bdry_gc(k) .eqv. .true.) cycle ! markers outside contribute nothing
       kernel = (charge_gk(ns)/qu)*w_gk(k)/npt0 !for computing charge density
       kernel2 = kernel*vpar_gk(k)*(vn_gk(ns)/vu) !for computing parallel current

       do kr = 1,npt0 !loop over the points on a gyro-ring
          i=floor((y_ring(kr,k)-yarray(1))/dtor+1) 
          cy1=(y_ring(kr,k)-yarray(i))/dtor
          cy2=one-cy1

          j=floor((x_ring(kr,k)-xarray(1))/dradcor+1)
          cx1= (x_ring(kr,k)-xarray(j))/dradcor
          cx2=one-cx1

          cz1=(z_ring(kr,k)-theta_start)/dtheta2
          cz2=one-cz1

          i_plus1=i+1
          if(i.eq.m) i_plus1=1 !periodic condition
          j_plus1=j+1
          if(j.eq.n) cycle !marker is out of radial computational region

          density_left(i,j) = density_left(i,j) +kernel*cz2*cy2*cx2
          density_left(i_plus1,j) = density_left(i_plus1,j)  +kernel*cz2*cy1*cx2
          density_left(i,j_plus1) = density_left(i,j_plus1) +kernel*cz2*cy2*cx1
          density_left(i_plus1,j_plus1) = density_left(i_plus1,j_plus1)+kernel*cz2*cy1*cx1

          density_right(i,j) = density_right(i,j)+kernel*cz1*cy2*cx2
          density_right(i_plus1,j) = density_right(i_plus1,j)+kernel*cz1*cy1*cx2
          density_right(i,j_plus1) = density_right(i,j_plus1)+kernel*cz1*cy2*cx1
          density_right(i_plus1,j_plus1) = density_right(i_plus1,j_plus1)+kernel*cz1*cy1*cx1

          jpar_left(i,j) = jpar_left(i,j) + kernel2*cz2*cy2*cx2
          jpar_left(i_plus1,j) = jpar_left(i_plus1,j) + kernel2*cz2*cy1*cx2
          jpar_left(i,j_plus1) = jpar_left(i,j_plus1) + kernel2*cz2*cy2*cx1
          jpar_left(i_plus1,j_plus1) = jpar_left(i_plus1,j_plus1)+kernel2*cz2*cy1*cx1

          jpar_right(i,j) = jpar_right(i,j) + kernel2*cz1*cy2*cx2
          jpar_right(i_plus1,j) = jpar_right(i_plus1,j)+ kernel2*cz1*cy1*cx2
          jpar_right(i,j_plus1) = jpar_right(i,j_plus1) + kernel2*cz1*cy2*cx1
          jpar_right(i_plus1,j_plus1) = jpar_right(i_plus1,j_plus1)+kernel2*cz1*cy1*cx1          
       enddo

!!$       do kr = 1, npt0 !nearest neighbour interpolation, turns out to give more noisy results
!!$          i = nint((y_ring(kr,k)-yarray(1))/dtor)+1
!!$          j = nint((x_ring(kr,k)-xarray(1))/dradcor)+1
!!$          iz = nint((z_ring(kr,k)-theta_start)/dtheta2)
!!$
!!$          if(iz==0) then
!!$             density_left(i,j) = density_left(i,j) + kernel
!!$             jpar_left(i,j) = jpar_left(i,j) + kernel2
!!$          else
!!$             density_right(i,j) = density_right(i,j) + kernel
!!$             jpar_right(i,j) = jpar_right(i,j) + kernel2
!!$          endif
!!$       enddo

    enddo
  end subroutine deposit_gk

end module deposit_gk_module
