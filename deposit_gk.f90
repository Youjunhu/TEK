module deposit_gk_module
contains
  subroutine deposit_gk(ns, touch_bdry_e, theta_e, vpar_e, w_e, x_ring,y_ring)
    !with markers' (v,x,w) given, do the deposition to get density and jpar on grids
    !use omp_lib
    use constants,only : p_, one
    use magnetic_coordinates,only:n=>nflux2,m=>mtor,xarray=>radcor_1d_array2,yarray=>tor_1d_array,dtor,dradcor
    use domain_decomposition,only: dtheta2,theta_start
    use gk_module,only: nm_gk, gk_flr, gyro_npt, charge_e, vn_e
    use normalizing, only : vu,  qu
    use perturbation_field_matrix,only: my_den_e_left, my_den_e_right, & !input and output
         & my_jpar_e_left, my_jpar_e_right !input and output
    implicit none
    integer,intent(in)  :: ns
    logical,intent(in)  :: touch_bdry_e(:)
    real(p_),intent(in) :: theta_e(:),  vpar_e(:), w_e(:), x_ring(:,:), y_ring(:,:)
    real(p_) :: cz1, cz2, cy1, cy2, cx1, cx2
    real(p_) :: kernel, kernel2
    integer  :: i,j,k, i_plus1, j_plus1, kr, npt0, nm
    
    nm= nm_gk(ns)
    if(gk_flr(ns) .eqv. .false.) then
       npt0=1
    else
       npt0=gyro_npt
    endif

    do k=1,nm
       if (touch_bdry_e(k) .eqv. .true.) cycle ! markers outside the computational region do not contribute anything to any grids
       kernel=(charge_e(ns)/qu)*w_e(k)/npt0 !for computing charge density
       kernel2=kernel*vpar_e(k)*(vn_e(ns)/vu) !for computing parallel current
       cz1=(theta_e(k)-theta_start)/dtheta2 !interpolating coefficient
       cz2=one-cz1
       do kr=1,npt0 !loop over the points on a gyro-ring
          i=floor((y_ring(kr,k)-yarray(1))/dtor+1) 
          cy1=(y_ring(kr,k)-yarray(i))/dtor
          cy2=one-cy1

          j=floor((x_ring(kr,k)-xarray(1))/dradcor+1)
          cx1= (x_ring(kr,k)-xarray(j))/dradcor
          cx2=one-cx1

          i_plus1=i+1
          if(i.eq.m) i_plus1=1 !periodic condition
          j_plus1=j+1
          if(j.eq.n) cycle !marker is out of radial computational region

          my_den_e_left(i,j)= my_den_e_left(i,j) +kernel*cz2*cy2*cx2
          my_den_e_right(i,j)=my_den_e_right(i,j)+kernel*cz1*cy2*cx2
          my_den_e_left(i_plus1,j)=my_den_e_left(i_plus1,j)  +kernel*cz2*cy1*cx2
          my_den_e_right(i_plus1,j)=my_den_e_right(i_plus1,j)+kernel*cz1*cy1*cx2
          my_den_e_left(i,j_plus1)= my_den_e_left(i,j_plus1) +kernel*cz2*cy2*cx1
          my_den_e_right(i,j_plus1)=my_den_e_right(i,j_plus1)+kernel*cz1*cy2*cx1
          my_den_e_left(i_plus1,j_plus1)=my_den_e_left(i_plus1,j_plus1)+kernel*cz2*cy1*cx1
          my_den_e_right(i_plus1,j_plus1)=my_den_e_right(i_plus1,j_plus1)+kernel*cz1*cy1*cx1

          my_jpar_e_left(i,j)= my_jpar_e_left(i,j) +kernel2*cz2*cy2*cx2
          my_jpar_e_right(i,j)=my_jpar_e_right(i,j)+kernel2*cz1*cy2*cx2
          my_jpar_e_left(i_plus1,j)=my_jpar_e_left(i_plus1,j)  +kernel2*cz2*cy1*cx2
          my_jpar_e_right(i_plus1,j)=my_jpar_e_right(i_plus1,j)+kernel2*cz1*cy1*cx2
          my_jpar_e_left(i,j_plus1)= my_jpar_e_left(i,j_plus1) +kernel2*cz2*cy2*cx1
          my_jpar_e_right(i,j_plus1)=my_jpar_e_right(i,j_plus1)+kernel2*cz1*cy2*cx1
          my_jpar_e_left(i_plus1,j_plus1)=my_jpar_e_left(i_plus1,j_plus1)+kernel2*cz2*cy1*cx1
          my_jpar_e_right(i_plus1,j_plus1)=my_jpar_e_right(i_plus1,j_plus1)+kernel2*cz1*cy1*cx1

       enddo
    enddo

  end subroutine deposit_gk

end module deposit_gk_module
