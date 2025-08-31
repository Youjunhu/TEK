module deposit_fk_module
  implicit none
contains
  subroutine deposit_fk(nmarker_i,active_i,radcor_i,theta_i,alpha_i,phi_i,w_i) !given (r,v,w) of markers, do the deposition to get perpendicular currents and number density on spatial grids
    !periodic toroidal boundary condition and the poloidal boundary condtion (connection condition along the field line) are taken into account
    !    use mpi
    use constants,only:p_
    use normalizing,only: nu
    use constants,only: one,two
    use fk_module,only: ni0, my_den_i_left, my_den_i_right !as output
    use magnetic_coordinates,only:m=>mtor,n=>nrad
    use magnetic_coordinates,only:xgrid,zgrid,ygrid,dtor,dradcor,dtheta,jacobian
    !use magnetic_coordinates,only:phi_grid_left,phi_grid_right
    use domain_decomposition,only: dtheta2,theta_start,ipol_eq, multi_eq_cells


    integer,intent(in)::nmarker_i
    logical,intent(in):: active_i(nmarker_i)
    real(p_),intent(in)::radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i),phi_i(nmarker_i)
    real(p_),intent(in):: w_i(nmarker_i)
    real(p_):: coeff_theta_1,coeff_theta_2,coeff_alpha_1,coeff_alpha_2,coeff_radcor_1,coeff_radcor_2
    real(p_):: dv1, dv2
    integer:: k,i,j, jeq,i_plus_one,j_plus_one

    !before the deposition, set the array to zero
    my_den_i_left=0._p_
    my_den_i_right=0._p_

    !$omp parallel do private(coeff_theta_1,coeff_theta_2, coeff_alpha_1, coeff_alpha_2, &
    !$omp & coeff_radcor_1, coeff_radcor_2,i,j,i_plus_one,j_plus_one,dphi_left00,dphi_left01,dphi_left10,dphi_left11,&
    !$omp & dphi_right00, dphi_right01, dphi_right10,dphi_right11,kernel) !!&
!!!!$omp & reduction(+:my_jr_left,my_jr_right,my_jphi_left,my_jphi_right,my_jz_left,my_jz_right,my_den_i_left,my_den_i_right)
    do k=1,nmarker_i !for each marker, deposit it to the corresponding grids
       if(active_i(k).eqv..false.) cycle ! markers outside the computational region do not contribute density or current to any grids
       !determine the interpolating coeefficeint according to the location of a marker
       coeff_theta_1=(theta_i(k)-theta_start)/dtheta2
       coeff_theta_2=one-coeff_theta_1

       !call location(m,ygrid,alpha_i(k),i)
       i=floor((alpha_i(k)-ygrid(1))/dtor+1) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
       coeff_alpha_1=(alpha_i(k)-ygrid(i))/dtor
       coeff_alpha_2=one-coeff_alpha_1
       !if(myid.eq.2) write(*,*) 'alpha, i=',i, 'k=',k
       !call location(n, xgrid, radcor_i(k),j)
       j=floor((radcor_i(k)-xgrid(1))/dradcor+1)
       coeff_radcor_1= (radcor_i(k)-xgrid(j))/dradcor
       coeff_radcor_2=one-coeff_radcor_1

       i_plus_one=i+1
       if(i.eq.m) i_plus_one=1 !periodic condition in the toroidal direction
       j_plus_one=j+1
       if(j.eq.n) j_plus_one=j !all the active markers are actually within the boundary flux surface labeled by n. If j.eq.n, then the marker is exactly on the boundary flux surface, this code line is needed to avoid exeeding array bounds in case that the marker is exactly on the boundary flux surface or slightly outside of it due to numerical truncation errors
       !if(j.eq.n) j_plus_one=1 !periodic condition. is disabled later in this subroutine, now using fixed zero boundary condition along the radial direction. In the past, I used the periodic bounary condition because I want to use Fourier transform along the radial direction. However this periodic radial boundary condition is not reasonable. The reason is as follows: Here the radial direction is d/dpsi in (psi,theta,alpha) coordinates. This direction is a combinition of the toroidal and the usual radial direction and thus is an artificial direction introduced when using (psi,theta,alpha) coordinates, and does not has any physical reason to justify that a peroidic condition should be satisfied along this artificial direction.

       !density, only used in isotropic fluid electrons model or adiabatic electrons model to provide delta_ne, which is assumed to be equal to delta_ni
       my_den_i_left(i,j)=my_den_i_left(i,j) +w_i(k)*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
       my_den_i_right(i,j)=my_den_i_right(i,j)+w_i(k)*coeff_theta_1*coeff_alpha_2*coeff_radcor_2
       my_den_i_left(i_plus_one,j)= my_den_i_left(i_plus_one,j)+w_i(k)*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
       my_den_i_right(i_plus_one,j)=my_den_i_right(i_plus_one,j)+w_i(k)*coeff_theta_1*coeff_alpha_1*coeff_radcor_2
       my_den_i_left(i,j_plus_one)=my_den_i_left(i,j_plus_one)+w_i(k)*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
       my_den_i_right(i,j_plus_one)=my_den_i_right(i,j_plus_one)+w_i(k)*coeff_theta_1*coeff_alpha_2*coeff_radcor_1
       my_den_i_left(i_plus_one,j_plus_one)=my_den_i_left(i_plus_one,j_plus_one)+w_i(k)*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
       my_den_i_right(i_plus_one,j_plus_one)=my_den_i_right(i_plus_one,j_plus_one)+w_i(k)*coeff_theta_1*coeff_alpha_1*coeff_radcor_1
    enddo
    !$omp end parallel do

    do j=1,n  !divided by the space volume of a cell, to give the current denisty, note that this 'cell' is the cell defined by PIC (i.e., grid is the center of the cell). (while grids are the boundaries of the "mpi cell" for grouping and pushing particles)
       jeq=j
       dv1=abs(jacobian(ipol_eq,jeq))*dradcor*dtheta2*dtor !volume of the cell (the center of the cell is the grid)
       dv2=abs(jacobian(ipol_eq+multi_eq_cells,jeq))*dradcor*dtheta2*dtor !volume of the cell
       my_den_i_left(:,j)  = my_den_i_left (:,j)/dv1
       my_den_i_right(:,j) = my_den_i_right(:,j)/dv2
    enddo

    my_den_i_left  = my_den_i_left /nu
    my_den_i_right = my_den_i_right/nu
  end subroutine deposit_fk
end module deposit_fk_module
