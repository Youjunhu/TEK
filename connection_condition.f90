module connection_condition
  use interpolate_module,only: linear_1d_interpolation
  implicit none

contains
  subroutine connection_condition_at_theta_cut(u) !interpolating data on theta=-pi plane to get values on theta=+pi plane
    !the two adjacent cells separated by theta=-pi or pi share the interface. However the grids on the interface of one cell is different from another,
    !specificly, the toroidal location of point indexed by (itor,jrad) on theta=+pi is different from the corresponding point indexed by (itor,jrad) on theta=-pi.
    !Therefore an interpolation over the toroidal angle is needed.
    use constants,only:p_
    use magnetic_coordinates,only:tor_1d_array,tor_shift_mc,mpol,j_low2,mtor
    use math,only: shift_to_specified_toroidal_range
    real(p_),intent(inout):: u(:,:)
    real(p_):: phi_old(mtor+1),u_old(mtor+1)
    real(p_):: phi_new
    integer:: i,j,jeq,m,n

    m=size(u,1)
    n=size(u,2)
    do j=1,n
       jeq=j-1+j_low2
       do i=1,mtor+1
          phi_old(i)=tor_1d_array(i)+tor_shift_mc(1,jeq) !at theta=-pi
          phi_old(i)=phi_old(i)-tor_shift_mc(1,jeq) !substracting by tor_shift_mc(1,jeq), which is a constant for the i=1,m+1 loop, is to re-define phi. If we do not re-define phi_old and use shift_to_specified_toroidal_rangle subroutine to shift the phi_old to the specified range, then phi_old array is changing in the way like 0.8pi->end_value-->starting_value--0.8pi, which is non-monotonical, which is not compatible with the interpolating function which requires the numerical table to be monotonical. The phi of the points that are to be interpolated should also be redefined in this way, i.e., subtracted by tor_shift_mc(1,jeq)
          !call shift_to_specified_toroidal_range(phi_old(i)) !this shift is not necessary because phi_old array is already in the specified range after we re-define it.
          if(i.eq.m+1) then
             u_old(i)=u(1,j) !using toroidal periodic condition, u(m+1,j)=u(1,j).
          else
             u_old(i)=u(i,j)
          endif
       enddo

       do i=1,m
          phi_new=tor_1d_array(i)+tor_shift_mc(mpol,jeq) !at theta=+pi
          phi_new=phi_new-tor_shift_mc(1,jeq) !re-define the toroidal angle in the same way that the above numerical table is defined, so that the following interpolating can work correctly
          call shift_to_specified_toroidal_range(phi_new) !this call is necessary
          call linear_1d_interpolation(mtor+1,phi_old,u_old,phi_new,u(i,j))  
       enddo
    enddo
  end subroutine connection_condition_at_theta_cut

  subroutine connection_condition_at_theta_cut3(u) !interpolating data on theta=+pi-dtheta2 plane to get data on grids defined on the plane of theta=-pi-dtheta2. The grids on theta=-pi-dtheta2 are determined by following the magnetic field line starting from the grids defined on theta=-pi plane to the plane of theta=-pi-dtheta2. !see the comments in subroutine "connection_condition_at_theta_cut"
    use constants,only:p_
    use constants,only: twopi
    use magnetic_coordinates,only:tor_1d_array,tor_shift_mc,mpol,j_low2,mtor
    use magnetic_coordinates,only: tor_shift_mc_left_bdry_minus_one
    use domain_decomposition, only : multi_eq_cells
    use math,only: shift_to_specified_toroidal_range
    implicit none
    real(p_),intent(inout):: u(:,:)
    real(p_):: phi_old(mtor+1),u_old(mtor+1)
    real(p_):: phi_at_midplane,phi_new
    integer:: i,j,jeq,m,n

    m=size(u,1)
    n=size(u,2)


    do j=1,n
       jeq=j-1+j_low2
       do i=1,mtor+1
          phi_old(i)=tor_1d_array(i)+tor_shift_mc(mpol-multi_eq_cells,jeq)
          phi_old(i)=phi_old(i)-tor_shift_mc(mpol-multi_eq_cells,jeq) !re-define phi_old by shifting phi_old by a contant, the reason is explained at the place labeled by mark2017-11-7-1 in this file
          if(i.eq.m+1) then
             u_old(i)=u(1,j) !using toroidal periodic condition, u(m+1,j)=u(1,j).
          else
             u_old(i)=u(i,j)
          endif
       enddo

       do i=1,m
          phi_at_midplane=tor_1d_array(i)+tor_shift_mc(1,jeq) !phi angle of point (i,j) on theta=-pi plane
          !phi_new=phi_at_midplane+(tor_shift_mc(mpol-multi_eq_cells,jeq)-tor_shift_mc(mpol,jeq)) !the phi angle of field-lines when it gos backward to theta=-pi-theata_interval
          phi_new=tor_1d_array(i)+tor_shift_mc_left_bdry_minus_one(jeq)  !the phi angle of field-lines when it gos backward to theta=-pi-theata_interval
          phi_new=phi_new-tor_shift_mc(mpol-multi_eq_cells,jeq) !shift by a constant, so it is consistent with the definition of phi_old
          call shift_to_specified_toroidal_range(phi_new)
          call linear_1d_interpolation(mtor+1,phi_old,u_old,phi_new,u(i,j))  
       enddo
    enddo
  end subroutine connection_condition_at_theta_cut3

  subroutine connection_condition_at_theta_cut0(u, ipol) 
    !Interpolate field at ipol poloidal plane. Handle 5 planes near theta=+-pi, namely, gclr=0,1,2, mpol2-1, mpol2-2. This is a generalized version of the above subroutines
    use constants, only : p_
    use magnetic_coordinates, only : tor_1d_array, tor_shift_mc, mpol, j_low2, mtor, mpol2
    use domain_decomposition, only : multi_eq_cells
    use math, only: shift_to_specified_toroidal_range
    integer, intent(in) :: ipol !uses zero-based index of the perturbation field, i.e., GCLR
    real(p_), intent(inout) :: u(:,:) !field values at grid of ipol plane, shape: (toroidal, radial)
    real(p_) :: phi_old(mtor+1), u_old(mtor+1)
    real(p_) :: phi_new, tor_shift1, tor_shift2
    integer :: i, j, jeq, m, n

    m=size(u,1)
    n=size(u,2)
    do j=1,n
       u_old(1:m)=u(:,j)
       u_old(m+1)=u_old(1) !toroidal periodic condition

       jeq=j-1+j_low2
       tor_shift1=tor_shift_mc(ipol*multi_eq_cells+1, jeq)
       phi_old(:)=tor_1d_array(:) + tor_shift1 - tor_shift1 !i.e., cancel the shift, why I did this way is explained above

       if((ipol==0) .or. (ipol==1) .or. (ipol==2)) then
          tor_shift2=tor_shift_mc(mpol,jeq)+(tor_shift_mc(ipol*multi_eq_cells+1,jeq) - tor_shift_mc(1,jeq))
       elseif(ipol==mpol2-1) then
          tor_shift2=tor_shift_mc(1,jeq) - (tor_shift_mc(mpol,jeq) - tor_shift_mc(mpol-1*multi_eq_cells,jeq))
       elseif(ipol==mpol2-2) then
          tor_shift2=tor_shift_mc(1,jeq) - (tor_shift_mc(mpol,jeq) - tor_shift_mc(mpol-2*multi_eq_cells,jeq))
       else
          stop "******toroidal shift interpolation at this poloidal location is not implemented****"
       endif

       do i=1,m
          phi_new=tor_1d_array(i)+ tor_shift2 -tor_shift1 !also cancel the shift1 to be consistent with phi_old defined above
          call shift_to_specified_toroidal_range(phi_new) 
          call linear_1d_interpolation(m+1,phi_old,u_old,phi_new,u(i,j))  
       enddo
    enddo
  end subroutine connection_condition_at_theta_cut0
  
  subroutine connection_condition_at_theta_cut_for_deposition(u) !interpolating data on theta=+pi plane to get values on theta=-pi plane, used in deposition
    use constants,only:p_
    use constants,only: twopi
    use magnetic_coordinates,only:m=>mtor, n=>nflux2, tor_1d_array,tor_shift_mc,mpol,j_low2
    use math,only: shift_to_specified_toroidal_range
    implicit none
    real(p_),intent(inout):: u(m,n)
    real(p_):: phi_old(m+1),u_old(m+1)
    real(p_):: phi_new
    integer:: i,j,jeq

    do j=1,n
       u_old(1:m) = u(:,j)
       u_old(m+1) = u(1,j) !toroidal periodic condition
       jeq=j-1+j_low2
       do i=1,m+1
          phi_old(i)=tor_1d_array(i)+tor_shift_mc(mpol,jeq) !at theta cut: theta=+pi
          phi_old(i)=phi_old(i)     -tor_shift_mc(mpol,jeq) !re-define phi_old, the reason is explained in "connection_condition_at_theta_cut"
       enddo

       do i=1,m
          phi_new=tor_1d_array(i)+tor_shift_mc(1,jeq) !at theta cut: theta=-pi.
          phi_new=phi_new        -tor_shift_mc(mpol,jeq) !The subtraction is to re-define phi,the reason is explained above
          call shift_to_specified_toroidal_range(phi_new)
          call linear_1d_interpolation(m+1,phi_old,u_old,phi_new,u(i,j))  
       enddo
    enddo
  end subroutine connection_condition_at_theta_cut_for_deposition

  subroutine components_transformation_at_theta_cut(ex,ey) !transform the components at theta=-pi to theta=pi
    use constants,only:p_
    use magnetic_coordinates,only: mpol,nflux2,j_low2,  &
     & grad_psi, grad_alpha,grad_psi_dot_grad_alpha, &
         & grad_psi_r,grad_psi_z, grad_alpha_r,grad_alpha_z,grad_alpha_phi
    implicit none
    real(p_),intent(inout):: ex(:,:),ey(:,:) !assumed-shape array
    real(p_)::source1(size(ex,1),size(ex,2)),source2(size(ex,1),size(ex,2)) !automatic array
    real(p_)::det_x(size(ex,1),size(ex,2)),det_y(size(ex,1),size(ex,2))
    real(p_)::det_c,a11,a12,a21,a22
    real(p_)::grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
    real(p_)::grad_alpha_dot_grad_alpha_positive,grad_psi_dot_grad_alpha_positive !_positive indicates value at theta=+pi, i.e., i=mpol
    integer:: j,jeq

    do j=1,nflux2
       jeq=j-1+j_low2
       grad_psi_val=grad_psi(1,jeq)
       grad_alpha_val=grad_alpha(1,jeq)
       grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha(1,jeq)

       source1(:,j)= ex(:,j)*grad_psi_val**2+ey(:,j)*grad_psi_dot_grad_alpha_val
       source2(:,j)= ex(:,j)*grad_psi_dot_grad_alpha_val+ey(:,j)*grad_alpha_val**2

       grad_psi_dot_grad_alpha_positive=grad_psi_r(1,jeq)*grad_alpha_r(mpol,jeq)+&
            & grad_psi_z(1,jeq)*grad_alpha_z(mpol,jeq)

       grad_alpha_dot_grad_alpha_positive=grad_alpha_r(1,jeq)*grad_alpha_r(mpol,jeq)&
            & +grad_alpha_z(1,jeq)*grad_alpha_z(mpol,jeq)&
            & +grad_alpha_phi(1,jeq)*grad_alpha_phi(mpol,jeq)

       a11=grad_psi_val**2
       a12=grad_psi_dot_grad_alpha_positive
       a21=grad_psi_dot_grad_alpha_val
       a22=grad_alpha_dot_grad_alpha_positive

       det_c=a11*a22-a12*a21 !determinant of the coefficient matrix
       det_x(:,j)=source1(:,j)*a22-source2(:,j)*a12
       det_y(:,j)=source2(:,j)*a11-source1(:,j)*a21
       ex(:,j)=det_x(:,j)/det_c !Cramer's rule to solve linear equation system
       ey(:,j)=det_y(:,j)/det_c !Cramer's rule to solve linear equation system
    enddo
  end subroutine components_transformation_at_theta_cut


  subroutine components_transformation_at_theta_cut_scalar_version(jeq, cxx,cxy,cyx,cyy) !transform the components at theta=-pi to theta=pi
    !(Ex, Ey) at theta=+pi is a function of (Ex, Ey) at theta=-pi
    !This routine returns the coefficient before (Ex, Ey) at theta=-pi in the function
    use constants,only:p_
    use magnetic_coordinates,only: mpol, &
      & grad_psi, grad_alpha,grad_psi_dot_grad_alpha, &
         & grad_psi_r,grad_psi_z, grad_alpha_r,grad_alpha_z,grad_alpha_phi
    implicit none
    integer,intent(in):: jeq
    real(p_),intent(out):: cxx,cxy,cyx,cyy
    real(p_)::det_c,a11,a12,a21,a22
    real(p_)::grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
    real(p_)::grad_alpha_dot_grad_alpha_positive,grad_psi_dot_grad_alpha_positive !_positive indicates value at theta=+pi, i.e., i=mpol

    grad_psi_val=grad_psi(1,jeq)
    grad_alpha_val=grad_alpha(1,jeq)
    grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha(1,jeq)

    grad_psi_dot_grad_alpha_positive=grad_psi_r(1,jeq)*grad_alpha_r(mpol,jeq)+&
         & grad_psi_z(1,jeq)*grad_alpha_z(mpol,jeq)

    grad_alpha_dot_grad_alpha_positive=grad_alpha_r(1,jeq)*grad_alpha_r(mpol,jeq)&
         & +grad_alpha_z(1,jeq)*grad_alpha_z(mpol,jeq)&
         & +grad_alpha_phi(1,jeq)*grad_alpha_phi(mpol,jeq)

    a11=grad_psi_val**2
    a12=grad_psi_dot_grad_alpha_positive
    a21=grad_psi_dot_grad_alpha_val
    a22=grad_alpha_dot_grad_alpha_positive

    det_c=a11*a22-a12*a21 !determinant of the coefficient matrix

    cxx=1._p_ !the coefficient before Ex_at_theta=pi in the function Ex_positive
    cxy=(a22*a21-a12*grad_alpha_val**2)/det_c !the coefficient before Ey_at_theta=pi  in the function Ex_positive
    cyx=0._p_ !the coefficient before Ex_at_theta=pi in the function Ey_positive
    cyy=(a11*grad_alpha_val**2-a21**2)/det_c !the coefficient before Ey_at_theta=pi in the function Ey_positive
  end subroutine components_transformation_at_theta_cut_scalar_version


  subroutine components_transformation_at_theta_cut2(ex,ey) 
    use constants,only:p_
    use magnetic_coordinates,only: mpol,j_low2,nflux2,dtheta, &
      & grad_psi, grad_alpha,grad_psi_dot_grad_alpha, &
      & grad_psi_r,grad_psi_z, grad_alpha_r,grad_alpha_z, &
    & grad_alpha_phi, grad_alpha_r_left_bdry_minus_one,grad_alpha_z_left_bdry_minus_one
    use domain_decomposition,only: dtheta2
    implicit none
    real(p_),intent(inout):: ex(:,:),ey(:,:) !assumed-shape array
    real(p_)::source1(size(ex,1),size(ex,2)),source2(size(ex,1),size(ex,2)) !automatic array
    real(p_)::det_x(size(ex,1),size(ex,2)),det_y(size(ex,1),size(ex,2))
    real(p_)::det_c,a11,a12,a21,a22
    real(p_)::gx0,grad_alpha1,grad_alpha2
    real(p_)::grad_psi_dot_grad_alpha1,grad_psi_dot_grad_alpha2
    real(p_)::grad_alpha1_dot_grad_alpha2
    integer:: j,jeq, ipol

    !ipol=mpol-1 !wrong, a bug found
    ipol=mpol-NINT(dtheta2/dtheta)
    do j=1,nflux2
       jeq=j-1+j_low2
       gx0=grad_psi(ipol,jeq)
       grad_alpha1=grad_alpha(ipol,jeq)
       grad_psi_dot_grad_alpha1=grad_psi_dot_grad_alpha(ipol,jeq)

       source1(:,j)= ex(:,j)*gx0**2+ey(:,j)*grad_psi_dot_grad_alpha1
       source2(:,j)= ex(:,j)*grad_psi_dot_grad_alpha1+ey(:,j)*grad_alpha1**2

       grad_psi_dot_grad_alpha2=grad_psi_r(ipol,jeq)*grad_alpha_r_left_bdry_minus_one(jeq) &
            & +grad_psi_z(ipol,jeq)*grad_alpha_z_left_bdry_minus_one(jeq)

       grad_alpha1_dot_grad_alpha2=grad_alpha_r(ipol,jeq)*grad_alpha_r_left_bdry_minus_one(jeq)&
            &+grad_alpha_z(ipol,jeq)*grad_alpha_z_left_bdry_minus_one(jeq)&
            & +grad_alpha_phi(ipol,jeq)**2

       a11=gx0**2
       a12=grad_psi_dot_grad_alpha2
       a21=grad_psi_dot_grad_alpha1
       a22=grad_alpha1_dot_grad_alpha2

       det_c=a11*a22-a12*a21 !determinant of the coefficient matrix
       !       write(*,*) 'det_c=',det_c
       det_x(:,j)=source1(:,j)*a22-source2(:,j)*a12
       det_y(:,j)=source2(:,j)*a11-source1(:,j)*a21
       ex(:,j)=det_x(:,j)/det_c !Cramer's rule to solve linear equation system
       ey(:,j)=det_y(:,j)/det_c !Cramer's rule to solve linear equation system
    enddo
  end subroutine components_transformation_at_theta_cut2


  subroutine components_transformation_at_theta_cut2_scalar_version(jeq, cxx,cxy,cyx,cyy)
    use constants,only:p_
    use magnetic_coordinates,only: mpol,dtheta, &
         & grad_psi, grad_alpha,grad_psi_dot_grad_alpha, grad_psi_r,&
         & grad_psi_z, grad_alpha_r,grad_alpha_z,grad_alpha_phi,&
         & grad_alpha_r_left_bdry_minus_one,grad_alpha_z_left_bdry_minus_one
    use domain_decomposition,only: dtheta2
    implicit none
    integer,intent(in):: jeq
    real(p_),intent(out):: cxx,cxy,cyx,cyy
    real(p_)::det_c,a11,a12,a21,a22
    real(p_)::gx0,grad_alpha1,grad_alpha2
    real(p_)::grad_psi_dot_grad_alpha1,grad_psi_dot_grad_alpha2
    real(p_)::grad_alpha1_dot_grad_alpha2
    integer:: j,ipol

    ipol=mpol-NINT(dtheta2/dtheta)

    gx0=grad_psi(ipol,jeq)
    grad_alpha1=grad_alpha(ipol,jeq)
    grad_psi_dot_grad_alpha1=grad_psi_dot_grad_alpha(ipol,jeq)

    grad_psi_dot_grad_alpha2=grad_psi_r(ipol,jeq)*grad_alpha_r_left_bdry_minus_one(jeq) &
         & +grad_psi_z(ipol,jeq)*grad_alpha_z_left_bdry_minus_one(jeq)

    grad_alpha1_dot_grad_alpha2=grad_alpha_r(ipol,jeq)*grad_alpha_r_left_bdry_minus_one(jeq)&
         &+grad_alpha_z(ipol,jeq)*grad_alpha_z_left_bdry_minus_one(jeq)&
         & +grad_alpha_phi(ipol,jeq)**2

    a11=gx0**2
    a12=grad_psi_dot_grad_alpha2
    a21=grad_psi_dot_grad_alpha1
    a22=grad_alpha1_dot_grad_alpha2

    det_c=a11*a22-a12*a21 !determinant of the coefficient matrix
    cxx=1.0_p_ !the coefficient before Ex_1 in the function Ex_2
    cxy=(a22*a21-a12*grad_alpha1**2)/det_c !the coefficient before Ey_1  in the function Ex_2
    !cyx=(a11*a21-a21*a11)/det_c !the coefficient before Ex_1 in the function Ey_2
    cyx=0.0_p_ !the coefficient before Ex_1 in the function Ey_2
    cyy=(a11*grad_alpha1**2-a21**2)/det_c !the coefficient before Ey_1 in the function Ey_2

  end subroutine components_transformation_at_theta_cut2_scalar_version

  subroutine check_components_transformation(ex,ey,pol_location)
    use constants,only:p_
    use magnetic_coordinates,only: mpol,j_low2,nflux2,dtheta, &
    & grad_psi, grad_alpha,grad_psi_dot_grad_alpha, &
         & grad_psi_r,grad_psi_z, grad_alpha_r,grad_alpha_z,grad_alpha_phi,&
         & grad_alpha_r_left_bdry_minus_one,grad_alpha_z_left_bdry_minus_one
    use domain_decomposition,only: dtheta2
    implicit none
    real(p_),intent(in):: ex(:,:),ey(:,:) !assumed-shape array
    integer,intent(in):: pol_location
    real(p_)::a(size(ex,1),size(ex,2)),b(size(ex,1),size(ex,2)),c(size(ex,1),size(ex,2)) !automatic array    
    real(p_)::er(size(ex,1),size(ex,2)),ez(size(ex,1),size(ex,2)),ephi(size(ex,1),size(ex,2))
    real(p_)::grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
    real(p_):: grad_psi_r_val,grad_psi_z_val,grad_alpha_r_val,grad_alpha_z_val
    real(p_):: grad_alpha_phi_val
    real(p_):: cos_alpha
    integer:: j,jeq

    do j=1,nflux2
       jeq=j-1+j_low2
       if(pol_location.eq.-1) then
          grad_psi_r_val=grad_psi_r(mpol-NINT(dtheta2/dtheta),jeq)
          grad_psi_z_val=grad_psi_z(mpol-NINT(dtheta2/dtheta),jeq)
          grad_psi_val=grad_psi(mpol-NINT(dtheta2/dtheta),jeq)

          grad_alpha_r_val=grad_alpha_r_left_bdry_minus_one(jeq)
          grad_alpha_z_val=grad_alpha_z_left_bdry_minus_one(jeq)
          grad_alpha_phi_val=grad_alpha_phi(mpol-NINT(dtheta2/dtheta),jeq)
          grad_alpha_val=sqrt(grad_alpha_r_val**2 +grad_alpha_z_val**2+grad_alpha_phi_val**2)
          grad_psi_dot_grad_alpha_val=grad_psi_r_val*grad_alpha_r_val+grad_psi_z_val*grad_alpha_z_val
       else
          grad_psi_r_val=grad_psi_r(pol_location,jeq)
          grad_psi_z_val=grad_psi_z(pol_location,jeq)
          grad_psi_val=grad_psi(pol_location,jeq)

          grad_alpha_r_val=grad_alpha_r(pol_location,jeq)
          grad_alpha_z_val=grad_alpha_z(pol_location,jeq)
          grad_alpha_phi_val=grad_alpha_phi(pol_location,jeq)
          grad_alpha_val=grad_alpha(pol_location,jeq)
          grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha(pol_location,jeq)
       endif
       cos_alpha=grad_psi_dot_grad_alpha_val/(grad_psi_val*grad_alpha_val)
       a(:,j)=ex(:,j)*grad_psi_val
       b(:,j)=ey(:,j)*grad_alpha_val

       c(:,j)=sqrt(a(:,j)**2+b(:,j)**2+2*a(:,j)*b(:,j)*cos_alpha) !using cosine theorem to calculate the magnitude

       er(:,j)=ex(:,j)*grad_psi_r_val+ey(:,j)*grad_alpha_r_val !componet in cylindreical coordinates
       ez(:,j)=ex(:,j)*grad_psi_z_val+ey(:,j)*grad_alpha_z_val !componet in cylindreical coordinates
       ephi(:,j)=ey(:,j)*grad_alpha_phi_val !componet in cylindreical coordinates
    enddo

    write(*,*) 'pol_location=',pol_location,'magnitude=', c(15,15),'R comp=',er(15,15),'Z comp=',ez(15,15),'Phi comp=',ephi(15,15) !the values across the theta cut should be equal to each other
  end subroutine check_components_transformation

end module connection_condition
