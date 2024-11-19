module derivatives_in_xyz
  use constants,only : p_, two
  implicit none

contains

  subroutine radial_derivative(field,field_x) !calculate the derivative of field with respect to psi in field-line-following-coordinates (psi,theta,alpha)
    use magnetic_coordinates,only: dradcor
    real(p_),intent(in)::field(:,:)
    real(p_),intent(out)::field_x(:,:)
    integer:: m,n,i,j

    m=size(field,1)
    n=size(field,2)
    do i=1,m
       do j=2,n-1
          field_x(i,j)= (field(i,j+1)-field(i,j-1))/(two*dradcor)
       enddo
!!$       field_x(i,1)=(field(i,2)-0._p_)/(two*dradcor) !zero boundary condition for the field is used
!!$       field_x(i,n)=(0._p_-field(i,n-1))/(two*dradcor)!zero boundary condition for the field is used
       field_x(i,1)=two*field_x(i,2)-field_x(i,3) !linear interpolation to obtain the derivative at the boundary point
       field_x(i,n)=two*field_x(i,n-1)-field_x(i,n-2) !linear interpolation
    enddo

  end subroutine radial_derivative

  subroutine toroidal_derivative(field,field_y) !calculate the derivative of pper_e with respect to y:
    use magnetic_coordinates,only: dtor
    real(p_),intent(in)::field(:,:)
    real(p_),intent(out)::field_y(:,:)
    integer:: m,n,i,j,i_left,i_right

    m=size(field,1)
    n=size(field,2)
    do i=1,m   
       do j=1,n
          i_left=i-1
          if(i.eq.1) i_left=m !periodic boundary condition
          i_right=i+1
          if(i.eq.m) i_right=1 !periodic boundary condition
          field_y(i,j)=(field(i_right,j)-field(i_left,j))/(two*dtor)
       enddo
    enddo

  end subroutine toroidal_derivative

  subroutine theta_derivative(a, a_theta) !partial derivative with respect to theta (with psi and alpha fixed), i.e., along the magnetic field line
    use constants,only : two
    use domain_decomposition,only : dtheta2,myid, my_right, my_left, Tube_COMM, gclr
    use magnetic_coordinates, only : mpol2
    use connection_condition, only : connection_condition_at_theta_cut0
    use mpi
    implicit none
    real(p_),intent(in) :: a(:,:,:)
    real(p_),intent(out) :: a_theta( size(a,1), size(a,2) )
    integer :: status(MPI_STATUS_SIZE), ierr, m, n
    real(p_) :: a_left(size(a,1),size(a,2))

    m = size(a,1)
    n = size(a,2)

    call MPI_Sendrecv(a(:,:,1), m*n, MPI_real8, my_right, 1, &
         &            a_left,   m*n, MPI_real8, my_left,  1, Tube_COMM, status, ierr)
    if(GCLR==0)  call connection_condition_at_theta_cut0(a_left, mpol2-1)

    a_theta(:,:) = (a(:,:,2) - a_left(:,:))/(two*dtheta2) !centered difference
  end subroutine theta_derivative

!!$  subroutine xy_components_theta_derivative(ex,ey,ex_theta,ey_theta) !partial derivative with respect to theta (with psi and alpha fixed), i.e., along the magnetic field line
!!$    use constants,only:two
!!$    use domain_decomposition,only: dtheta2,myid
!!$    use domain_decomposition,only: GCLR,ntube,TUBE_COMM,my_left,my_right,GCLR_cut
!!$    use connection_condition,only: connection_condition_at_theta_cut, connection_condition_at_theta_cut3,&
!!$         & components_transformation_at_theta_cut, components_transformation_at_theta_cut2
!!$    use mpi
!!$    implicit none
!!$    real(p_),intent(in):: ex(:,:),ey(:,:)
!!$    real(p_),intent(out):: ex_theta(:,:),ey_theta(:,:)
!!$    real(p_):: ex_left(size(ex,1),size(ex,2)),ey_left(size(ey,1),size(ey,2))
!!$    real(p_):: ex_right(size(ex,1),size(ex,2)),ey_right(size(ey,1),size(ey,2))
!!$    integer:: status(MPI_STATUS_SIZE),ierr
!!$    integer:: m,n
!!$
!!$
!!$    if(any(shape(ex).ne.shape(ex_theta)) .or. &
!!$         & any(shape(ey).ne.shape(ey_theta)) .or. &
!!$         & any(shape(ex).ne.shape(ey))) stop
!!$
!!$    m=size(ex,1)
!!$    n=size(ex,2)
!!$
!!$    call MPI_Sendrecv(ex, m*n, MPI_real8, my_right, 1,&
!!$         &            ex_left,m*n, MPI_real8, my_left, 1,Tube_COMM,status,ierr)
!!$
!!$    call MPI_Sendrecv(ex, m*n, MPI_real8, my_left, 2,&
!!$         &            ex_right,m*n, MPI_real8, my_right, 2,Tube_COMM,status,ierr)
!!$
!!$    call MPI_Sendrecv(ey, m*n, MPI_real8, my_right, 3,&
!!$         &            ey_left,m*n, MPI_real8, my_left, 3,Tube_COMM,status,ierr)
!!$
!!$    call MPI_Sendrecv(ey, m*n, MPI_real8, my_left, 4,&
!!$         &            ey_right,m*n, MPI_real8, my_right, 4,Tube_COMM,status,ierr)
!!$
!!$    if(GCLR.eq.GCLR_cut) then 
!!$       call connection_condition_at_theta_cut(ex_right) 
!!$       call connection_condition_at_theta_cut(ey_right) 
!!$       call components_transformation_at_theta_cut(ex_right,ey_right)
!!$    elseif(GCLR.eq.GCLR_cut+1) then !It is a little difficult for me to obtain ex_left for No. GCLR_cut+1 pro., so I use the non-centered difference scheme
!!$       call connection_condition_at_theta_cut3(ex_left) 
!!$       call connection_condition_at_theta_cut3(ey_left) 
!!$       !call check_components_transformation(ex_left,ey_left,mpol-1)
!!$       call components_transformation_at_theta_cut2(ex_left,ey_left) !transform the components at theta=pi-dtheta2 to theta=-pi-dtheta2
!!$       !call check_components_transformation(ex_left,ey_left,-1)
!!$    endif
!!$
!!$    ex_theta=(ex_right-ex_left)/(two*dtheta2) !centered difference
!!$    ey_theta=(ey_right-ey_left)/(two*dtheta2)
!!$
!!$    !write(*,*) a_left(10,10),a_right(10,10),a_theta(10,10),'myid=',myid
!!$  end subroutine xy_components_theta_derivative

end module derivatives_in_xyz

