module derivatives_in_xyz
  use constants,only : p_, two
  implicit none

contains
  
subroutine x_derivative0(field, field_x)
    use magnetic_coordinates, only: dradcor
    complex(p_), intent(in) :: field(:,:)
    complex(p_), intent(out) :: field_x(:,:)
    integer :: m, n, i, j

    m=size(field,1)
    n=size(field,2)
    do i = 1, m
       do j = 2, n-1
          field_x(i,j) = (field(i,j+1)-field(i,j-1))/(two*dradcor)
       enddo
       field_x(i,1) = two*field_x(i,2)-field_x(i,3) !linear interpolation to obtain the derivative at the boundary point
       field_x(i,n) = two*field_x(i,n-1)-field_x(i,n-2) !linear interpolation
    enddo
  end subroutine x_derivative0
  
  subroutine x_derivative(field, field_x)
    !calculate the derivative of field with respect to psi in field-line-following-coordinates (psi,theta,alpha)
    use magnetic_coordinates, only: dradcor
    real(p_), intent(in) :: field(:,:)
    real(p_), intent(out) :: field_x(:,:)
    integer :: m, n, i, j

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

  end subroutine x_derivative

  subroutine y_derivative(field, field_y)
    !partial derivative of field with respect to y
    use magnetic_coordinates,only: dtor
    real(p_), intent(in) :: field(:,:)
    real(p_), intent(out) :: field_y(:,:)
    integer :: m,n,i,j,i_left,i_right

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

  end subroutine y_derivative

  subroutine z_derivative(a, a_theta) ! calculating derivative along a magnetic field line
    use domain_decomposition, only: dtheta2, myid, my_right, my_left, Tube_COMM, gclr
    use magnetic_coordinates, only: mpol2
    use communication_connection, only: connect_condition_across_theta_cut
    use mpi
    implicit none
    real(p_), intent(in) :: a(:,:,:)
    real(p_), intent(out) :: a_theta( size(a,1), size(a,2) )
    integer :: m, n, status(MPI_STATUS_SIZE), ierr
    real(p_) :: a_left(size(a,1),size(a,2))

    m = size(a,1)
    n = size(a,2)

    call MPI_Sendrecv(a(:,:,1), m*n, MPI_real8, my_right, 1, &
         &            a_left,   m*n, MPI_real8, my_left,  1, Tube_COMM, status, ierr)

    if(GCLR == 0) call connect_condition_across_theta_cut(a_left, -1)
    a_theta(:,:) = (a(:,:,2) - a_left(:,:))/(two*dtheta2) ! centered finite difference
  end subroutine z_derivative

  subroutine z_derivative0(a, a_theta) ! calculating derivative along a magnetic field line
    use domain_decomposition, only: dtheta2
    implicit none
    real(p_), intent(in) :: a(:, :, :)
    real(p_), intent(out) :: a_theta(:, :)

    a_theta(:,:) = (a(:,:,2) - a(:,:,1))/dtheta2 ! forward difference
  end subroutine z_derivative0

end module derivatives_in_xyz
