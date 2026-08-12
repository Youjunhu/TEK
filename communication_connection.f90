module communication_connection
  ! Communicate field value between neighbour theta cells, and handle the connecting condition across the theta cut
  use constants, only: p_, twopi
  implicit none
contains

  subroutine connect_condition_across_theta_cut(u, direction)
    use interpolate_module, only: linear_1d_interpolate
    use magnetic_coordinates, only: ygrid, toroidal_range !tor_shift_mc, mpol
    use radial_module, only: qrad
    use math, only: shift_toroidal
    real(p_), intent(inout) :: u(:,:)
    integer, intent(in) :: direction ! alowable values: +1 or -1
    !If direction == +1, follow magnetic field line along the positive direction of theta
    !If direction == -1, follow magnetic field line along the negative direction of theta
    real(p_), allocatable :: u_old(:)
    real(p_) :: y0, twopiq
    integer :: i, j, m, n

    m = size(u,1) ! toroidal
    n = size(u,2) ! radial
    allocate(u_old(m+1))

    do j = 1, n !radial
       u_old(1:m) = u(:, j)
       u_old(m+1) = u(1, j) ! toroidal periodic condition

       !twopiq = tor_shift_mc(mpol,j) - tor_shift_mc(1,j)
       twopiq = twopi*qrad(j)
       do i = 1, m
          y0 = ygrid(i) + twopiq*direction
          call shift_toroidal(y0, toroidal_range) 
          call linear_1d_interpolate(m+1, ygrid, u_old, y0, u(i,j))  
       enddo
    enddo
  end subroutine connect_condition_across_theta_cut


  subroutine merge_source(my_jpar_left, my_jpar_right, jpar)
    use magnetic_coordinates, only: m=>mtor, n=>nrad
    use domain_decomposition, only: GRID_COMM, TUBE_COMM, GCLR, my_left, my_right
    use mpi
    real(p_), intent(in) :: my_jpar_left(m,n), my_jpar_right(m,n)
    real(p_), intent(out) :: jpar(m,n)
    real(p_) :: jpar_left(m,n), jpar_right(m,n), jpar_left0(m,n)
    integer :: status(MPI_STATUS_SIZE), ierr

    !summing over all those procs in a z-cell
    call MPI_ALLREDUCE(my_jpar_left,  jpar_left,  m*n, MPI_REAL8, MPI_SUM, GRID_COMM, ierr)
    call MPI_ALLREDUCE(my_jpar_right, jpar_right, m*n, MPI_REAL8, MPI_SUM, GRID_COMM, ierr)

    call MPI_Sendrecv(jpar_right, m*n, MPI_real8, my_right, 3, &
         &            jpar_left0, m*n, MPI_real8, my_left,  3, Tube_comm, status, ierr)

    if(GCLR==0) call connect_condition_across_theta_cut(jpar_left0, -1) 
    jpar(:,:) = jpar_left(:,:) + jpar_left0(:,:) !add the contribution from the neighbour cell

  end subroutine merge_source

  subroutine update_scalar_at_right_boundary(s) 
    !Every proc is response for one cell which has two boundaries.
    !Only the field on left-boundary is computed by the present proc.
    !The field on the right-boundary is received from the neighbour proc. This subroutine handle this communication and the edge case for gclr = mpol2-1
    !Note that the definition of cell here is different from the the definition of cell in PIC.
    !In PIC,the grids are the centers of the cells while here the grids are the boundaries of the cells
    use mpi
    use domain_decomposition, only: myid, GCLR, Tube_comm, my_left, my_right
    use magnetic_coordinates, only : mpol2
    real(p_), intent(inout) :: s(:,:,:)
    integer :: status(MPI_STATUS_SIZE), ierr, m, n

    m = size(s,1) !toroidal
    n = size(s,2) !radial
    call MPI_Sendrecv(s(:,:,1),  m*n,  MPI_real8, my_left,  41, &
         &            s(:,:,2),  m*n,  MPI_real8, my_right, 41, Tube_COMM, status, ierr)
    !special treatment at theta cut, to handle the phi-grid mismatch
    if(GCLR == mpol2-1) call connect_condition_across_theta_cut(s(:,:,2), 1)

  end subroutine update_scalar_at_right_boundary


  subroutine update_derivatives_at_right_boundary(ax, ay, az)
    ! Communication between neighbour cells.
    ! Every proc is response for one z cell, which has two boundaries.
    ! Only the field on left-boundary is computed by the present proc.
    ! The value on the right boundary is obtained by communicating with neighbour process.
    ! (Field value at right-boundary of the present cell is needed when pushing particle weights)
    use mpi
    use magnetic_coordinates, only: m=>mtor, n=>nrad, mpol2
    use domain_decomposition, only: myid, GCLR, Tube_comm, ntube, my_left, my_right
    real(p_), intent(inout) ::  ax(:,:,:), ay(:,:,:), az(:,:,:)
    integer :: status(MPI_STATUS_SIZE), ierr

    call MPI_Sendrecv(ax(:,:,1), m*n,  MPI_real8, my_left,  4, &
         &            ax(:,:,2), m*n,  MPI_real8, my_right, 4, Tube_COMM,status,ierr)

    call MPI_Sendrecv(ay(:,:,1), m*n,  MPI_real8, my_left,  5, &
         &            ay(:,:,2), m*n,  MPI_real8, my_right, 5, Tube_COMM,status,ierr)

    call MPI_Sendrecv(az(:,:,1), m*n,  MPI_real8, my_left,  6, &
         &            az(:,:,2), m*n,  MPI_real8, my_right, 6, Tube_COMM,status,ierr)

    if(GCLR == mpol2-1) then !special treatment at theta cut
       call connect_condition_across_theta_cut(ax(:,:,2), direction=1)
       call connect_condition_across_theta_cut(ay(:,:,2), direction=1) 
       call connect_condition_across_theta_cut(az(:,:,2), direction=1) 
       call derivative_transform_at_theta_cut(ax(:,:,2), ay(:,:,2)) !transform the components at theta=-pi to theta=+pi
    endif
  end subroutine update_derivatives_at_right_boundary


  subroutine derivative_transform_at_theta_cut(ax, ay)
    !Transform the components of gradient of a scalar field from theta=-pi to theta=pi
    !The basis vectors grad_alpha, in terms of which the gradients of phi and apara are decomposed, is discontinuous across the theta-cut. Therefore the components at one side need to be transformed to the components on the other side
    use magnetic_coordinates, only: mpol, nrad, grad_psi, grad_psi_dot_grad_alpha, &
         & grad_psi_r, grad_psi_z, grad_alpha_r, grad_alpha_z
    real(p_), intent(in) :: ay(:,:)
    real(p_), intent(inout) :: ax(:,:)
    real(p_) :: grad_psi0, gxdgy0, gxdgy1 !gxdgy1 is the value of gxdgy at theta=+pi
    integer :: j

    do j = 1, nrad
       grad_psi0 = grad_psi(1,j)
       gxdgy0 = grad_psi_dot_grad_alpha(1,j)
       gxdgy1 = grad_psi_r(1,j) * grad_alpha_r(mpol,j)  &
            & + grad_psi_z(1,j) * grad_alpha_z(mpol,j)

       ax(:,j) = ax(:,j) + ay(:,j)*(gxdgy0-gxdgy1)/grad_psi0**2
    enddo
  end subroutine derivative_transform_at_theta_cut


  subroutine communicate_field_value_between_neighbour_cells2() !use cylindrical components of the electric field
    use mpi
    use perturbation_field,only: ef_cyl_r_left, ef_cyl_z_left, ef_cyl_phi_left !already known before entering this subroutine
    use perturbation_field,only: ef_cyl_r_right, ef_cyl_z_right, ef_cyl_phi_right !as output
    use magnetic_coordinates,only: m=>mtor, n=>nrad, mpol2
    use domain_decomposition,only: myid, GCLR,Tube_comm,ntube, my_left,my_right

    integer:: status(MPI_STATUS_SIZE),ierr
    !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc, the field on the right-boundary is received from the neighbour proc. Note that the definition of cell here is different from the the definition of cell in PIC: the grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.

    call MPI_Sendrecv(ef_cyl_r_left,  (m+1)*n,  MPI_real8, my_left,  4,&
         &            ef_cyl_r_right, (m+1)*n,  MPI_real8, my_right, 4,Tube_COMM,status,ierr)
    call MPI_Sendrecv(ef_cyl_z_left,  (m+1)*n,  MPI_real8, my_left,  5,&
         &            ef_cyl_z_right, (m+1)*n,  MPI_real8, my_right, 5,Tube_COMM,status,ierr)
    call MPI_Sendrecv(ef_cyl_phi_left, (m+1)*n,  MPI_real8, my_left,  6,&
         &            ef_cyl_phi_right,(m+1)*n,  MPI_real8, my_right, 6,Tube_COMM,status,ierr)

    if(GCLR == mpol2-1) then !special treatment at theta cut
       call connect_condition_across_theta_cut(ef_cyl_r_right, direction=1) 
       call connect_condition_across_theta_cut(ef_cyl_z_right, direction=1) 
       call connect_condition_across_theta_cut(ef_cyl_phi_right, direction=1) 
    endif
  end subroutine communicate_field_value_between_neighbour_cells2


  subroutine get_nearby_field_along_field_line(a, a_left, a_right, a_left2, a_right2)
    !get value of field on the two grids that are to the left/right of the present grid
    use constants, only: p_
    use magnetic_coordinates, only: mpol2
    use domain_decomposition, only: GCLR, TUBE_COMM, my_left, my_right, my_left2, my_right2
    use mpi
    implicit none
    real(p_), intent(in) :: a(:,:)
    real(p_), intent(out) :: a_left(:,:), a_right(:,:)
    real(p_), optional, intent(out) :: a_left2(:,:), a_right2(:,:)
    integer :: status(MPI_STATUS_SIZE), ierr, m, n

    m=size(a,1)
    n=size(a,2)

    call MPI_Sendrecv(a,      m*n, MPI_real8, my_right, 1, &
         &            a_left, m*n, MPI_real8, my_left, 1, Tube_COMM, status, ierr)

    call MPI_Sendrecv(a,       m*n, MPI_real8, my_left, 2, &
         &            a_right, m*n, MPI_real8, my_right, 2,Tube_COMM, status,ierr)


    if(GCLR==mpol2-1) then
       call connect_condition_across_theta_cut(a_right, 1)
    endif

    if(GCLR==0) then
       call connect_condition_across_theta_cut(a_left, -1)
    endif

    if(present(a_left2) .neqv. present(a_right2)) stop "a_left2 and a_right2 must be both present or both absent"
    if( .not. present(a_left2)) return
    call MPI_Sendrecv(a,       m*n, MPI_real8, my_right2, 3, &
         &            a_left2, m*n, MPI_real8, my_left2,  3,Tube_COMM,status,ierr)

    call MPI_Sendrecv(a,        m*n, MPI_real8, my_left2,  4, &
         &            a_right2, m*n, MPI_real8, my_right2, 4,Tube_COMM,status,ierr)

    if(GCLR==mpol2-2) then
       call connect_condition_across_theta_cut(a_right2, 1)
    endif

    if(GCLR==mpol2-1) then
       call connect_condition_across_theta_cut(a_right2, 1)
    endif

    if(GCLR==1) then
       call connect_condition_across_theta_cut(a_left2, -1)
    endif
    if(GCLR==0) then
       call connect_condition_across_theta_cut(a_left2, -1)
    endif

  end subroutine get_nearby_field_along_field_line
  
end module communication_connection
