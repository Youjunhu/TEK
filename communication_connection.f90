module communication_connection !communicate field value between neighbour cells and handle the connection condition near the theta cut
contains

  subroutine communicate_field_value_between_neighbour_cells(s) 
    use mpi
    use constants,only:p_
    use domain_decomposition,only:myid,GCLR,Tube_comm,ntube,GCLR_cut,GCLR_cut_left,my_left,my_right
    use connection_condition,only:  connection_condition_at_theta_cut
    use magnetic_coordinates, only : mpol2
    implicit none
    real(p_),intent(inout) :: s(:,:,:)
    integer:: status(MPI_STATUS_SIZE),ierr, m, n
    !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc, the field on the right-boundary is received from the neighbour proc. Note that the definition of cell here is different from the the definition of cell in PIC: the grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.
    m=size(s,1)
    n=size(s,2)
    call MPI_Sendrecv(s(:,:,1),  m*n,  MPI_real8, my_left,  41,&
         &            s(:,:,2),  m*n,  MPI_real8, my_right, 41,Tube_COMM,status,ierr)
    !special treatment at theta cut, to handle the phi_grid mismatch
    if(GCLR.eq.mpol2-1) call connection_condition_at_theta_cut(s(:,:,2)) 

  end subroutine communicate_field_value_between_neighbour_cells

  subroutine communicate_field_value_between_neighbour_cells2() !use cylindrical components of the electric field
    use mpi
    use constants,only:p_
    use perturbation_field_matrix,only: ef_cyl_r_left, ef_cyl_z_left, ef_cyl_phi_left !already known before entering this subroutine
    use perturbation_field_matrix,only: ef_cyl_r_right, ef_cyl_z_right, ef_cyl_phi_right !as output
    use magnetic_coordinates,only: m=>mtor,n=>nflux2
    use domain_decomposition,only:myid,GCLR,Tube_comm,ntube,GCLR_cut,GCLR_cut_left,my_left,my_right
    use connection_condition,only:  connection_condition_at_theta_cut
    implicit none

    integer:: status(MPI_STATUS_SIZE),ierr
    !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc, the field on the right-boundary is received from the neighbour proc. Note that the definition of cell here is different from the the definition of cell in PIC: the grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.

    call MPI_Sendrecv(ef_cyl_r_left,  (m+1)*n,  MPI_real8, my_left,  4,&
         &            ef_cyl_r_right, (m+1)*n,  MPI_real8, my_right, 4,Tube_COMM,status,ierr)
    call MPI_Sendrecv(ef_cyl_z_left,  (m+1)*n,  MPI_real8, my_left,  5,&
         &            ef_cyl_z_right, (m+1)*n,  MPI_real8, my_right, 5,Tube_COMM,status,ierr)
    call MPI_Sendrecv(ef_cyl_phi_left, (m+1)*n,  MPI_real8, my_left,  6,&
         &            ef_cyl_phi_right,(m+1)*n,  MPI_real8, my_right, 6,Tube_COMM,status,ierr)

    !if(GCLR.eq.GCLR_cut) then !special treatment at theta cut
    if(GCLR.eq.GCLR_cut_left) then !special treatment at theta cut
       call connection_condition_at_theta_cut(ef_cyl_r_right) 
       call connection_condition_at_theta_cut(ef_cyl_z_right) 
       call connection_condition_at_theta_cut(ef_cyl_phi_right) 
    endif

  end subroutine communicate_field_value_between_neighbour_cells2
  
  subroutine update_field_at_right_boundary_of_present_cell(ax,ay,az) !field value at right-boundary of the present cell is needed when pushing particle weights
    use mpi
    use constants,only:p_
    use magnetic_coordinates,only: m=>mtor,n=>nflux2
    use domain_decomposition,only:myid,GCLR,Tube_comm,ntube,GCLR_cut,GCLR_cut_left,my_left,my_right
    use connection_condition,only:connection_condition_at_theta_cut,components_transformation_at_theta_cut
    implicit none
    real(p_), intent(inout) ::  ax(:,:,:), ay(:,:,:), az(:,:,:)
    integer:: status(MPI_STATUS_SIZE),ierr
    !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc

    call MPI_Sendrecv(ax(:,:,1), m*n,  MPI_real8, my_left,  4,&
         &            ax(:,:,2), m*n,  MPI_real8, my_right, 4,Tube_COMM,status,ierr)

    call MPI_Sendrecv(ay(:,:,1), m*n,  MPI_real8, my_left,  5,&
         &            ay(:,:,2), m*n,  MPI_real8, my_right, 5,Tube_COMM,status,ierr)

    call MPI_Sendrecv(az(:,:,1), m*n,  MPI_real8, my_left,  6,&
         &            az(:,:,2), m*n,  MPI_real8, my_right, 6,Tube_COMM,status,ierr)

    if(GCLR.eq.GCLR_cut_left) then !special treatment at theta cut
       call connection_condition_at_theta_cut(ax(:,:,2)) !need interpolation because of the toroidal grid mismatch
       call connection_condition_at_theta_cut(ay(:,:,2)) 
       call connection_condition_at_theta_cut(az(:,:,2)) 
       call derivative_transform_at_theta_cut(ax(:,:,2),ay(:,:,2)) !transform the components at theta=-pi to theta=+pi
    endif
  end subroutine update_field_at_right_boundary_of_present_cell

  subroutine derivative_transform_at_theta_cut(ax,ay) !transform the components at theta=-pi to theta=pi
    !the basis vectors grad_alpha, in terms of which the gradient of apara/potential are decomposed, is discontinuous across the theta-cut. Therefore the components at one side need to be transformed to the components on the other side

    use constants,only:p_
    use magnetic_coordinates,only: mpol,nflux2,j_low2, &
     & grad_psi, grad_alpha,grad_psi_dot_grad_alpha, &
         & grad_psi_r,grad_psi_z, grad_alpha_r,grad_alpha_z,grad_alpha_phi
    implicit none
    real(p_),intent(in) :: ay(:,:)
    real(p_),intent(inout) :: ax(:,:)
    real(p_)::grad_psi_val,grad_alpha_val,grad_psi_dot_grad_alpha_val
    real(p_)::grad_alpha_dot_grad_alpha_positive,grad_psi_dot_grad_alpha_positive !_positive indicates value at theta=+pi, i.e., i=mpol
    integer:: j,jeq

    do j=1,nflux2
       jeq=j-1+j_low2
       grad_psi_val=grad_psi(1,jeq)
       grad_alpha_val=grad_alpha(1,jeq)
       grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha(1,jeq)
       grad_psi_dot_grad_alpha_positive=grad_psi_r(1,jeq)*grad_alpha_r(mpol,jeq)+&
            & grad_psi_z(1,jeq)*grad_alpha_z(mpol,jeq)

       grad_alpha_dot_grad_alpha_positive=grad_alpha_r(1,jeq)*grad_alpha_r(mpol,jeq)&
            & +grad_alpha_z(1,jeq)*grad_alpha_z(mpol,jeq)&
            & +grad_alpha_phi(1,jeq)*grad_alpha_phi(mpol,jeq)
       ax(:,j) = ax(:,j) + ay(:,j)*(grad_psi_dot_grad_alpha_val-grad_psi_dot_grad_alpha_positive)/grad_psi_val**2
    enddo
  end subroutine derivative_transform_at_theta_cut
  
  





!!$  subroutine update_efield_at_second_left_boundary() !get the value from the left neighbour
!!$    use constants,only:p_
!!$    use perturbation_field_matrix,only: ex=>ex_left,ey=>ey_left,epar=>epar_left !already known before entering this subroutine
!!$    use perturbation_field_matrix,only: ex_left2,ey_left2,epar_left2 !output
!!$    use magnetic_coordinates,only: m=>mtor,n=>nflux2,mpol,dtheta
!!$    use domain_decomposition,only:myid,GCLR,Tube_comm,ntube,GCLR_cut,my_left,my_right,dtheta2
!!$    use connection_condition
!!$    use mpi
!!$    implicit none
!!$    integer:: status(MPI_STATUS_SIZE),ierr
!!$
!!$    call MPI_Sendrecv(epar, (m+1)*n, MPI_real8, my_right,     10,&
!!$         &            epar_left2,(m+1)*n, MPI_real8, my_left, 10,Tube_COMM,status,ierr)
!!$
!!$    call MPI_Sendrecv(ex, (m+1)*n, MPI_real8, my_right,     11,&
!!$         &            ex_left2,(m+1)*n, MPI_real8, my_left, 11,Tube_COMM,status,ierr)
!!$
!!$    call MPI_Sendrecv(ey, (m+1)*n, MPI_real8, my_right,     12,&
!!$         &            ey_left2,(m+1)*n, MPI_real8, my_left, 12,Tube_COMM,status,ierr)
!!$
!!$    !  if(GCLR.eq.(GCLR_cut+1)) then
!!$    if(GCLR.eq.GCLR_cut) then
!!$       call connection_condition_at_theta_cut3(epar_left2) 
!!$       call connection_condition_at_theta_cut3(ex_left2) 
!!$       call connection_condition_at_theta_cut3(ey_left2) 
!!$       !call check_components_transformation(ex_left2,ey_left2,mpol-NINT(dtheta2/dtheta))
!!$       call components_transformation_at_theta_cut2(ex_left2,ey_left2) !transform the components at theta=pi-dtheta to theta=-pi-dtheta
!!$       !call check_components_transformation(ex_left2,ey_left2,-1)
!!$    endif
!!$
!!$  end subroutine update_efield_at_second_left_boundary





end module communication_connection
