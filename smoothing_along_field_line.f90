module nearby_field_module
contains
  subroutine get_nearby_field_along_field_line(a, a_left,a_right, a_left2, a_right2)
    use constants, only : p_
    use magnetic_coordinates, only : mpol2
    use domain_decomposition, only : GCLR,ntube,TUBE_COMM,my_left,my_right, my_left2,my_right2
    use connection_condition
    use mpi
    implicit none
    real(p_), intent(in) :: a(:,:)
    real(p_), intent(out) :: a_left(:,:), a_right(:,:)
    real(p_), optional, intent(out):: a_left2(:,:), a_right2(:,:)
    integer :: status(MPI_STATUS_SIZE), ierr, m, n

    m=size(a,1)
    n=size(a,2)

    call MPI_Sendrecv(a,      m*n, MPI_real8, my_right, 1, &
         &            a_left, m*n, MPI_real8, my_left, 1, Tube_COMM, status, ierr)

    call MPI_Sendrecv(a,       m*n, MPI_real8, my_left, 2, &
         &            a_right, m*n, MPI_real8, my_right, 2,Tube_COMM,status,ierr)


    if(GCLR==mpol2-1) then
       !call connection_condition_at_theta_cut(a_right) !a_right array received by GCLR_cut from its right neighbour can not be directly used by GCLR_cut processor because the a_right array grid-points are not at the same field lines expected by gclr_cut. Interpolation is needed. !Interpolating the numerical table defined on the theta=-pi plane to get values on the grids defined on the theta=+pi plane
       call connection_condition_at_theta_cut0(a_right, 0)
    endif


    if(GCLR==0) then
       !call connection_condition_at_theta_cut3(a_left)
       call connection_condition_at_theta_cut0(a_left, mpol2-1)
    endif

    if(present(a_left2) .neqv. present(a_right2)) stop "a_left2 and a_right2 must be both present or both absent"
    if( .not. present(a_left2)) return
    call MPI_Sendrecv(a,       m*n, MPI_real8, my_right2, 3,&
         &            a_left2, m*n, MPI_real8, my_left2,  3,Tube_COMM,status,ierr)

    call MPI_Sendrecv(a,        m*n, MPI_real8, my_left2,  4,&
         &            a_right2, m*n, MPI_real8, my_right2, 4,Tube_COMM,status,ierr)

    if(GCLR==mpol2-2) then
       call connection_condition_at_theta_cut0(a_right2, 0)
    endif

    if(GCLR==mpol2-1) then
       call connection_condition_at_theta_cut0(a_right2, 1)
    endif

    if(GCLR==1) then
       call connection_condition_at_theta_cut0(a_left2, mpol2-1)
    endif
    if(GCLR==0) then
       call connection_condition_at_theta_cut0(a_left2, mpol2-2)
    endif

  end subroutine get_nearby_field_along_field_line

end module nearby_field_module


module smoothing_module
  use constants,only:p_
  use constants,only:one,two,six
  use nearby_field_module
  implicit none
  real(p_),parameter:: weight1=0.5_p_,weight2=-one/six !got to know these values from Yang Chen (GEM code)

contains

subroutine smoothing_along_field_line_core(a) !smoothing along theta with psi and alpha fixed, i.e., along the magnetic field line
  real(p_),intent(inout):: a(:,:)
  real(p_):: a_left(size(a,1),size(a,2)),a_right(size(a,1),size(a,2))

  call get_nearby_field_along_field_line(a,a_left,a_right) !get value of field on the two grids that are to the left/rightof the present grid

  a=(weight1*a_right + a + weight1*a_left)/(one+two*weight1) !smoothing using weight1, get to know this smoothing scheme from Yang Chen
!the smoothing can be appled multiple times (possibley with different weights):
  call get_nearby_field_along_field_line(a,a_left,a_right) !get value of field on the two grids that are to the left/right of the present grid

  a=(weight2*a_right+a+weight2*a_left)/(one+two*weight2) !smoothing using weight2

end subroutine smoothing_along_field_line_core


subroutine smoothing_along_field_line_core5(a) !smoothing along theta with psi and alpha fixed, i.e., along the magnetic field line
  real(p_),intent(inout):: a(:,:)
  real(p_):: a_left(size(a,1),size(a,2)),a_right(size(a,1),size(a,2))
  real(p_):: a_left2(size(a,1),size(a,2)),a_right2(size(a,1),size(a,2))

  call get_nearby_field_along_field_line(a,a_left,a_right, a_left2, a_right2) !get value of field on the two grids that are to the left/rightof the present grid

  a=(-a_left2 + 4*a_left + 10*a + 4*a_right -a_right2)/16.

end subroutine smoothing_along_field_line_core5


!!$subroutine smoothing_em_field_along_field_line(smoothing_times)
!!$  use perturbation_field_matrix,only: epar=>epar_left,ex=>ex_left,ey=>ey_left
!!$  use perturbation_field_matrix,only:mf_par=>mf_par_left,mf_x=>mf_x_left,mf_y=>mf_y_left
!!$  use perturbation_field_matrix,only: epar_right,ex_right,ey_right,mf_par_right,mf_x_right,mf_y_right
!!$  use perturbation_field_matrix,only: epar_left=>epar_left2,ex_left=>ex_left2,ey_left=>ey_left2
!!$  use perturbation_field_matrix,only: mf_par_left=>mf_par_left2,mf_x_left=>mf_x_left2,mf_y_left=>mf_y_left2
!!$  use communication_connection
!!$
!!$  implicit none
!!$  integer,intent(in):: smoothing_times
!!$  integer:: k
!!$
!!$  do k=1,smoothing_times
!!$     epar=(weight1*epar_right+epar+weight1*epar_left)/(one+two*weight1) !smoothing using weight1, get to know this smoothing scheme from Dr. Yang Chen
!!$     ex=(weight1*ex_right+ex+weight1*ex_left)/(one+two*weight1) !smoothing using weight1, get to know this smoothing scheme from Dr. Yang Chen
!!$     ey=(weight1*ey_right+ey+weight1*ey_left)/(one+two*weight1) !smoothing using weight1
!!$     mf_par=(weight1*mf_par_right+mf_par+weight1*mf_par_left)/(one+two*weight1) !smoothing using weight1
!!$     mf_x=(weight1*mf_x_right+mf_x+weight1*mf_x_left)/(one+two*weight1) !smoothing using weight1
!!$     mf_y=(weight1*mf_y_right+mf_y+weight1*mf_y_left)/(one+two*weight1) !smoothing using weight1
!!$
!!$     call update_efield_at_right_boundary_of_present_cell() !electric field value at right-boundary of the present cell is needed when pushing particles
!!$     call update_bfield_at_right_boundary_of_present_cell() !magnetic field value at right-boundary of the present cell is needed when pushing particles
!!$     call update_efield_at_second_left_boundary() !electric field value at the second left-boundary of the present cell is needed when smoothing and taking z derivatives
!!$     call update_bfield_at_second_left_boundary() !magnetic field value at the second left-boundary of the present cell is needed when smoothing and taking z derivatives
!!$
!!$     epar=(weight2*epar_right+epar+weight2*epar_left)/(one+two*weight2) !smoothing using weight2
!!$     ex=(weight2*ex_right+ex+weight2*ex_left)/(one+two*weight2) !smoothing using weight2
!!$     ey=(weight2*ey_right+ey+weight2*ey_left)/(one+two*weight2) !smoothing using weight2
!!$     mf_par=(weight2*mf_par_right+mf_par+weight2*mf_par_left)/(one+two*weight2) !smoothing using weight2
!!$     mf_x=(weight2*mf_x_right+mf_x+weight2*mf_x_left)/(one+two*weight2) !smoothing using weight2
!!$     mf_y=(weight2*mf_y_right+mf_y+weight2*mf_y_left)/(one+two*weight2) !smoothing using weight2
!!$
!!$     call update_efield_at_right_boundary_of_present_cell() !electric field value at right-boundary of the present cell is needed when pushing particles
!!$     call update_bfield_at_right_boundary_of_present_cell() !magnetic field value at right-boundary of the present cell is needed when pushing particles
!!$     call update_efield_at_second_left_boundary() !electric field value at the second left-boundary of the present cell is needed when smoothing and taking z derivatives
!!$     call update_bfield_at_second_left_boundary() !magnetic field value at the second left-boundary of the present cell is needed when smoothing and taking z derivatives
!!$  end do
!!$end subroutine smoothing_em_field_along_field_line


subroutine smoothing_along_field_line_for_adiabatic_electron_model()
  use perturbation_field_matrix,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
  use magnetic_coordinates,only: m=>mtor,n=>nflux2 
  use communication_connection
  
  call smoothing_along_field_line_core(ef_cyl_r_left) 
  call smoothing_along_field_line_core(ef_cyl_z_left) 
  call smoothing_along_field_line_core(ef_cyl_phi_left)
  call communicate_field_value_between_neighbour_cells2() !for adiabatic electrons model
end subroutine smoothing_along_field_line_for_adiabatic_electron_model

end module smoothing_module
