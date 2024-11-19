module ampere
  use constants, only : p_
  implicit none
  save
  private
  complex(p_),allocatable :: mpara(:,:,:) !coefficient matrix for parallel Ampere's law
  integer, allocatable :: ipiv(:,:)
  public :: prepare_ampere_matrix, solve_ampere, apara_s_evolution, apara_resplit_and_weight_pullback
contains
  subroutine prepare_ampere_matrix
    use constants,only:zero,one,two,pi,twopi,kev,epsilon0, c
    use normalizing, only : tu,qu, nu
    use gk_module,only: mass_e,charge_e, nsm, gk_flr
    use magnetic_coordinates,only: nflux, nflux2, toroidal_range, j_low2, &
         & radcor_1d_array, radcor_low2, radcor_upp2, grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use domain_decomposition,only: ipol_eq
    use density_temperature_profile_mod, only : te_func, ne_func
    use control_parameters, only : nh
    implicit none
    real(p_) :: gy0, gx0, gxdgy0
    integer :: j, jp,m, n, nx, jeq, ns, info
    real(p_) :: part1, lx, ly, x
    real(p_), allocatable :: skin(:), wp2(:, :)
    complex(p_) :: s, part2
    complex(p_),parameter :: ii=(0._p_,1._p_)

    nx=nflux2-1
    allocate(mpara(nx-1,nx-1,0:nh), source=(zero,zero))
    allocate(ipiv(nx-1,0:nh))
    allocate(skin(nx-1))
    allocate(wp2(nflux, nsm))
    do ns=1,nsm
       do j=1,nflux
          x=radcor_1d_array(j)
          wp2(j,ns)=ne_func(x,ns)*charge_e(ns)**2/(mass_e(ns)*epsilon0) !plasma frequency**2 for each species
       enddo
    enddo

    do j=1,nx-1
       skin(j) = sum(wp2(j_low2+j,:))/c**2  !skin current coefficient
    enddo

    lx=radcor_upp2-radcor_low2
    ly=toroidal_range
    do j=1,nx-1
       jeq=j+j_low2
       gx0=grad_psi(ipol_eq,jeq)
       gy0=grad_alpha(ipol_eq,jeq)
       gxdgy0=grad_psi_dot_grad_alpha(ipol_eq,jeq)
       do n=0,nh !toroidal harmonics
          do jp=1,nx-1
             s = 0._p_
             do m=1,nx-1
                part1=(m*pi/lx)**2*gx0**2+(n*twopi/ly)**2*gy0**2
                part2=-ii*n*twopi/ly*(m*pi/lx)*2*gxdgy0
                s = s + sin(jp*m*pi/nx)*(part1*sin(j*m*pi/nx)+part2*cos(j*m*pi/nx))
             enddo
             mpara(j,jp,n) = s*two/nx  
          enddo
       enddo
    enddo

    do j=1,nx-1
       mpara(j,j,:) = mpara(j,j,:) + skin(j)
    enddo

    do n=0,nh !for each toroidal Fourier component
       call ZGETRF(nx-1, nx-1, mpara(:,:,n), nx-1, ipiv(:,n),info) !LU factorize the radial coeff matrix
    enddo
  end subroutine prepare_ampere_matrix


  subroutine solve_ampere()
    use constants, only : p_, c, kev, epsilon0
    use normalizing, only : vu, tu, qu, nu
    use magnetic_coordinates, only : m=>mtor, n=>nflux2
    use  perturbation_field_matrix, only : apara, apara_s, apara_h, &
         & ax, ay, az, ahx, ahy, ahz
    use control_parameters, only : ismooth, nh
    use derivatives_in_xyz, only : radial_derivative, &
         & toroidal_derivative, theta_derivative
    use communication_connection, only : communicate_field_value_between_neighbour_cells,&
         & update_field_at_right_boundary_of_present_cell
    use transform_module
    use smoothing_module, only : smoothing_along_field_line_core5
    use math, only : ZGETRS_wrapper
    implicit none
    real(p_) :: jpar(m,n), lap(m,n)
    real(p_) :: rhs(m,n-2), rhs0(m,n-2), solution(m,n-2)
    complex(p_) :: rhs_dft(0:m-1,n-2), solution_dft(0:m-1,n-2)
    real(p_), parameter :: debye = (epsilon0*tu*kev)/(nu*qu**2)
    integer :: kn, i
    integer,parameter :: niter=1

    call prepare_jpar_term(jpar)
    call laplacian(apara_s(:,:,1), lap)

    rhs0(:,:) = lap(:,2:n-1) + jpar(:,2:n-1)/debye*(vu/c)**2

    do i=1,niter
       rhs(:,:) = rhs0(:,:) +0
       call oned_fourier_transform1(rhs, rhs_dft, m,n-2) !calculating DFT of source1(:,:) along the first dimension
       solution_dft=0._p_ !initialized to zero
       do kn=nh,nh
          call ZGETRS_wrapper(kn, mpara, IPIV, rhs_dft(kn,:), solution_dft(kn,:)) !solve the field equation for a single toroidal harmonic
          solution_dft(m-kn,:) = conjg(solution_dft(kn,:)) !corresponding to the negative toroidal mode number
       enddo
       call oned_backward_fourier_transform1(solution_dft, solution, m, n-2)
       apara_h(:, 2:n-1, 1) = solution(:,:)
       apara_h(:, 1, 1) = 0 !zero boundary condition
       apara_h(:, n, 1) = 0
       call communicate_field_value_between_neighbour_cells(apara_h)
    enddo

    do i=1, ismooth
       call smoothing_along_field_line_core5(apara_h(:,:,1))
       call smoothing_along_field_line_core5(apara_s(:,:,1))
    enddo

    call radial_derivative  (apara_h(:,:,1),ahx(:,:,1))
    call toroidal_derivative(apara_h(:,:,1),ahy(:,:,1))
    call theta_derivative   (apara_h(:,:,:),ahz(:,:,1))
    call update_field_at_right_boundary_of_present_cell(ahx,ahy,ahz)

    apara = apara_h + apara_s 
    call communicate_field_value_between_neighbour_cells(apara) 
    call radial_derivative  (apara(:,:,1),ax(:,:,1))
    call toroidal_derivative(apara(:,:,1),ay(:,:,1))
    call theta_derivative   (apara(:,:,:),az(:,:,1))
    call update_field_at_right_boundary_of_present_cell(ax,ay,az)

!!$    ax=0 !testing
!!$    ay=0
!!$    az=0
!!$    apara_s=0
!!$    apara_h=0
!!$    ahx=0
!!$    ahy=0
!!$    ahz=0
  end subroutine solve_ampere

  subroutine apara_s_evolution(aOld,aNew,dtao)
    use constants, only : p_
    use magnetic_coordinates, only : m=>mtor, n=>nflux2, j_low2
    use domain_decomposition, only : ipol_eq
    use table_in_mc, only : w2
    use perturbation_field_matrix, only : phiz
    implicit none
    real(p_),intent(in) :: aOld(:,:), dtao
    real(p_),intent(out) :: aNew(:,:)
    real(p_) :: rate
    integer :: i, j, jeq

    do i = 1,m
       do j = 1,n
          jeq = j_low2 + (j-1)
          rate = -phiz(i,j,1) * w2(ipol_eq, jeq)
          aNew(i,j) = aOld(i,j) + rate * dtao
       enddo
    enddo

  end subroutine apara_s_evolution

  subroutine apara_resplit_and_weight_pullback()
    use constants, only : p_
    use gk_module, only :  nsm, vn_e, charge_e, w_e, radcor_e, vpar_e, nm_gk, &
         & theta_e, x_ring, y_ring,  touch_bdry_e, gk_flr, ptcl_num0_e
    use normalizing, only : qu, tu, vu
    use density_temperature_profile_mod,only : te_func
    use perturbation_field_matrix, only : apara_s, apara_h, ahx, ahy, ahz
    use gyro_average_mod, only : gyro_average0
    implicit none
    real(p_) :: x, te, ah_av
    integer :: k, ns

    do ns = 1, nsm
       do k = 1, nm_gk(ns)
          x = radcor_e(k,ns)
          te = te_func(x,ns)
          call gyro_average0(gk_flr(ns), theta_e(k,ns), x_ring(:,k,ns), y_ring(:,k,ns),&
               &   touch_bdry_e(k,ns), apara_h, ah_av)
          w_e(k,ns) = w_e(k,ns) -(charge_e(ns)/qu)/(te/tu)*(vn_e(ns)/vu)*vpar_e(k,ns)*ah_av*ptcl_num0_e(k,ns)  !pullback
       enddo
    enddo
    apara_s = apara_s + apara_h !collect apara_h into apara_s
    apara_h = 0
    ahx = 0
    ahy = 0
    ahz = 0
  end subroutine apara_resplit_and_weight_pullback


  subroutine laplacian(apara,out)
    use constants, only : p_
    use magnetic_coordinates, only : m=>mtor,n=>nflux2, j_low2, &
         & grad_psi, grad_psi_dot_grad_alpha, grad_alpha
    use domain_decomposition, only : ipol_eq
    use derivatives_in_xyz, only: radial_derivative, toroidal_derivative
    
    implicit none
    real(p_),intent(in) :: apara(m,n)
    real(p_),intent(out) :: out(m,n)
    real(p_),dimension(m,n) :: apara_x, apara_xx, apara_y, apara_yy, apara_xy
    integer :: i,j,jeq
    call radial_derivative(apara,apara_x)
    call radial_derivative(apara_x,apara_xx)
    call toroidal_derivative(apara,apara_y)
    call toroidal_derivative(apara_y,apara_yy)
    call toroidal_derivative(apara_x,apara_xy)
    do i=1,m
       do j=1,n
          jeq=j_low2+(j-1)
          out(i,j)=grad_psi  (ipol_eq, jeq)**2*apara_xx(i,j)  &
               & + grad_alpha(ipol_eq, jeq)**2*apara_yy(i,j)  &
               & + grad_psi_dot_grad_alpha(ipol_eq, jeq)*apara_xy(i,j)

       enddo
    enddo
  end subroutine laplacian

  subroutine prepare_jpar_term(jpar)
    use constants,only: p_
    use magnetic_coordinates,only:m=>mtor,n=>nflux2, j_low2, dtor,dradcor, jacobian
    use perturbation_field_matrix,only : my_jpar_e_left, my_jpar_e_right
    use domain_decomposition,only: GRID_COMM,TUBE_COMM, GCLR, GCLR_cut,my_left,my_right,&
         &  dtheta2,ipol_eq, multi_eq_cells
    use gk_module, only : w_unit
    use normalizing, only : nu
    use connection_condition
    use communication_connection
    use mpi
    implicit none
    real(p_), intent(out) :: jpar(m,n)
    real(p_) :: jpar_left(m,n), jpar_right(m,n), jpar_left0(m,n)
    integer :: j, jeq, ipol_eq1, status(MPI_STATUS_SIZE),ierr
    real(p_) :: dv1, dv2

    !summing over all those procs in a z-cell
    call MPI_ALLREDUCE(my_jpar_e_left,  jpar_left,  m*n, MPI_REAL8, MPI_SUM, GRID_COMM,ierr)
    call MPI_ALLREDUCE(my_jpar_e_right, jpar_right, m*n, MPI_REAL8, MPI_SUM, GRID_COMM,ierr)

    ipol_eq1=ipol_eq+multi_eq_cells
    do j=1,n  !divided by the spatial volume of a cell (a gridpoint is the center of the cell)
       jeq=j_low2+(j-1)
       dv1=abs(jacobian(ipol_eq, jeq))*dradcor*dtheta2*dtor !volume of the cell
       dv2=abs(jacobian(ipol_eq1,jeq))*dradcor*dtheta2*dtor 
       jpar_left(:,j)  = jpar_left(:,j)*(w_unit/dv1/nu)
       jpar_right(:,j) = jpar_right(:,j)*(w_unit/dv2/nu)
    enddo

    call MPI_Sendrecv(jpar_right, m*n, MPI_real8, my_right, 3,&
         &            jpar_left0, m*n, MPI_real8, my_left,  3, Tube_comm, status, ierr)

    if(GCLR.eq.GCLR_cut) call connection_condition_at_theta_cut_for_deposition(jpar_left0) 
    jpar_left(:,:) = jpar_left(:,:) + jpar_left0(:,:) !add the contribution from the neighbour cell
    jpar(:,:) = jpar_left(:,:) 

  end subroutine prepare_jpar_term

end module ampere
