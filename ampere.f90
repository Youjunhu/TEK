module ampere
  use constants, only: p_
  implicit none
  complex(p_), allocatable :: mpara(:,:,:) !coefficient matrix for parallel Ampere's law
  integer, allocatable :: ipiv(:,:)

contains
  subroutine prepare_ampere_matrix()
    use constants, only: zero, one, two, pi, twopi, kev, epsilon0, c
    use gk_module, only: mass_gk, charge_gk, nsm, gk_flr
    use magnetic_coordinates, only: nrad, nrad, toroidal_range, &
         & xgrid, xlow, xupp, grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use domain_decomposition, only: ipol_eq
    use gk_profile_funcs, only : gkn_func
    use control_parameters, only : nh_max
    implicit none
    real(p_) :: gy0, gx0, gxdgy0
    integer :: j, jp, m, n, nx, jeq, ns, info
    real(p_) :: part1, lx, ly, x
    real(p_) :: skin(nrad-2), wp2(nrad, nsm)
    complex(p_) :: s, part2
    complex(p_), parameter :: ii = (0._p_, 1._p_)

    nx = nrad - 1
    allocate(mpara(nx-1,nx-1,0:nh_max))
    allocate(ipiv(nx-1,0:nh_max))

    do ns = 1, nsm
       do j = 1, nrad
          x = xgrid(j)
          wp2(j,ns) = gkn_func(x,ns)*charge_gk(ns)**2/(mass_gk(ns)*epsilon0) !plasma frequency**2 
       enddo
    enddo

    do j = 1, nrad-2
       skin(j) = -sum(wp2(1+j,:))/c**2  !sum over species, skin current coefficient
    enddo

    lx = xupp - xlow
    ly = toroidal_range
    do n = 0, nh_max !toroidal harmonics
       do j = 1, nx-1
          jeq = j + 1
          gx0 = grad_psi(ipol_eq, jeq)
          gy0 = grad_alpha(ipol_eq, jeq)
          gxdgy0 = grad_psi_dot_grad_alpha(ipol_eq, jeq)
          do jp = 1, nx-1
             s = 0._p_
             do m = 1, nx-1
                part1 = (m*pi/lx)**2*gx0**2 + (n*twopi/ly)**2*gy0**2
                part2 = -ii*n*twopi/ly*(m*pi/lx)*2*gxdgy0
                s = s + sin(jp*m*pi/nx) * (part1*sin(j*m*pi/nx) + part2*cos(j*m*pi/nx))
             enddo
             mpara(j,jp,n) = s*two/nx  
          enddo
       enddo
    enddo

    do j=1,nx-1
       mpara(j,j,:) = mpara(j,j,:) - skin(j)
    enddo

    do n=0, nh_max !for each toroidal Fourier component
       call ZGETRF(nx-1, nx-1, mpara(:,:,n), nx-1, ipiv(:,n), info) !LU factorize the radial coeff matrix
    enddo
  end subroutine prepare_ampere_matrix


  subroutine solve_ampere(isolve, jpar_left, jpar_right, apara_s, apara_h, apara, ax, ay, az, ahx, ahy, ahz)
    use constants, only: p_, c, kev, epsilon0
    use normalizing, only: vu, tu, qu, nu
    use magnetic_coordinates, only: m=>mtor, n=>nrad, dtor,dradcor, jacobian
    use gk_module, only: w_unit
    use control_parameters, only: ismooth, nh_min, nh_max
    use derivatives_in_xyz, only: radial_derivative, toroidal_derivative, theta_derivative
    use communication_connection, only: communicate_between_neighbour_cells,&
         & update_field_at_right_boundary_of_present_cell
    use transform_module
    use poisson, only: merge_source
    use domain_decomposition, only: myid, numprocs, dtheta2, dvol
    use smoothing_module, only: smoothing_along_field_line_core5
    use math, only: ZGETRS_wrapper
    implicit none
    integer, intent(in) :: isolve
    real(p_), intent(in) :: jpar_left(:,:), jpar_right(:,:)
    real(p_), intent(in) :: apara_s(:,:,:)
    real(p_), intent(inout), dimension(:,:,:) :: apara_h
    real(p_), intent(out), dimension(:,:,:) :: apara, ax, ay, az, ahx, ahy, ahz
    real(p_) :: jpar(m,n), lap_As(m,n)
    real(p_) :: rhs(m,n-2), rhs0(m,n-2), solution(m,n-2), residual(m,n-2)
    complex(p_) :: rhs_dft(0:m-1,n-2), solution_dft(0:m-1,n-2)
    real(p_) :: debye_sq
    integer :: kn, i, j, jeq

    debye_sq = (epsilon0*tu*kev)/(nu*qu**2)

    call merge_source(jpar_left, jpar_right, jpar)
    do j = 1, n  
       jpar(:,j)  = jpar(:,j)*(w_unit/dvol(j)/nu)
    enddo

    call laplacian(apara_s(:,:,1), lap_As)

    rhs0(:,:) = lap_As(:,2:n-1) + jpar(:,2:n-1)/debye_sq*(vu/c)**2

    do i = 1, 1 !cancellation scheme
       !call skin_current_residual(i, isolve, apara_h, lap_As, jpar, residual)
       rhs(:,:) = rhs0(:,:) !+ residual
       call oned_fourier_transform1(rhs, rhs_dft, m,n-2) !calculating DFT along the first dimension
       solution_dft = 0._p_ !only some harmonics are solved, others are assumed to be zero
       do kn = nh_min, nh_max ! Solve the field equation for each toroidal harmonic
          call ZGETRS_wrapper(kn, mpara, IPIV, rhs_dft(kn,:), solution_dft(kn,:)) 
       enddo
       do kn = 1, nh_max ! negative toroidal mode number
          solution_dft(m-kn,:) = conjg(solution_dft(kn,:)) 
       enddo
       call oned_backward_fourier_transform1(solution_dft, solution, m, n-2)
       apara_h(:, 2:n-1, 1) = solution(:,:)
       apara_h(:, 1, 1) = 0 !zero boundary condition
       apara_h(:, n, 1) = 0
       call communicate_between_neighbour_cells(apara_h)
    enddo

    do i=1, ismooth
       call smoothing_along_field_line_core5(apara_h(:,:,1))
    enddo

    if(isolve==1) then
       call radial_derivative  (apara_h(:,:,1), ahx(:,:,1))
       call toroidal_derivative(apara_h(:,:,1), ahy(:,:,1))
       call theta_derivative   (apara_h(:,:,:), ahz(:,:,1))
       call update_field_at_right_boundary_of_present_cell(ahx,ahy,ahz)
    endif

    apara = apara_h + apara_s 
    call communicate_between_neighbour_cells(apara) 
    call radial_derivative  (apara(:,:,1), ax(:,:,1))
    call toroidal_derivative(apara(:,:,1), ay(:,:,1))
    call theta_derivative   (apara(:,:,:), az(:,:,1))
    call update_field_at_right_boundary_of_present_cell(ax,ay,az)
  end subroutine solve_ampere

  
  subroutine apara_s_evolution(aOld, aNew, phiz, dtao)
    use constants, only: p_
    use magnetic_coordinates, only: m=>mtor, n=>nrad
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: w2
    implicit none
    real(p_),intent(in) :: aOld(:,:), dtao
    real(p_), intent(in) :: phiz(:,:,:)
    real(p_),intent(out) :: aNew(:,:)
    real(p_) :: rate
    integer :: i, j, jeq
    ! use dapara_s/dt =-b\cdot\grad_phi
    do i = 1,m
       do j = 1,n
          jeq = j
          rate = -phiz(i,j,1) * w2(ipol_eq, jeq)
          aNew(i,j) = aOld(i,j) + rate * dtao
       enddo
    enddo
    !use dapara_s/dt =0, worse than the above, tested
    !aNew = aOld 
  end subroutine apara_s_evolution

  subroutine apara_resplit_and_weight_pullback(w_gk, apara_s, apara_h, ahx, ahy, ahz)
    use constants, only: p_, kev
    use gk_module, only:  nsm, vn_gk, mass_gk, v_gk, charge_gk, xgc, vpar_gk, nm_gk, &
         & zgc, x_ring, y_ring, z_ring, touch_bdry_gc, gk_flr, ptcl_num0_gk, w_unit, ps_vol_gk
    use normalizing, only: qu, tu, vu
    use load_gk_mod, only: maxwellian
    use gk_profile_funcs, only: gkt_func
    use gyro_average_mod, only: gyro_average0
    implicit none
    real(p_), intent(inout) :: w_gk(:,:), apara_s(:,:,:), apara_h(:,:,:)
    real(p_), intent(out) :: ahx(:,:,:), ahy(:,:,:), ahz(:,:,:)
    real(p_) :: x, te, ah_av, E, factor
    integer :: k, ns
    do ns = 1, nsm 
       do k = 1, nm_gk(ns)
          x = xgc(k,ns)
          te = gkt_func(x,ns)
          E = 0.5*mass_gk(ns)*(v_gk(k,ns)*vn_gk(ns))**2/kev
          call gyro_average0(gk_flr(ns), x_ring(:,k,ns), y_ring(:,k,ns), z_ring(:,k,ns), &
               &   touch_bdry_gc(k,ns), apara_h, ah_av)
          factor = (charge_gk(ns)/qu)/(te/tu)*(vn_gk(ns)/vu)
          !weight pullback
!w_gk(k,ns) = w_gk(k,ns) -factor*vpar_gk(k,ns)*ah_av*ptcl_num0_gk(k,ns)
w_gk(k,ns) = w_gk(k,ns) -factor*vpar_gk(k,ns)*ah_av*maxwellian(x,E,ns)* ps_vol_gk(k,ns)*vn_gk(ns)**3/w_unit
       enddo
    enddo
    !Apar resplit
    apara_s = apara_s + apara_h !collect apara_h into apara_s
    apara_h = 0
    ahx = 0
    ahy = 0
    ahz = 0
  end subroutine apara_resplit_and_weight_pullback


  subroutine laplacian(apara, out)
    use constants, only: p_
    use magnetic_coordinates, only: m=>mtor, n=>nrad,  &
         & grad_psi, grad_psi_dot_grad_alpha, grad_alpha
    use domain_decomposition, only: ipol_eq
    use derivatives_in_xyz, only: radial_derivative, toroidal_derivative
    implicit none
    real(p_), intent(in) :: apara(m,n)
    real(p_), intent(out) :: out(m,n)
    real(p_), dimension(m,n) :: apara_x, apara_xx, apara_y, apara_yy, apara_xy
    integer :: j, jeq

    call radial_derivative(apara, apara_x)
    call radial_derivative(apara_x, apara_xx)
    call toroidal_derivative(apara, apara_y)
    call toroidal_derivative(apara_y, apara_yy)
    call toroidal_derivative(apara_x, apara_xy)

    do j = 1, n
       jeq = j
       out(:,j) = grad_psi  (ipol_eq, jeq)**2*apara_xx(:,j)  &
            &   + grad_alpha(ipol_eq, jeq)**2*apara_yy(:,j)  &
            &   + 2*grad_psi_dot_grad_alpha(ipol_eq, jeq)*apara_xy(:,j)

    enddo
  end subroutine laplacian


  subroutine skin_current_residual(iter, isolve, apara_h, lap_As, jpar_ref, residual)
    use constants,only:zero,one,two,pi,twopi,kev,epsilon0, c, mu0
    use normalizing, only: qu, vu, nu, tu
    use magnetic_coordinates,only: m=>mtor, n=>nrad, xgrid, dtor, dradcor, &
         & ygrid, jacobian, xgrid
    use gk_profile_funcs, only : gkn_func, gkt_func
    use gk_module, only: mass_gk, charge_gk, nsm, vn_gk,  nm_gk, w_unit, ps_vol_gk, &
         & xgc, zgc, ygc, vpar_gk, xgc_mid, zgc_mid, ygc_mid, vpar_gk_mid, &
         & v_gk,  touch_bdry_gc, ptcl_num0_gk
    use load_gk_mod, only : maxwellian
    use poisson, only: merge_source
    use domain_decomposition,only: dtheta2,theta_start, ipol_eq, multi_eq_cells, myid, tclr, dvol
    use gyro_average_mod, only : field_at_particle0
    implicit none
    integer, intent(in) :: iter, isolve
    real(p_), intent(in) :: apara_h(:,:,:), lap_As(:,:), jpar_ref(:,:)
    real(p_), intent(out) :: residual(m,n-2)
    real(p_)  :: debye_sq, wp2, ah_e
    real(p_) :: jpar_left(m,n), jpar_right(m,n), jpar(m,n)
    real(p_) :: x, y, z, te, E, delta_f, w, vpar
    integer :: ns, i, j, k, i_plus1, j_plus1, jeq
    real(p_) :: cz1, cz2, cy1, cy2, cx1, cx2, kernel

    debye_sq = (epsilon0*tu*kev)/(nu*qu**2)
    do ns=1,nsm
       if(charge_gk(ns)<0) exit !select out electrons
    enddo

    jpar_left = 0; jpar_right = 0;
    do k = 1, nm_gk(ns)     ! Monte-Carlo integration
       if (touch_bdry_gc(k,ns) .eqv. .true.) cycle
       if(isolve==1) then
          x = xgc_mid(k,ns)
          y = ygc_mid(k,ns)
          z = zgc_mid(k,ns)
          vpar = vpar_gk_mid(k,ns)
       else
          x = xgc(k,ns)
          y = ygc(k,ns)
          z = zgc(k,ns)
          vpar = vpar_gk(k,ns)
       endif
       te = gkt_func(x,ns)
       E = 0.5*mass_gk(ns)*(v_gk(k,ns)*vn_gk(ns))**2/kev
       call field_at_particle0(x, y, z, touch_bdry_gc(k,ns), apara_h, ah_e)
       !delta_f = -charge_gk(ns)/(te*kev)*vpar*vn_gk(ns)*ah_e*(Tu*kev/(qu*vu))*maxwellian(x,E,ns)
       !w = ps_vol_gk(k,ns)*vn_gk(ns)**3 * delta_f / w_unit
       w = -charge_gk(ns)/(te*kev)*vpar*vn_gk(ns)*ah_e*(Tu*kev/(qu*vu))*ptcl_num0_gk(k,ns)
       kernel = w*(charge_gk(ns)/qu)*vpar*(vn_gk(ns)/vu)

       cz1 = (z-theta_start)/dtheta2 
       cz2 = one-cz1

       i = floor((y-ygrid(1))/dtor+1) 
       cy1 = (y-ygrid(i))/dtor
       cy2 = one-cy1

       j = floor((x-xgrid(1))/dradcor+1)
       cx1 = (x-xgrid(j))/dradcor
       cx2 = one-cx1

       i_plus1=i+1
       if(i.eq.m) i_plus1=1 !periodic condition
       j_plus1=j+1
       if(j.eq.n) cycle !marker is out of radial computational region

       jpar_left(i,j) = jpar_left(i,j) + kernel*cz2*cy2*cx2
       jpar_left(i_plus1,j) = jpar_left(i_plus1,j) + kernel*cz2*cy1*cx2
       jpar_left(i,j_plus1) = jpar_left(i,j_plus1) + kernel*cz2*cy2*cx1
       jpar_left(i_plus1,j_plus1) = jpar_left(i_plus1,j_plus1)+kernel*cz2*cy1*cx1

       jpar_right(i,j) = jpar_right(i,j) + kernel*cz1*cy2*cx2
       jpar_right(i_plus1,j) = jpar_right(i_plus1,j)+ kernel*cz1*cy1*cx2
       jpar_right(i,j_plus1) = jpar_right(i,j_plus1) + kernel*cz1*cy2*cx1
       jpar_right(i_plus1,j_plus1) = jpar_right(i_plus1,j_plus1)+kernel*cz1*cy1*cx1
    enddo

    call merge_source(jpar_left, jpar_right, jpar)
    do j=1,n  !divided by cell volume (gridpoint is the center of the cell)
       jpar(:,j)  = jpar(:,j)*(w_unit/dvol(j)/nu)
    enddo

    do j =1, n-2
       x = xgrid(j+1)
       wp2= gkn_func(x,ns)*charge_gk(ns)**2/(mass_gk(ns)*epsilon0)
       residual(:,j) = jpar(:,j+1)/debye_sq*(vu/c)**2 -(-wp2/c**2*apara_h(:,j+1, 1))
    enddo
    ! if(tclr==0) then
    !    if(isolve==1 .and. iter==2 ) write(*,'(A40, 2i4, 10ES14.4)') 'ipol_eq, iter, jpar1, jpar2=', iter, ipol_eq, &
    !         & jpar(m/2,n/2)/debye_sq*(vu/c)**2, -wp2/c**2*apara_h(m/2,n/2, 1), lap_As(m/2,n/2), jpar_ref(m/2,n/2)/debye_sq*(vu/c)**2
    ! endif
  end subroutine skin_current_residual
  
end module ampere
