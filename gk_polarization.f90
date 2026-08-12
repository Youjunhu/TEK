module gk_polarization
  !use constants,only: p_
  implicit none
!!$  real(p_),allocatable:: sigma(:,:)
!!$  complex(p_),dimension(:,:,:),allocatable::  u,  vt
!!$  complex(p_),dimension(:,:,:),allocatable::  ut, v
!!$  real(p_),parameter:: singular_value_threshold=6d-4 !singular value s smaller than this will be revmoved, (1/s be replace by zero in the inverse of sigma matrix)
contains
  subroutine prepare_polarization_matrix(ns, mmm) 
    use constants, only: p_, zero,one,two,pi,twopi,kev,epsilon0, ii
    use control_parameters, only: space_charge_switch, nh_max
    use normalizing, only : tu,qu, nu
    use gk_module, only: mass_gk,charge_gk, gk_flr
    use magnetic_coordinates, only: nrad,mtor,toroidal_range, dradcor, dtheta,&
         & xgrid, xlow, xupp, grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use table_in_mc, only: b_mc
    use domain_decomposition, only: ipol_eq
    use math, only: srcbes, bessi0
    use gk_profile_funcs, only: gkt_func, gkn_func
    complex(p_), intent(out) :: mmm(:,:,0:)
    integer :: nx, j, jp, jeq, m, ns, kn
    real(p_) :: ky, kx1, kx2, kper1_sq, kper2_sq, b1, b2
    real(p_) :: radcor, bval, gx0, gy0, gxdgy0
    real(p_) :: lx, x, ly
    real(p_) :: omega, gyro_radius_sq, gamma1,gamma2, trash !omega, cycontron angular frequency
    complex(p_) :: part1, part2, sum, space_charge1, space_charge2
    real(p_) :: coefficient, debye_length_sq

    mmm = (0._p_,0._p_)
    if(gk_flr(ns) .eqv. .false.) return !no polarization density

    debye_length_sq=tu*kev*epsilon0/(nu*qu**2)
    !write(*,*) 'debye_length=',sqrt(debye_length_sq)
    nx = nrad - 1
    lx = xupp - xlow
    ly = toroidal_range

    do j = 1, nx-1
       jeq = j+1
       radcor = xgrid(jeq)
       x = radcor - xlow
       gx0 = grad_psi(ipol_eq, jeq)
       gy0 = grad_alpha(ipol_eq, jeq)
       gxdgy0 = grad_psi_dot_grad_alpha(ipol_eq, jeq)
       bval = b_mc(ipol_eq, jeq)
       omega = bval*abs(charge_gk(ns))/mass_gk(ns)
       gyro_radius_sq = (gkt_func(radcor,ns)*kev/mass_gk(ns))/omega**2
       coefficient = (charge_gk(ns)/qu)**2/(gkt_func(radcor,ns)/tu)*(gkn_func(radcor,ns)/nu)
       do kn = 0, nh_max !toroidal harmonics
          ky = kn*twopi/ly !the toroidal wavenumber
          do jp = 1, nx-1
             sum = 0._p_
             do m = 1, nx-1
                kx1 = m*pi/lx
                kx2 = -kx1
                kper1_sq = kx1**2*gx0**2 + two*kx1*ky*gxdgy0 + ky**2*gy0**2
                kper2_sq = kx2**2*gx0**2 + two*kx2*ky*gxdgy0 + ky**2*gy0**2
                b1 = kper1_sq*gyro_radius_sq
                b2 = kper2_sq*gyro_radius_sq
                call srcbes(b1, gamma1, trash)
                call srcbes(b2, gamma2, trash)
                part1 =  exp(ii*x*kx1)*(one-gamma1)
                part2 = -exp(ii*x*kx2)*(one-gamma2)
                space_charge1 =  exp(ii*x*kx1)*kper1_sq*debye_length_sq
                space_charge2 = -exp(ii*x*kx2)*kper2_sq*debye_length_sq
                !write(*,*) kper1_sq*debye_length_sq,kper2_sq*debye_length_sq
                !----use Pade approximation: gamma=I0(b)e^(-b)=1/(1+b). The results agree with the above more accurate treatment
!!$                part1= exp(ii*x*kx1)*(one-one/(one+b1)) 
!!$                part2=-exp(ii*x*kx2)*(one-one/(one+b2))
                !----pade approximation end----

                !--- bessi0(b1)*exp(-b1), sometime gave NaN, so not reliable, do not use it.
!!$                part1= exp( ii*m*pi*x/lx)*(one-bessi0(b1)*exp(-b1)) 
!!$                part2=-exp(-ii*m*pi*x/lx)*(one-bessi0(b2)*exp(-b2))
                !------
                sum = sum + sin(jp*m*pi/nx)*(coefficient*(part1+part2) &
                     & + space_charge_switch*(space_charge1+space_charge2))
             enddo
             mmm(j, jp, kn) = sum/(two*ii)*two/nx
          enddo
       enddo
    enddo

  end subroutine prepare_polarization_matrix

  subroutine prepare_polarization_matrix2(ns, mmm, nx) 
    ! direct numerical computation of the poloial density matrix,
    !refer to my note "nonlinear_gyrokinetic_equation.tm"
    use constants, only: p_, zero,one,two,pi,twopi,kev,epsilon0, four, ii
    use normalizing, only: bn,ln,tu,qu, nu
    use gk_module, only: gk_flr, mass_gk, charge_gk
    use magnetic_coordinates, only: nrad, mtor, toroidal_range, dradcor, zgrid, dtheta, &
         & xgrid,  r_mc, z_mc, tor_shift_mc, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: bphi_mc
    use magnetic_coordinates, only: grad_alpha_r, grad_alpha_z, grad_psi_r, grad_psi_z
    use gk_profile_funcs, only : gkt_func, gkn_func
    use func_in_mc, only: tor_shift_func
    use control_parameters, only: nh_min, nh_max
    use magnetic_field, only: pfn_func
    use map_to_mc, only: interpolate_from_cylindrical_to_magnetic_coordinates1
    integer, intent(in) :: ns, nx
    complex(p_), intent(out) :: mmm(nx, nx, 0:nh_max)
    integer ::  j, jeq, jp0, jp1, ivper, j1, j2, kn
    real(p_) :: x, bphi
    real(p_) :: dvper,vper_range, vper, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg
    real(p_) :: coefficient, local_contribution, vt, xlow, xupp
    real(p_) :: R0, Z0, R, Z, dr, dz, dxdr, dxdz
    real(p_) :: dydr, dydz, delta_y, tor_shift0, angle1, angle2, c0, c1, xring, theta_tmp
    integer, parameter :: nvper = 100, n1 = 16, n2 = 32

    mmm=(0._p_, 0._p_) !mmm matrix corresponds to -np/nu, where np is the poloarization number density
    if(gk_flr(ns) .eqv. .false.) return

    xlow = xgrid(2) 
    xupp = xgrid(nrad-1)
    vper_range = 3._p_ !in unit of local sqrt(T/m)
    dvper = vper_range/(nvper-1) 

    do j = 1, nx
       jeq = j + 1
       x = xgrid(jeq)
       bphi = bphi_mc(ipol_eq, jeq)
       omega = abs(bphi*charge_gk(ns))/mass_gk(ns)
       coefficient = (charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)*(gkn_func(x,ns)/nu)
       R0 = r_mc(ipol_eq, jeq)
       Z0 = z_mc(ipol_eq, jeq)
       tor_shift0 = tor_shift_mc(ipol_eq, jeq)
       dydr = grad_alpha_r(ipol_eq, jeq)
       dydz = grad_alpha_z(ipol_eq, jeq)
       dxdr = grad_psi_r(ipol_eq, jeq)
       dxdz = grad_psi_z(ipol_eq, jeq)
       vt = sqrt(gkt_func(x,ns)*kev/mass_gk(ns))
       do ivper = 1, nvper
          vper = 0 + dvper*(ivper-1)
          gyro_radius = vper*vt/omega
          do j1 = 1, n1
             angle1 = 0. + twopi/(n1-1)*(j1-1)
             do j2 = 1, n2
                angle2 = 0.0 + twopi/(n2-1)*(j2-1)
                dr = gyro_radius*cos(angle1) + gyro_radius*cos(angle2)
                dz = gyro_radius*sin(angle1) + gyro_radius*sin(angle2)
                xring = pfn_func(R0+dr, Z0+dz)
                !xring = x + dxdr*dr + dxdz*dz
                if((xring >= xupp) .or. (xring < xlow)) cycle !no contribution
                jp0 = 1 + floor((xring-xlow)/dradcor)
                jp1 = jp0 + 1
                jeq = jp0 + 1 !eauilibrium xgrid index corresponding to jp0
                c0 = (xring - xgrid(jeq))/dradcor 
                c1 = one - c0
                !call interpolate_from_cylindrical_to_magnetic_coordinates1(R, Z, theta_tmp) !wrong at theta cut
                !delta_y = tor_shift0 - tor_shift_func(theta_tmp, xring) !wrong at theta cut
                delta_y = dydr*dr + dydz*dz
                do kn = nh_min, nh_max !toroidal harmonics
                   intg = dvper*vper*exp(-vper**2/two)*exp(ii*kn*twopi/toroidal_range*delta_y)
                   mmm(j, jp0, kn) = mmm(j, jp0, kn) + intg*c1
                   mmm(j, jp1, kn) = mmm(j, jp1, kn) + intg*c0
                enddo
             enddo
          enddo
       enddo
       mmm(j,:,:)= -coefficient*mmm(j,:,:)/(n1*n2)       
    enddo

    do j = 1, nx
       jeq = j + 1
       x = xgrid(jeq)
       local_contribution = (gkn_func(x,ns)/nu)*(charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)
       mmm(j, j, :) = mmm(j, j, :) + local_contribution !the non-gyro-averaging contribution
    enddo

  end subroutine prepare_polarization_matrix2


  subroutine prepare_polarization_matrix20(ns, mmm, nx) 
    ! direct numerical computation of the poloial density matrix using linear interpolation
    use constants, only: p_, zero, one, two, pi, twopi, kev, ii
    use normalizing, only: tu, qu, nu
    use gk_module, only: gk_flr, mass_gk, charge_gk
    use magnetic_coordinates, only: nrad, mtor, nsegment, dradcor, zgrid, xgrid
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: b_mc, bphi_mc
    use gk_profile_funcs, only : gkt_func, gkn_func
    use control_parameters, only: nh_min, nh_max
    use gyro_ring_mod, only: gyro_ring_core2
    integer, intent(in) :: ns, nx
    complex(p_), intent(out) :: mmm(nx, nx, 0:nh_max)
    real(p_) :: x, y, z, b
    real(p_) :: dvper, vper_range, vper, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg
    real(p_) :: coefficient(nx) 
    real(p_) :: vt, xlow, xupp, delta_y, c0, c1
    integer, parameter :: nvper = 500, n1 = 32, n2 = 33
    real(p_) :: xg(n1), yg(n1) ! guiding-center location
    real(p_) :: xp(n2), yp(n2) ! particle location
    integer ::  i, j, jeq, jp0, jp1, kg, kp, kn

    mmm = (0._p_, 0._p_) !mmm matrix corresponds to -(np/nu)*(q/qu), where np is the poloarization number density
    if(gk_flr(ns) .eqv. .false.) return

    xlow = xgrid(2) 
    xupp = xgrid(nrad-1)
    vper_range = 5._p_ !in unit of local sqrt(T/m)
    dvper = vper_range/(nvper-1) 

    z = zgrid(ipol_eq)
    y = 0
    do j = 1, nx
       jeq = j + 1
       x = xgrid(jeq)
       b = b_mc(ipol_eq, jeq)
       omega = abs(b*charge_gk(ns))/mass_gk(ns)
       coefficient(j) = (charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)*(gkn_func(x,ns)/nu)
       vt = sqrt(gkt_func(x,ns)*kev/mass_gk(ns))
       do i = 1, nvper
          vper = 0 + dvper*(i-1)
          gyro_radius = vper*vt/omega
          call gyro_ring_core2(n1, x, y, z, gyro_radius, xg, yg)
          do kg = 1, n1
             call gyro_ring_core2(n2, xg(kg), yg(kg), z, gyro_radius, xp, yp)
             do kp = 1, n2
                jp0 = 1 + floor((xp(kp)-xlow)/dradcor)
                jp1 = jp0 + 1
                if(jp1 >= nx .or. jp0<=1) cycle !no contribution
                jeq = jp0 + 1 !eauilibrium xgrid index corresponding to jp0
                c0 = (xp(kp) - xgrid(jeq))/dradcor 
                c1 = one - c0
                delta_y = yp(kp) - y
                do kn = nh_min, nh_max !toroidal harmonics
                   intg = dvper*vper*exp(-vper**2/two)*exp(ii*kn*nsegment*delta_y)
                   mmm(j, jp0, kn) = mmm(j, jp0, kn) + intg*c1
                   mmm(j, jp1, kn) = mmm(j, jp1, kn) + intg*c0
                enddo
             enddo
          enddo
       enddo
       mmm(j,:,:) = -coefficient(j)*mmm(j,:,:)/(n1*n2)       
    enddo

    do j = 1, nx
       mmm(j, j, :) = mmm(j, j, :) + coefficient(j) ! the non-gyro-averaging contribution
    enddo
  end subroutine prepare_polarization_matrix20


  subroutine prepare_polarization_matrix3(ns, mmm) 
    ! direct numerical computation of the poloial density matrix,
    use constants, only: p_, zero,one,two,pi,twopi,kev,epsilon0, four, ii
    use normalizing, only: bn,ln,tu,qu, nu
    use gk_module, only: gk_flr, mass_gk, charge_gk
    use magnetic_coordinates, only: nrad, mtor, toroidal_range, dradcor, zgrid, dtheta, &
         & xgrid,  r_mc, z_mc, tor_shift_mc, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: bphi_mc
    use magnetic_coordinates, only: xlow, xupp, grad_alpha_r, grad_alpha_z, grad_psi_r, grad_psi_z
    use gk_profile_funcs, only : gkt_func, gkn_func
    use func_in_mc, only: tor_shift_func
    use control_parameters, only: nh_min, nh_max
    use magnetic_field, only: pfn_func
    use map_to_mc, only: interpolate_from_cylindrical_to_magnetic_coordinates1
    integer, intent(in) :: ns
    complex(p_), intent(out) :: mmm(:,: , 0:)
    integer ::  j, jeq, m, ivper, j1, j2, kn, k, nx
    real(p_) :: x, bphi
    real(p_) :: dvper, vper_range, vper, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg, sum
    real(p_) :: coefficient, local_contribution, vt, lx, ly
    real(p_) :: R0, Z0, R, Z, dr, dz, dxdr, dxdz
    real(p_) :: dydr, dydz, delta_y, tor_shift0, angle1, angle2, xring
    integer, parameter :: nvper = 50, n1 = 8, n2 = 8

    mmm=(0._p_, 0._p_) !mmm matrix corresponds to -np/nu, where np is the poloarization number density
    if(gk_flr(ns) .eqv. .false.) return

    nx = nrad -1
    lx = xupp - xlow
    ly = toroidal_range
    vper_range = 3._p_ !in unit of local sqrt(T/m)
    dvper = vper_range/(nvper-1) 
    do j = 1, nx-1
       jeq = j + 1
       x = xgrid(jeq)
       bphi = bphi_mc(ipol_eq, jeq)
       omega = abs(bphi*charge_gk(ns))/mass_gk(ns)
       coefficient = (charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)*(gkn_func(x,ns)/nu)
       R0 = r_mc(ipol_eq, jeq)
       Z0 = z_mc(ipol_eq, jeq)
       tor_shift0 = tor_shift_mc(ipol_eq, jeq)
       dydr = grad_alpha_r(ipol_eq, jeq)
       dydz = grad_alpha_z(ipol_eq, jeq)
       dxdr = grad_psi_r(ipol_eq, jeq)
       dxdz = grad_psi_z(ipol_eq, jeq)
       vt = sqrt(gkt_func(x,ns)*kev/mass_gk(ns))
       do ivper = 1, nvper
          vper = 0 + dvper*(ivper-1)
          gyro_radius = vper*vt/omega
          do j1 = 1, n1
             angle1 = 0. + twopi/(n1-1)*(j1-1)
             do j2 = 1, n2
                angle2 = 0. + twopi/(n2-1)*(j2-1)
                dr = gyro_radius*cos(angle1) + gyro_radius*cos(angle2)
                dz = gyro_radius*sin(angle1) + gyro_radius*sin(angle2)
                !xring = pfn_func(R0+dr, Z0+dz) - xlow
                xring = x + dxdr*dr + dxdz*dz - xlow
                if((xring > lx) .or. (xring < 0)) cycle !no contribution
                delta_y = dydr*dr + dydz*dz
                do kn = nh_min, nh_max !toroidal harmonics
                   intg = dvper*vper*exp(-vper**2/two)*exp(ii*kn*twopi/ly*delta_y)
                   do k = 1, nx-1
                      sum = 0
                      do m = 1, nx-1
                         sum = sum + intg*2/nx*sin(k*m*pi/nx)*sin(m*pi/lx*xring)
                      enddo
                      mmm(j, k, kn) = mmm(j, k, kn) + sum
                   enddo
                enddo
             enddo
          enddo
       enddo
       mmm(j,:,:) = -coefficient*mmm(j,:,:)/(n1*n2)       
    enddo

    do j = 1, nx-1
       jeq = j + 1
       x = xgrid(jeq)
       local_contribution = (gkn_func(x,ns)/nu)*(charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)
       mmm(j, j, :) = mmm(j, j, :) + local_contribution !the non-gyro-averaging contribution
    enddo

  end subroutine prepare_polarization_matrix3


  pure  subroutine prepare_slowing_down_polarization_matrix(mmm, mass, charge, e_cut) 
    ! direct numerical computation of the poloial density matrix for slowing-down distribution
    use constants,only: p_, zero,one,two,pi,twopi, Mev,epsilon0, four, kev
    use normalizing,only:bn,ln,tu,qu, nu
    use magnetic_coordinates,only: nrad,mtor,toroidal_range, dradcor,zgrid,dtheta,&
         & xgrid, xlow, xupp, r_mc, z_mc, &
         & grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use domain_decomposition,only: ipol_eq
    use table_in_mc,only: b_mc
    use func_in_mc, only : tor_shift_func
    use magnetic_field, only : radcor_as_func_of_pfn
    use control_parameters, only : nh_max
    use magnetic_field, only : pfn_func
    use map_to_mc, only : interpolate_from_cylindrical_to_magnetic_coordinates1

    complex(p_), intent(out) :: mmm(:,:,0:)
    real(p_), intent(in) :: mass, charge, e_cut
    complex(p_), parameter :: ii = (0._p_, 1._p_)
    integer ::  nx, jp, j, jeq, m, ierr, ivper, ivpar, j1, j2, ns
    integer :: kn
    real(p_) :: ky, x, bval, gx0, energy, v, v_cut
    real(p_) :: dvper, dvpar, vper, vpar, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg
    real(p_) :: local_contribution
    real(p_) :: R0, Z0, R, Z, delta_y, tor_shift0, angle1, angle2, c0, c1, x_tmp, theta_tmp
    integer, parameter :: nvper=50, nvpar=50, n1=8, n2=8

    nx = nrad - 2 !perturbations at the two end points are fixed at zero.
    !allocate(mmm(nx,nx,0:nh_max)) !the matrix corresponds to -np/nu, where np is the poloarization number density

    v_cut = sqrt(2*e_cut/mass)
    dvper = v_cut/(nvper-1)
    dvpar = 2*v_cut/(nvpar-1)

    mmm = (0._p_, 0._p_)
    do j = 1, nx
       jeq = j+1
       x = xgrid(jeq)
       gx0 = grad_psi(ipol_eq,jeq)
       bval = b_mc(ipol_eq,jeq)
       omega = bval*abs(charge)/mass
       R0 = r_mc(ipol_eq, jeq)
       Z0 = z_mc(ipol_eq, jeq)
       tor_shift0 = tor_shift_func(zgrid(ipol_eq), x)

       do ivpar = 1, nvpar
          vpar =  -v_cut + dvpar*(ivpar-1)
          do ivper = 1, nvper
             vper = 0 + dvper*(ivper-1)
             v = sqrt(vpar**2 + vper**2)
             energy = 0.5*mass*v**2
             if(v > v_cut) cycle
             gyro_radius = vper/omega
             do j1 = 1, n1
                angle1 = 0.+twopi/(n1-1)*(j1-1)
                do j2 = 1, n2
                   angle2 = 0.+twopi/(n2-1)*(j2-1)
                   R = R0+gyro_radius*cos(angle1)+gyro_radius*cos(angle2)
                   Z = Z0+gyro_radius*sin(angle1)+gyro_radius*sin(angle2)
                   x_tmp = radcor_as_func_of_pfn(pfn_func(R,Z))
                   if((x_tmp.ge.xupp) .or. (x_tmp.le.xlow)) cycle !no contribution to the polarization
                   call interpolate_from_cylindrical_to_magnetic_coordinates1(R,Z,theta_tmp)
                   delta_y = tor_shift0 - tor_shift_func(theta_tmp, x_tmp)

                   do kn = 0, nh_max !toroidal harmonics
                      ky = twopi*kn/toroidal_range !toroidal mode number
                      !intg = dvper*vper*exp(-vper**2/two)*exp(ii*ky*delta_y)
                      intg = dvpar * dvper * vper * slowing_down_Ed(x, energy) *exp(ii*ky*delta_y)
                      jp = int((x_tmp-xlow)/dradcor) !locating
                      c1=(x_tmp-xgrid(jp+1))/dradcor !linear interpolating coefficient
                      c0=one-c1
                      if(jp == 0) then
                         mmm(j, jp+1, kn)= mmm(j, jp+1, kn) + intg*c1
                      elseif(jp == nx) then
                         mmm(j, jp, kn)= mmm(j,jp,kn) + intg*c0
                      elseif( (jp .gt. 0) .and. (jp .lt. nx)) then
                         mmm(j,jp,  kn)= mmm(j, jp, kn) + intg*c0   
                         mmm(j,jp+1,kn)= mmm(j, jp+1, kn) + intg*c1
                      else
                         !do nothing
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       mmm(j,:,:)= - mmm(j,:,:)*twopi/(n1*n2)*(charge/qu)**2*tu*kev/nu
    enddo

    do j = 1, nx
       jeq = j+1
       x = xgrid(jeq)
       local_contribution = (charge/qu)**2*(f_Ed_integral(x))*tu*kev/nu
       mmm(j,j, :) = mmm(j,j, :) + 1*local_contribution !add the gyro-averaging-free part to the diagonal elements
    enddo

  end subroutine prepare_slowing_down_polarization_matrix


  pure real(p_) function slowing_down(x, E) result(f) !not used
    use constants,only: p_, kev
    use gk_radial_profiles, only : nalpha_object, alpha_normc, alpha_ecrit
    real(p_), intent(in) :: x, E
    real(p_) :: N, Ec

    N = alpha_normc%func(x)
    Ec = alpha_ecrit%func(x)*kev
    f = N/(1+sqrt(E/Ec)**3)

  end function slowing_down

  elemental real(p_) function slowing_down_Ed(x, E) result(f) 
    use constants,only: p_, kev
    use gk_radial_profiles, only : nalpha_object, alpha_normc, alpha_ecrit
    real(p_), intent(in) :: x, E
    real(p_) :: N, Ec

    N = alpha_normc%func(x)
    Ec = alpha_ecrit%func(x)*kev

    f = -N/(1+sqrt(E/Ec)**3)**2*1.5*sqrt(E/Ec)/Ec

  end function slowing_down_Ed

  pure real(p_) function f_Ed_integral(x) result(z)
    use constants,only: p_, kev, Mev, fourpi, atom_mass_unit
    implicit none
    real(p_), intent(in) :: x
    integer, parameter :: n = 100
    real(p_), parameter :: mass = 4*atom_mass_unit
    real(p_) :: vmax, v(n), E(n), dv
    integer :: i

    vmax = sqrt(2*3.5*Mev/mass)
    dv = vmax/(n-1)
    do i = 1, n
       v(i) = 0 + dv*(i-1)
    enddo
    E = 0.5*mass*v**2
    z =sum(slowing_down_Ed(x,E)*v**2)*fourpi*dv

  end function f_Ed_integral


  pure real(p_) function dNdx(x) result (z)
    use constants, only: p_
    use gk_radial_profiles, only : alpha_normc
    real(p_), intent(in) :: x
    real(p_), parameter :: dx = 1d-3

    z = (alpha_normc%func(x+dx) - alpha_normc%func(x-dx))/(2*dx)

  end function dNdx

  pure real(p_) function dEcdx(x) result (z) !SI
    use constants,only: p_, kev
    use gk_radial_profiles, only : alpha_ecrit
    real(p_), intent(in) :: x
    real(p_), parameter :: dx = 1d-3

    z = (alpha_ecrit%func(x+dx) - alpha_ecrit%func(x-dx))*kev/(2*dx)

  end function dEcdx


  pure real(p_) function slowing_down_xd(x, E) result(f) 
    use constants,only: p_, kev
    use gk_radial_profiles, only : alpha_normc, alpha_ecrit
    real(p_), intent(in) :: x, E
    real(p_) :: N, Ec

    N = alpha_normc%func(x)
    Ec = alpha_ecrit%func(x)*kev

    f = dNdx(x)/(1+sqrt(E/Ec)**3) - N/((1+sqrt(E/Ec)**3))**2*1.5*sqrt(E/Ec)*(-E/Ec**2*dEcdx(x))

  end function slowing_down_xd


end module gk_polarization
