module gk_polarization
  use constants,only: p_
  implicit none
!!$  real(p_),allocatable:: sigma(:,:)
!!$  complex(p_),dimension(:,:,:),allocatable::  u,  vt
!!$  complex(p_),dimension(:,:,:),allocatable::  ut, v
!!$  real(p_),parameter:: singular_value_threshold=6d-4 !singular value s smaller than this will be revmoved, (1/s be replace by zero in the inverse of sigma matrix)
contains
  subroutine prepare_polarization_matrix(ns, mmm) 
    use constants,only:zero,one,two,pi,twopi,kev,epsilon0
    use control_parameters,only: space_charge_switch, nh_max
    use normalizing, only : tu,qu, nu
    use gk_module,only: mass_gk,charge_gk, gk_flr
    use magnetic_coordinates,only: nrad,mtor,toroidal_range, dradcor, dtheta,&
         & xgrid, xlow, xupp, grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use table_in_mc, only:     b_mc
    use domain_decomposition,only: ipol_eq
    use math,only: srcbes, bessi0
    use gk_profile_funcs, only: gkt_func, gkn_func
    !use density_temperature_funcs, only : te_func, ne_func
    complex(p_), intent(out) :: mmm(:,:,0:)
    complex(p_), parameter :: ii=(0._p_,1._p_)
    integer :: nx, j, jp, jeq, m, ierr, ns
    integer :: i0,kn0,kn, info
    real(p_) :: ktor, kx1,kx2, kper1_sq, kper2_sq, b1,b2
    real(p_) :: radcor,bval,gx0, gy0,gxdgy0
    real(p_) :: lx,x, ly
    real(p_) :: omega, gyro_radius_sq, gamma1,gamma2, trash !omega, cycontron angular frequency, gyro_radius_sq, square of the gyro-radius
    complex(p_):: part1, part2, sum, space_charge1, space_charge2
    real(p_):: coefficient,   debye_length_sq

    mmm = (0._p_,0._p_)
    if(gk_flr(ns) .eqv. .false.) return !no polarization density
    
    debye_length_sq=tu*kev*epsilon0/(nu*qu**2)
    !write(*,*) 'debye_length=',sqrt(debye_length_sq)
    nx=nrad-1

    lx=xupp-xlow
    ly=toroidal_range

    do j=1,nx-1
       jeq=j+1
       radcor=xgrid(jeq)
       gx0=grad_psi(ipol_eq,jeq)
       gy0=grad_alpha(ipol_eq,jeq)
       gxdgy0=grad_psi_dot_grad_alpha(ipol_eq,jeq)
       bval=b_mc(ipol_eq,jeq)
       omega=bval*abs(charge_gk(ns))/mass_gk(ns)
       gyro_radius_sq=(gkt_func(radcor,ns)*kev/mass_gk(ns))/omega**2
       x=radcor-xlow
       coefficient=(charge_gk(ns)/qu)**2/(gkt_func(radcor,ns)/tu)*(gkn_func(radcor,ns)/nu)
       do kn=0,nh_max !toroidal harmonics
          ktor=twopi*kn/ly !the toroidal wavenumber
          do jp=1,nx-1
             sum=0._p_
             do m=1,nx-1
                kx1= m*pi/lx
                kx2=-kx1
                kper1_sq=kx1**2*gx0**2+two*kx1*ktor*gxdgy0+ktor**2*gy0**2
                kper2_sq=kx2**2*gx0**2+two*kx2*ktor*gxdgy0+ktor**2*gy0**2
                b1=kper1_sq*gyro_radius_sq
                b2=kper2_sq*gyro_radius_sq
                call srcbes(b1,gamma1,trash)
                call srcbes(b2,gamma2,trash)
                part1=  exp(ii*x*kx1)*(one-gamma1)
                part2= -exp(ii*x*kx2)*(one-gamma2)
                space_charge1=  exp(ii*x*kx1)*kper1_sq*debye_length_sq
                space_charge2= -exp(ii*x*kx2)*kper2_sq*debye_length_sq
                !write(*,*) kper1_sq*debye_length_sq,kper2_sq*debye_length_sq
                !----use Pade approximation: gamma=I0(b)e^(-b)=1/(1+b). The results agree with the above more accurate treatment
!!$                part1= exp(ii*x*kx1)*(one-one/(one+b1)) 
!!$                part2=-exp(ii*x*kx2)*(one-one/(one+b2))
                !----pade approximation end----

                !----! bessi0(b1)*exp(-b1) sometime give NaN, so not reliable, do not use it.
!!$                part1= exp( ii*m*pi*x/lx)*(one-bessi0(b1)*exp(-b1)) 
!!$                part2=-exp(-ii*m*pi*x/lx)*(one-bessi0(b2)*exp(-b2))
                !------
                sum = sum + sin(jp*m*pi/nx)*(coefficient*(part1+part2) &
                     & + space_charge_switch*(space_charge1+space_charge2))
             enddo
             mmm(j,jp,kn)=sum/(two*ii)*two/nx  !the matrix corresponds to -np/n0, where np is the poloarization number density
          enddo
       enddo
    enddo

  end subroutine prepare_polarization_matrix

  subroutine prepare_polarization_matrix2(ns, mmm) 
    ! direct numerical computation of the poloial density matrix,
    !refer to my note "nonlinear_gyrokinetic_equation.tm"
    use constants, only: zero,one,two,pi,twopi,kev,epsilon0, four
    use normalizing, only: bn,ln,tu,qu, nu
    use gk_module, only: gk_flr,mass_gk,charge_gk
    use magnetic_coordinates, only: nrad,mtor,toroidal_range, dradcor,zgrid,dtheta,&
         & xgrid, xlow, xupp, r_mc, z_mc, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: b_mc
    use gk_profile_funcs, only : gkt_func, gkn_func
    use func_in_mc, only: tor_shift_func
    use magnetic_field, only: radcor_as_func_of_pfn
    use control_parameters, only: nh_max
    use magnetic_field, only: pfn_func
    use map_to_mc, only: interpolate_from_cylindrical_to_magnetic_coordinates1
    integer, intent(in) :: ns
    complex(p_), intent(out) :: mmm(:,:,0:)
    complex(p_), parameter:: ii = (0._p_,1._p_)
    integer ::  nx, jp, j, jeq, m, ierr, ivper, j1, j2, kn
    real(p_) :: ktor, x,bval,gx0
    real(p_) :: dvper,vper_range, vper, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg
    real(p_) :: coefficient, local_contribution, vthermal
    real(p_) :: R0,Z0, R, Z, delta_y, tor_shift0, angle1, angle2, c_left, c_right, x_tmp, theta_tmp
    integer, parameter :: nvper=150, n1=32, n2=32
    
    mmm=(0._p_, 0._p_) !mmm matrix corresponds to -np/nu, where np is the poloarization number density
    if(gk_flr(ns) .eqv. .false.) return
    
    nx = nrad-2 !perturbations at the two end points are fixed at zero.
    vper_range = 3d0 !in unit of local vte=sqrt(2*T/m)
    dvper = vper_range/(nvper-1) 

       do j = 1, nx
          jeq = j + 1
          x = xgrid(jeq)
          gx0 = grad_psi(ipol_eq, jeq)
          bval = b_mc(ipol_eq,jeq)
          omega = bval*abs(charge_gk(ns))/mass_gk(ns)
          coefficient = (charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)*(gkn_func(x,ns)/nu)
          R0 = r_mc(ipol_eq, jeq)
          Z0 = z_mc(ipol_eq, jeq)
          tor_shift0 = tor_shift_func(zgrid(ipol_eq), x)
          vthermal = sqrt(gkt_func(x,ns)*kev/mass_gk(ns))
          do ivper = 1, nvper
             vper = 0 + dvper*(ivper-1)
             gyro_radius = vper*vthermal/omega
             do j1 = 1, n1
                angle1 = 0. + twopi/(n1-1)*(j1-1)
                do j2 = 1, n2
                   angle2 = 0. + twopi/(n2-1)*(j2-1)
                   R = R0 + gyro_radius*cos(angle1)+gyro_radius*cos(angle2)
                   Z = Z0 + gyro_radius*sin(angle1)+gyro_radius*sin(angle2)
                   x_tmp = pfn_func(R,Z)
                   if((x_tmp>=xupp) .or. (x_tmp<=xlow)) cycle !no contribution to the polarization matrix
                   call interpolate_from_cylindrical_to_magnetic_coordinates1(R, Z, theta_tmp)
                   delta_y = tor_shift0 - tor_shift_func(theta_tmp, x_tmp)
                   do kn = 0, nh_max !toroidal harmonics
                      ktor = twopi*kn/toroidal_range !toroidal mode number
                      intg = dvper*vper*exp(-vper**2/two)*exp(ii*ktor*delta_y)
                      jp = int((x_tmp-xlow)/dradcor) !locating
                      c_right = (x_tmp-xgrid(jp+1))/dradcor !interpolating coefficient
                      c_left = one - c_right
                      if(jp == 0) then
                         mmm(j, jp+1, kn) = mmm(j, jp+1, kn) + intg*c_right
                      elseif(jp == nx) then
                         mmm(j, jp, kn) = mmm(j, jp, kn) + intg*c_left
                      elseif( (jp .gt. 0) .and. (jp .lt. nx)) then
                         mmm(j,jp,  kn) = mmm(j, jp,  kn) + intg*c_left   
                         mmm(j,jp+1,kn) = mmm(j, jp+1,kn) + intg*c_right
                      else
                         !do nothing
                      endif
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


  pure  subroutine prepare_slowing_down_polarization_matrix(mmm, mass, charge, e_cut) 
    ! direct numerical computation of the poloial density matrix for slowing-down distribution
    use constants,only:zero,one,two,pi,twopi, Mev,epsilon0, four, kev
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
    real(p_) :: ktor, x, bval, gx0, energy, v, v_cut
    real(p_) :: dvper, dvpar, vper, vpar, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg
    real(p_) :: local_contribution
    real(p_) :: R0, Z0, R, Z, delta_y, tor_shift0, angle1, angle2, c_left, c_right, x_tmp, theta_tmp
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
                      ktor = twopi*kn/toroidal_range !toroidal mode number
                      !intg = dvper*vper*exp(-vper**2/two)*exp(ii*ktor*delta_y)
                      intg = dvpar * dvper * vper * slowing_down_Ed(x, energy) *exp(ii*ktor*delta_y)
                      jp = int((x_tmp-xlow)/dradcor) !locating
                      c_right=(x_tmp-xgrid(jp+1))/dradcor !linear interpolating coefficient
                      c_left=one-c_right
                      if(jp == 0) then
                         mmm(j, jp+1, kn)= mmm(j, jp+1, kn) + intg*c_right
                      elseif(jp == nx) then
                         mmm(j, jp, kn)= mmm(j,jp,kn) + intg*c_left
                      elseif( (jp .gt. 0) .and. (jp .lt. nx)) then
                         mmm(j,jp,  kn)= mmm(j, jp, kn) + intg*c_left   
                         mmm(j,jp+1,kn)= mmm(j, jp+1, kn) + intg*c_right
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
