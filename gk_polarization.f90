module gk_polarization
  use constants,only: p_
  implicit none
  save
  complex(p_),dimension(:,:,:),allocatable:: polarization
!!$  real(p_),allocatable:: sigma(:,:)
!!$  complex(p_),dimension(:,:,:),allocatable::  u,  vt
!!$  complex(p_),dimension(:,:,:),allocatable::  ut, v
!!$  real(p_),parameter:: singular_value_threshold=6d-4 !singular value s smaller than this will be revmoved, (1/s be replace by zero in the inverse of sigma matrix)
contains
  subroutine prepare_polarization_matrix() !output: polarization
    use constants,only:zero,one,two,pi,twopi,kev,epsilon0
    use control_parameters,only: fk_ions_switch, gk_species_switch, space_charge_switch, nh
    use normalizing, only : tu,qu, nu
    use gk_module,only: mass_e,charge_e, nsm, gk_flr
    use magnetic_coordinates,only: nflux2,mtor,toroidal_range, dradcor,theta_1d_array,dtheta,j_low2,&
         & radcor_1d_array, radcor_low2, radcor_upp2, grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use table_in_mc, only:     b_mc

    use domain_decomposition,only: ipol_eq
    
    use math,only: srcbes, bessi0
    use density_temperature_profile_mod, only : te_func, ne_func
    complex(p_),parameter:: ii=(0._p_,1._p_)
    integer:: nx, j, jp, jeq, m, ierr, ns
    integer:: i0,kn0,kn, info
    real(p_):: ktor, kx1,kx2, kper1_sq, kper2_sq, b1,b2
    real(p_):: radcor,bval,gx0, gy0,gxdgy0
    real(p_):: lx,x, ly
    real(p_):: omega, gyro_radius_sq, gamma1,gamma2, trash !omega, cycontron angular frequency, gyro_radius_sq, square of the gyro-radius
    complex(p_):: part1, part2, sum, space_charge1, space_charge2
    real(p_):: coefficient,   debye_length_sq
    complex(p_),dimension(:,:,:),allocatable:: mmm
    debye_length_sq=tu*kev*epsilon0/(nu*qu**2)
    !write(*,*) 'debye_length=',sqrt(debye_length_sq)
    nx=nflux2-1
    allocate(mmm(nx-1,nx-1,0:nh))
    allocate(polarization (nx-1,nx-1,0:nh))
    polarization=(0._p_,0._p_)

    !allocate(tmp0(nx-1,nx-1,0:nh))
!!$    allocate(sigma(nx-1,-nh:nh))
!!$    allocate(u(nx-1,nx-1,-nh:nh))
!!$    allocate(vt(nx-1,nx-1,-nh:nh))
!!$    allocate(ut(nx-1,nx-1,-nh:nh))
!!$    allocate(v(nx-1,nx-1,-nh:nh))

    lx=radcor_upp2-radcor_low2
    ly=toroidal_range
    do ns=1,nsm
       if(gk_flr(ns) .eqv. .false.) cycle !no polarization density
       mmm=(0._p_,0._p_)
       do j=1,nx-1
          jeq=j+j_low2
          radcor=radcor_1d_array(jeq)
          gx0=grad_psi(ipol_eq,jeq)
          gy0=grad_alpha(ipol_eq,jeq)
          gxdgy0=grad_psi_dot_grad_alpha(ipol_eq,jeq)
          bval=b_mc(ipol_eq,jeq)
          omega=bval*abs(charge_e(ns))/mass_e(ns)
          gyro_radius_sq=(te_func(radcor,ns)*kev/mass_e(ns))/omega**2
          x=radcor-radcor_low2
          coefficient=(charge_e(ns)/qu)**2/(te_func(radcor,ns)/tu)*(ne_func(radcor,ns)/nu)
          do kn=0,nh !toroidal harmonics
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
                   part1=  exp(ii*x*kx1)*(one-gamma1*gk_species_switch)
                   part2= -exp(ii*x*kx2)*(one-gamma2*gk_species_switch)
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
                   sum=sum+sin(jp*m*pi/nx)*(coefficient*(part1+part2) &
                        & +space_charge_switch*(space_charge1+space_charge2))
                enddo
                mmm(j,jp,kn)=sum/(two*ii)*two/nx  !the matrix corresponds to -np/n0, where np is the poloarization number density
             enddo
          enddo
       enddo
       polarization = polarization + mmm
    enddo
  end subroutine prepare_polarization_matrix

  subroutine prepare_polarization_matrix2() !output: polarization
    ! direct numerical computation of the poloial density matrix, refer to my note: nonlinear_gyrokinetic_equation.tm
    use constants,only:zero,one,two,pi,twopi,kev,epsilon0, four
    use normalizing,only:bn,ln,vn_i,beta_ni,tu,qu, nu
    use fk_module,only: mass_i,charge_i 
    use gk_module,only: gk_flr,mass_e,charge_e, nsm
    use magnetic_coordinates,only: nflux2,mtor,toroidal_range, dradcor,theta_1d_array,dtheta,j_low2,&
         & radcor_1d_array, radcor_low2, radcor_upp2, r_mc, z_mc, &
         & grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use domain_decomposition,only: ipol_eq
    use table_in_mc,only: b_mc
    use density_temperature_profile_mod, only : te_func, ne_func, ti_func, ni_func
    use func_in_mc, only : tor_shift_func
    use magnetic_field, only : radcor_as_func_of_pfn
    use control_parameters, only : nh
    use magnetic_field, only : pfn_func
    complex(p_),parameter:: ii=(0._p_,1._p_)
    integer::  nx, jp, j, jeq, m, ierr, ivper, j1, j2, ns
    integer:: kn
    real(p_):: ktor, radcor,bval,gx0
    real(p_):: dvper,vper_range, vper, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_):: intg
    real(p_):: coefficient, local_contribution,vthermal
    real(p_) :: R0, Z0, R, Z, delta_y, tor_shift0, angle1, angle2, c_left, c_right, radcor_tmp, theta_tmp
    integer,parameter :: nvper=150, n1=32, n2=32
    complex(p_),dimension(:,:,:),allocatable:: mmm

    nx=nflux2-2 !perturbations at the two end points are fixed at zero.
    allocate(mmm(nx,nx,0:nh)) !the matrix corresponds to -np/nu, where np is the poloarization number density
    allocate(polarization(nx,nx,0:nh)) !the matrix corresponds to -np/nu, where np is the poloarization number density
    vper_range = 3d0 !in unit of local vte=sqrt(2*T/m)
    dvper=vper_range/(nvper-1) 
    polarization=(0._p_,0._p_)    
    do ns=1,nsm
       if(gk_flr(ns) .eqv. .false.) cycle !no polarization density
       mmm=(0._p_,0._p_)
       do j=1,nx
          jeq=j+j_low2
          radcor=radcor_1d_array(jeq)
          gx0=grad_psi(ipol_eq,jeq)
          bval=b_mc(ipol_eq,jeq)
          omega=bval*abs(charge_e(ns))/mass_e(ns)
          coefficient=(charge_e(ns)/qu)**2/(te_func(radcor,ns)/tu)*(ne_func(radcor,ns)/nu)
          R0=r_mc(ipol_eq, jeq)
          Z0=z_mc(ipol_eq, jeq)
          tor_shift0=tor_shift_func(theta_1d_array(ipol_eq), radcor)
          vthermal=sqrt(2*te_func(radcor,ns)*kev/mass_e(ns))
          do ivper=1, nvper
             vper=0 + dvper*(ivper-1)
             gyro_radius=vper*vthermal/omega
             do j1=1, n1
                angle1=0.+twopi/(n1-1)*(j1-1)
                do j2=1, n2
                   angle2=0.+twopi/(n2-1)*(j2-1)
                   R=R0+gyro_radius*cos(angle1)+gyro_radius*cos(angle2)
                   Z=Z0+gyro_radius*sin(angle1)+gyro_radius*sin(angle2)
                   radcor_tmp=radcor_as_func_of_pfn(pfn_func(R,Z))
                   if((radcor_tmp.ge.radcor_upp2) .or. (radcor_tmp.le.radcor_low2)) cycle !no contribution to the polarization matrix
                   call interpolate_from_cylindrical_to_magnetic_coordinates1(R,Z,theta_tmp)
                   delta_y=tor_shift0-tor_shift_func(theta_tmp, radcor_tmp)
                   do kn=0,nh !toroidal harmonics
                      ktor=twopi*kn/toroidal_range !toroidal mode number
                      intg =dvper*vper*exp(-vper**2/two)*exp(ii*ktor*delta_y)
                      jp=int((radcor_tmp-radcor_low2)/dradcor) !locating
                      c_right=(radcor_tmp-radcor_1d_array(jp+j_low2))/dradcor !linear interpolating coefficient
                      c_left=one-c_right
                      if(jp == 0) then
                         mmm(j,jp+1,kn)= mmm(j,jp+1,kn)+intg*c_right
                      elseif(jp == nx) then
                         mmm(j,jp,kn)= mmm(j,jp,kn)+intg*c_left
                      elseif( (jp .gt. 0) .and. (jp .lt. nx)) then
                         mmm(j,jp,  kn)= mmm(j,jp,  kn)+intg*c_left   
                         mmm(j,jp+1,kn)= mmm(j,jp+1,kn)+intg*c_right
                      else
                         !do nothing
                      endif
                   enddo
                enddo
             enddo
          enddo
          mmm(j,:,:)= -coefficient*mmm(j,:,:)/(n1*n2)       
       enddo

       do j=1,nx
          jeq=j+j_low2
          radcor=radcor_1d_array(jeq)
          local_contribution=(ne_func(radcor,ns)/nu)*(charge_e(ns)/qu)**2/(te_func(radcor,ns)/tu)
          mmm(j,j, :) = mmm(j,j, :) + 1*local_contribution !add the non-gyro-averaging contribution to the diagonal elemets of the matrix
       enddo

       polarization = polarization + mmm !sum over species
    enddo
  end subroutine prepare_polarization_matrix2
end module gk_polarization
