module gk_weight_pusher
contains
  subroutine push_gk_weight(ns, dtao, touch_bdry_e, radcor_e, vpar_e, v_e, &
       & xdrift0, ydrift0, zdrift0, zdrift0m, mirror_force, xdrift1,&
       & phix_e, phiy_e, phiz_e, ahx_e, ahy_e, ahz_e, ah_e, w_e, w_e_new) 
    use constants,only: p_, two, kev
    use normalizing,only: vu, qu, tu
    use gk_module,only: nm_gk, charge_e, mass_e, vn_e, ptcl_num0_e
    use density_temperature_profile_mod,only : kappa_ne_func, kappa_te_func, te_func
    use magnetic_field, only : minor_r_prime
    implicit none
    integer,intent(in)   :: ns
    real(p_),intent(in)  :: dtao
    logical,intent(in)   :: touch_bdry_e(:)
    real(p_),intent(in)  :: radcor_e(:), vpar_e(:), v_e(:),  w_e(:)
    real(p_),intent(in)  :: xdrift0(:), ydrift0(:), zdrift0(:), zdrift0m(:), mirror_force(:), xdrift1(:)
    real(p_),intent(in)  :: phix_e(:), phiy_e(:), phiz_e(:), ahx_e(:), ahy_e(:), ahz_e(:), ah_e(:)
    real(p_),intent(out) :: w_e_new(:) !weight at t_{n}+dtao
    real(p_) :: gradient,  temperature, norm, norm2, rate, x, n0
    integer :: k

!!$omp parallel do private(rate)
    do k=1,nm_gk(ns)
       if(touch_bdry_e(k) .eqv. .true.) then
          w_e_new(k) = 0._p_
       else
          x = radcor_e(k)
          temperature = te_func(x,ns)
          gradient = kappa_ne_func(x,ns) &
               & + kappa_te_func(x,ns)*(v_e(k)**2/(two*temperature*kev/(mass_e(ns)*vn_e(ns)**2))-1.5_p_)
          norm = (charge_e(ns)/qu)/(temperature/tu)
          norm2 = norm*vn_e(ns)/vu
          n0=ptcl_num0_e(k,ns)
          rate = xdrift1(k)*minor_r_prime(x)*gradient*n0 &
               & - norm*(phix_e(k)*xdrift0(k)+phiy_e(k)*ydrift0(k)+phiz_e(k)*zdrift0m(k))*n0  &
               & + norm2*vpar_e(k)*(ahx_e(k)*xdrift0(k)+ahy_e(k)*ydrift0(k)+ahz_e(k)*zdrift0(k))*n0 &
               & - norm2*(-mirror_force(k))*ah_e(k)*n0
          w_e_new(k) = w_e(k) + rate * dtao
       endif
    enddo
!!$omp end parallel do
  end subroutine push_gk_weight

end module gk_weight_pusher
