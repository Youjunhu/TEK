module gk_weight_pusher
contains
pure  subroutine push_gk_weight(ns, nm, dtao, touch_bdry, radcor, vpar, v, &
       & xdrift0, ydrift0, zdrift0, zdrift00, mirror_force, xdrift1, ydrift1, zdrift1, &
       & phix, phiy, phiz, ahx, ahy, ahz, ah, w, w_new) 
    use constants, only: p_, two, kev
    use normalizing, only: vu, qu, tu
    use gk_module, only: charge_gk, mass_gk, vn_gk, ps_vol_gk, w_unit, ptcl_num0_gk, gk_nonlinear
    use gk_profile_funcs, only: gkt_func, gkn_func, gkdndx_func, gkdtdx_func
    use load_gk_mod, only: maxwellian
    use gk_polarization, only: slowing_down_xd, slowing_down_Ed
    implicit none
    integer, intent(in) :: ns, nm
    real(p_), intent(in) :: dtao
    logical, intent(in)  :: touch_bdry(nm)
    real(p_), intent(in) :: radcor(nm), vpar(nm), v(nm)
    real(p_), intent(in) :: xdrift0(nm), ydrift0(nm), zdrift0(nm), zdrift00(nm)
    real(p_), intent(in) :: mirror_force(nm)
    real(p_), intent(in) :: xdrift1(nm), ydrift1(nm), zdrift1(nm)
    real(p_), intent(in) :: phix(nm), phiy(nm), phiz(nm), ahx(nm), ahy(nm), ahz(nm), ah(nm)
    real(p_), intent(in) :: w(nm) !weight at t_n
    real(p_), intent(out) :: w_new(nm) !weight at t_n+dtao
    real(p_) :: gradient,  t0,n0, norm, norm2, rate, x, dfdx, dfdE, E
    integer :: k

!    if(ns==3) goto 1000

    do k = 1, nm  !for Maxwellian distribution
       if(touch_bdry(k) .eqv. .true.) then
          w_new(k) = 0._p_
       else
          x = radcor(k)
          E = 0.5*mass_gk(ns)*(v(k)*vn_gk(ns))**2/kev
          t0 = gkt_func(x, ns)
          n0 = gkn_func(x, ns)
          gradient = -gkdndx_func(x,ns)/n0 - gkdtdx_func(x,ns)/t0 &
               & *(v(k)**2/(two*t0*kev/(mass_gk(ns)*vn_gk(ns)**2))-1.5_p_)
          norm = (charge_gk(ns)/qu)/(t0/tu)
          norm2 = norm*vn_gk(ns)/vu

          rate = xdrift1(k)*gradient &
               & - norm*(phix(k)*xdrift0(k) +phiy(k)*ydrift0(k) +phiz(k)*zdrift00(k))  &
               & - norm*(phix(k)*xdrift1(k) +phiy(k)*ydrift1(k) +phiz(k)*zdrift1(k))*gk_nonlinear  & 
               & + norm2*vpar(k)*(ahx(k)*xdrift0(k)+ ahy(k)*ydrift0(k) + ahz(k)*zdrift0(k)) &
               & + norm2*vpar(k)*(ahx(k)*xdrift1(k)+ ahy(k)*ydrift1(k) + ahz(k)*zdrift1(k))*gk_nonlinear & 
               & - norm2*(-mirror_force(k))*ah(k)

          w_new(k) = w(k) + rate* maxwellian(x,E,ns)* ps_vol_gk(k,ns)*vn_gk(ns)**3/w_unit*dtao
          !w_new(k) = w(k) + rate* ptcl_num0_gk(k,ns)*dtao
       endif
    enddo

    return

1000 do k = 1, nm ! for slowing-down distribution
       if(touch_bdry(k) .eqv. .true.) then
          w_new(k) = 0._p_
       else
          x = radcor(k)
          E = 0.5*mass_gk(ns)*v(k)**2*vn_gk(ns)**2
          dfdx = slowing_down_xd(x, E)
          dfdE = slowing_down_Ed(x, E)
          norm = (charge_gk(ns)/qu)*tu*kev
          norm2 = norm*vn_gk(ns)/vu

          rate = - xdrift1(k) * dfdx &
               & + norm*(phix(k)*xdrift0(k)+phiy(k)*ydrift0(k)+phiz(k)*zdrift0(k))*dfdE  &
               & + norm*(phix(k)*xdrift1(k)+phiy(k)*ydrift1(k)+phiz(k)*zdrift1(k))*dfdE*gk_nonlinear  & 
               & - norm2*vpar(k)*(ahx(k)*xdrift0(k)+ ahy(k)*ydrift0(k) + ahz(k)*zdrift0(k))*dfdE &
               & - norm2*vpar(k)*(ahx(k)*xdrift1(k)+ ahy(k)*ydrift1(k) + ahz(k)*zdrift1(k))*dfdE*gk_nonlinear & 
               & + norm2*(-mirror_force(k))*ah(k)*dfdE
          w_new(k) = w(k) + rate*ps_vol_gk(k,ns)*vn_gk(ns)**3/w_unit*dtao
       endif
    enddo

  end subroutine push_gk_weight

end module gk_weight_pusher
