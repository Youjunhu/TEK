module density_temperature_profile_mod
  use constants,only: p_
  implicit none
  real(p_),parameter :: a=0.6_p_ !meter, minor radius
  real(p_),parameter :: r0=a*0.5_p_ !meter, radial center of the simulation box
  real(p_),parameter :: te_scale=0.3_p_ !dimension-less
  real(p_),parameter :: ne_scale=0.3_p_ !dimension-less
  real(p_),parameter :: ti_scale=0.3_p_ !dimension-less
  real(p_),parameter :: ni_scale=0.3_p_ !dimension-less
contains

pure function te_func(radcor,ns) result (z) !unit: kev
    use constants,only: one
    use gk_module,only: te0,kappa_te
    use magnetic_field,only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
        integer,intent(in) :: ns
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=te0(ns)*exp(-kappa_te(ns)*a*te_scale*tanh((r-r0)/(a*te_scale)))
  end function te_func

 function ti_func(radcor) result (z) !in the unit of kev
    use constants,only: one
    use fk_module,only: ti0,kappa_ti
    use magnetic_field,only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=ti0*exp(-kappa_ti*a*ti_scale*tanh((r-r0)/(a*ti_scale)))
  end function ti_func


pure function ne_func(radcor,ns) result (z) !in the SI unit
    use constants,only: one
    use gk_module,only: ne0,kappa_ne
    use magnetic_field,only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
        integer,intent(in) :: ns
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=ne0(ns)*exp(-kappa_ne(ns)*a*ne_scale*tanh((r-r0)/(a*ne_scale)))
  end function ne_func

  
function ni_func(radcor) result (z) !in the SI unit
    use constants,only: one
    use fk_module,only: ni0,kappa_ni
    use magnetic_field,only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=ni0*exp(-kappa_ni*a*ni_scale*tanh((r-r0)/(a*ni_scale)))
  end function ni_func


 function kappa_te_func(radcor,ns) result (z)
    use constants,only: one
    use gk_module,only: kappa_te
    use magnetic_field,only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
        integer,intent(in) :: ns
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=kappa_te(ns)*(one-tanh((r-r0)/(a*te_scale))**2)
  end function kappa_te_func

 function kappa_ti_func(radcor) result (z)
    use constants,only: one
    use fk_module,only: kappa_ti
    use magnetic_field,only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=kappa_ti*(one-tanh((r-r0)/(a*ti_scale))**2)
  end function kappa_ti_func

  function kappa_ne_func(radcor,ns) result (z)
    use constants,only: one
    use gk_module,only : kappa_ne
    use magnetic_field, only : minor_r_radcor !function
    implicit none
    real(p_),intent(in) :: radcor
    integer,intent(in) :: ns
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=kappa_ne(ns)*(one-tanh((r-r0)/(a*ne_scale))**2)
  end function kappa_ne_func

  function kappa_ni_func(radcor) result (z)
    use constants,only: one
    use fk_module,only : kappa_ni
    use magnetic_field, only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=kappa_ni*(one-tanh((r-r0)/(a*ni_scale))**2)
  end function kappa_ni_func

end module density_temperature_profile_mod
