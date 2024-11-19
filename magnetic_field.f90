module magnetic_field
contains
  !all independent and dependent variables are in SI units
  function psi_func(x,z) result(f)!SI units, poloidal magnetic flux function
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z,f
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi,x,z,f)  
  end function psi_func

  function pfn_func(x,z) result(f)!normalized poloidal magnetic flux
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi !as input
    use radial_module,only:psi_axis,psi_lcfs
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z,f
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi,x,z,f)  
    f=(f-psi_axis)/(psi_lcfs-psi_axis)
  end function pfn_func

  function psi_r_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_x
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z, psi_r_func
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi_x,x,z,psi_r_func)  
  end function psi_r_func


  function psi_z_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_z !as input
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z,psi_z_func
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi_z,x,z,psi_z_func)  
  end function psi_z_func


  function psi_rr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xx !as input
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z,psi_rr_func
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi_xx,x,z,psi_rr_func)  
  end function psi_rr_func

  function psi_zz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zz !as input
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z,psi_zz_func
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi_zz,x,z,psi_zz_func)  
  end function psi_zz_func


  function psi_rz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xz 
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z,psi_rz_func
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi_xz,x,z,psi_rz_func)  
  end function psi_rz_func


  function psi_zr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zx 
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: x,z,psi_zr_func
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi_zx,x,z,psi_zr_func)  
  end function psi_zr_func

  function psi_gradient_func(xval,zval) result (f)
    use constants,only:p_
    use poloidal_flux_2d,only: xarray,zarray,psi_gradient,y2a_gradient,nx,nz
    use interpolate_module,only: linear_2d_interpolation
    implicit none
    real(p_):: f,xval,zval
    call linear_2d_interpolation(nx,nz,xarray,zarray,psi_gradient,xval,zval,f)
  end function psi_gradient_func


  function br(r,z) !R component of magnetic field
    use constants,only:p_
    implicit none
    real(p_):: br,r,z
    br=-psi_z_func(r,z)/r
  end function br

  function bz(r,z) !Z component of magnetic field
    use constants,only:p_
    implicit none
    real(p_):: bz,r,z
    bz=psi_r_func(r,z)/r
  end function bz

  function bphi(r,z) !phi component of magnetic field
    use constants,only:p_
    implicit none
    real(p_):: bphi,r,z
    bphi=g_func(psi_func(r,z))/r
  end function bphi


  function b(r,z) !strength of magnetic field
    use constants,only:p_
    implicit none
    real(p_):: b,r,z
    b=sqrt(br(r,z)**2+bz(r,z)**2+bphi(r,z)**2)
  end function b

  function br_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_):: br_r,r,z
    br_r=psi_z_func(r,z)/r**2-psi_rz_func(r,z)/r
  end function br_r


  function br_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_):: br_z,r,z
    br_z=-psi_zz_func(r,z)/r
  end function br_z

  function bz_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_):: bz_r,r,z
    bz_r=-psi_r_func(r,z)/r**2+psi_rr_func(r,z)/r
  end function bz_r


  function bz_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_):: bz_z,r,z
    bz_z=psi_rz_func(r,z)/r
  end function bz_z


  function bphi_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_):: bphi_r,r,z
    bphi_r=gprime(psi_func(r,z))*psi_r_func(r,z)/r-g_func(psi_func(r,z))/r**2
  end function bphi_r

  function bphi_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one,zero
    implicit none
    real(p_):: bphi_z,r,z
    bphi_z=gprime(psi_func(r,z))*psi_z_func(r,z)/r
  end function bphi_z

  function b_r(r,z) !partial derivative of  magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_):: b_r,r,z
    b_r=-b(r,z)/r+one/(b(r,z)*r*r)*(psi_r_func(r,z)*psi_rr_func(r,z) &
         & +psi_z_func(r,z)*psi_rz_func(r,z)+g_func(psi_func(r,z))*g_r(r,z))
  end function b_r

  function b_z(r,z) !partial derivative of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_):: b_z,r,z
    b_z=one/(b(r,z)*r*r)*(psi_r_func(r,z)*psi_rz_func(r,z)+&
         & psi_z_func(r,z)*psi_zz_func(r,z)+g_func(psi_func(r,z))*g_z(r,z))
  end function b_z

  function b_phi(r,z) result(f)
    use constants,only:p_
    implicit none
    real(p_):: f,r,z
    f=0._p_
  end function b_phi

  function unitbr_r(r,z)  result(f)!partial derivative of the component of magnetic unit vector
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_):: f,r,z
    real(p_):: bval
    bval=b(r,z)
    f=(br_r(r,z)*bval-b_r(r,z)*br(r,z))/bval**2
  end function unitbr_r

  function unitbr_z(r,z) result(f)
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_):: f,r,z
    real(p_):: bval

    bval=b(r,z)
    f=(br_z(r,z)*bval-b_z(r,z)*br(r,z))/bval**2
  end function unitbr_z


  function unitbr_phi(r,z) result(f)
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_):: f,r,z
    real(p_):: bval

    !bval=b(r,z)
    !f=(0*bval-b_phi(r,z)*br(r,z))/bval**2
    f=0.
  end function unitbr_phi


  function unitbz_r(r,z) result(f) 
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_):: f,r,z
    real(p_)::bval
    bval=b(r,z)
    f=(bz_r(r,z)*bval-b_r(r,z)*bz(r,z))/bval**2
  end function unitbz_r

  function unitbz_z(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_):: f,r,z
    real(p_):: bval

    bval=b(r,z)
    f=(bz_z(r,z)*bval-b_z(r,z)*bz(r,z))/bval**2
  end function unitbz_z

  function unitbz_phi(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_):: f,r,z
    real(p_):: bval
    !bval=b(r,z)
    !f=(0*bval-b_phi(r,z)*bz(r,z))/bval**2
    f=0.
  end function unitbz_phi


  function unitbphi_r(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_):: f,r,z
    real(p_):: bval
    bval=b(r,z)
    f=(bphi_r(r,z)*bval-b_r(r,z)*bphi(r,z))/bval**2
  end function unitbphi_r


  function unitbphi_z(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_):: f,r,z
    real(p_):: bval
    bval=b(r,z)
    f=(bphi_z(r,z)*bval-b_z(r,z)*bphi(r,z))/bval**2
  end function unitbphi_z


  function g_r(r,z)
    use constants,only:p_
    implicit none
    real(p_):: r,z,g_r
    real(p_):: psival
    psival=psi_func(r,z)
    g_r=gprime(psival)*psi_r_func(r,z)
  end function g_r

  function g_z(r,z)
    use constants,only:p_
    implicit none
    real(p_):: r,z,g_z
    real(p_):: psival
    psival=psi_func(r,z)
    g_z=gprime(psival)*psi_z_func(r,z)
  end function g_z

  function g_func(psival) result(f) !all quantities are in S.I units.
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fpsi
    use interpolate_module,only: linear_1d_interpolation
    implicit none
    real(p_):: f,psival
    call linear_1d_interpolation(npsi,psi_1d,fpsi,psival,f)  
  end function g_func

  function gprime(psival) result(f)
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fprime
    use interpolate_module,only: linear_1d_interpolation
    implicit none
    real(p_):: f,psival
    call linear_1d_interpolation(npsi,psi_1d,fprime,psival,f)  
  end function gprime

pure  function qfunc(pfn) result(f) !safety factor, with correct sign
    use constants,only:p_
    use radial_module,only:npsi,pfn_npsi,q_with_sign
    use interpolate_module,only: linear_1d_interpolation
    implicit none
    real(p_),intent(in) :: pfn
    real(p_) :: f
    call linear_1d_interpolation(npsi,pfn_npsi,q_with_sign,pfn,f)  
  end function qfunc

  function qfunc0(psival) result(f) !safety factor
    !not used in computing the orbits, only used to do analytical estimation of some quantities, such as bounce frequency
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,qpsi
    use interpolate_module,only: linear_1d_interpolation
    implicit none
    real(p_):: f,psival
    call linear_1d_interpolation(npsi,psi_1d,qpsi,psival,f)  
  end function qfunc0


  function radcor_as_func_of_pfn(pfn) result (z)
  use constants,only: p_
  implicit none
  real(p_):: pfn,z
  z=pfn
end function radcor_as_func_of_pfn


pure function minor_r_radcor(radcor) result (z)
  use constants,only: p_
  use constants,only: two,twopi
  use magnetic_coordinates,only: nflux,radcor_1d_array,minor_r_array
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_),intent(in):: radcor
  real(p_) :: z

  call linear_1d_interpolation(nflux,radcor_1d_array,minor_r_array,radcor,z)  
end function minor_r_radcor

 function radcor_minor_r(minor_r) result (z)
  use constants,only: p_
  use magnetic_coordinates,only: nflux,radcor_1d_array,minor_r_array
  use interpolate_module,only: linear_1d_interpolation_general_case
  implicit none
  real(p_):: minor_r,z

  call linear_1d_interpolation_general_case(nflux,minor_r_array,radcor_1d_array,minor_r,z)

end function radcor_minor_r


function minor_r_prime(radcor) result (z) !derivative of minor_r with respect to the radial coordinate
  use constants,only: p_
  use constants,only: two
  use magnetic_coordinates,only: nflux,radcor_1d_array,minor_r_prime_array
  use interpolate_module,only: linear_1d_interpolation
  implicit none
  real(p_) :: radcor, z

  call linear_1d_interpolation(nflux, radcor_1d_array, minor_r_prime_array, radcor, z)  
end function minor_r_prime
  
end module magnetic_field
