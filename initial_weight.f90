module delta_ne_mod
  use constants,only: p_
  use constants,only: twopi,two,pi
  use gk_module,only:ne0 !1/m^3
  use magnetic_coordinates,only: radcor_low2,radcor_upp2, ntor=>nsegment
  use control_parameters, only : nh
    use magnetic_field, only : qfunc

  implicit none
  private
  public delta_ne,delta_ne_theta,delta_ne_psi,delta_ne_alpha
  real(p_):: kradcor0,delta_ne0
  real(p_):: radcor_center,radcor_width
  real(p_)::mprime
  real(p_):: q, qmax,qmin
  integer:: integer_val


contains
  function delta_ne(radcor,theta,alpha) result (z) ! 1/m^3
     use math, only : random_yj
    real(p_)::radcor,theta,alpha,z
    integer:: i
    delta_ne0=ne0(1)*1.d-4 !amplitude of the initial perturbation in terms of the equilibrium density
    !kradcor0=twopi/(radcor_upp2-radcor_low2) !radial wave number
    kradcor0=pi/(radcor_upp2-radcor_low2) !the fundament radial wave number in sine expansion
    radcor_center=(radcor_upp2+radcor_low2)/two
    radcor_width=(radcor_upp2-radcor_low2)*0.25_p_
    q=qfunc(radcor)
    !mprime=mod(ntor*q,1.0_p_) !this value is chosen in order to satisfy perioidic condition on (theta,alpha) plane when using the following spatial distribution
    qmax=qfunc(radcor_upp2)
    qmin=qfunc(radcor_low2)
    !integer_val=nint(ntor*(qmax+qmin)/2)
    integer_val=ntor*nh*nint((qmax+qmin)/2)
    mprime=ntor*nh*q-integer_val

    !  write(*,*) 'mprime=',mprime
    !  if(abs(mprime)>0.5) mprime=mprime-sign(1._p_,mprime) !choose the smallest abs(mprime) that is possible, can this help solve the radial resolution problem? to be verified, seems unlikely.
    !  write(*,*) 'mprime=',mprime, 'after'

    !z=delta_ne0 !for testing
    !z=delta_ne0*cos(kradcor0*radcor+mprime*theta+ntor*alpha) !incorrect
    !z=delta_ne0*sin(kradcor0*(radcor-radcor_low2))*cos(mprime*theta+ntor*alpha) !for ITG instabilities
    !z=delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*cos(mprime*theta+ntor*alpha) !set exponential radial dependence so that the perturbation near the radial boundary is small
    !z=delta_ne0*cos(ntor*alpha) !for testing
    !z=delta_ne0 !for testing
    z=0._p_
    do i=0,nh !set exponential radial dependence so that the perturbation near the radial boundary is small
       z = z + delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2) &
       &   *(cos(mprime*theta+i*ntor*alpha)+sin(mprime*theta+i*ntor*alpha))
    enddo
    z= delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*(random_yj(0)-0.5)
  end function delta_ne

  function delta_ne_theta_orgin(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    delta_ne0=ne0(1)*0.01_p_ !amplitude of the initial perturbation is chosen 0.05 of the equilibrium density
    kradcor0=twopi/(radcor_upp2-radcor_low2) !radial wave number
    radcor_center=(radcor_upp2+radcor_low2)/two
    radcor_width=(radcor_upp2-radcor_low2)*0.1_p_
    q=qfunc(radcor)
    !mprime=mod(ntor*q,1.0_p_) !this value is chosen in order to satisfy perioidic condition on (theta,alpha) plane when using the following spatial distribution
    qmax=qfunc(radcor_upp2)
    qmin=qfunc(radcor_low2)
    integer_val=nint(ntor*(qmax+qmin)/2)
    mprime=ntor*q-integer_val

    z=-delta_ne0*sin(kradcor0*(radcor-radcor_low2))*sin(mprime*theta+ntor*alpha)*mprime
    !z=-delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*sin(mprime*theta+ntor*alpha)*mprime
  end function delta_ne_theta_orgin


  function delta_ne_theta(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    real(p_):: del_theta
    del_theta=twopi/8._p_
    z=(delta_ne(radcor,theta+del_theta,alpha)-delta_ne(radcor,theta-del_theta,alpha))/(2*del_theta) !using numerical difference
  end function delta_ne_theta

  function delta_ne_psi(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    real(p_):: del_psi

    !z=kradcor0*delta_ne0*cos(kradcor0*(radcor-radcor_low2))*cos(mprime*theta+ntor*alpha) !wrong, there is a radial dependence in mprime
    del_psi=0.01*(radcor_upp2-radcor_low2)
    z=(delta_ne(radcor+del_psi,theta,alpha)-delta_ne(radcor-del_psi,theta,alpha))/(2*del_psi) !using numerical difference
  end function delta_ne_psi

  function delta_ne_alpha(radcor,theta,alpha) result (z) ! temporary, for testing
    real(p_)::radcor,theta,alpha,z
    delta_ne0=ne0(1)*0.01_p_ !amplitude of the initial perturbation is chosen 0.05 of the equilibrium density
    kradcor0=twopi/(radcor_upp2-radcor_low2) !radial wave number
    radcor_center=(radcor_upp2+radcor_low2)/two
    radcor_width=(radcor_upp2-radcor_low2)*0.1_p_
    q=qfunc(radcor)
    !mprime=mod(ntor*q,1.0_p_) !this value is chosen in order to satisfy perioidic condition on (theta,alpha) plane when using the following spatial distribution
    qmax=qfunc(radcor_upp2)
    qmin=qfunc(radcor_low2)
    integer_val=nint(ntor*(qmax+qmin)/2)
    mprime=ntor*q-integer_val

    z=-delta_ne0*sin(kradcor0*(radcor-radcor_low2))*sin(mprime*theta+ntor*alpha)*ntor
    !z=-delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*sin(mprime*theta+ntor*alpha)*ntor
  end function delta_ne_alpha

  function delta_ne_alpha_tmp(radcor,theta,alpha) result (z) ! temporary, for testing
    use magnetic_coordinates,only:toroidal_range
    real(p_)::radcor,theta,alpha,z
    real(p_):: del_alpha
    del_alpha=toroidal_range*0.01_p_
    !z=-delta_ne0*sin(kradcor0*(radcor-radcor_low2))*sin(mprime*theta+ntor*alpha)*ntor
    !z=-delta_ne0*exp(-(radcor-radcor_center)**2/radcor_width**2)*sin(mprime*theta+ntor*alpha)*ntor
    z=(delta_ne(radcor,theta,alpha+del_alpha)-delta_ne(radcor,theta,alpha-del_alpha))/(2*del_alpha) !using numerical difference
  end function delta_ne_alpha_tmp
end module delta_ne_mod



module initialization_mod
contains

  subroutine initialize_weight_i(w_i)
    use constants,only:p_
    use fk_module,only: nmarker_i, ps_vol_i 
    use fk_module,only: radcor_i,theta_i,alpha_i,v_i,ni0
    use domain_decomposition,only:myid
    implicit none
    real(p_),intent(out):: w_i(nmarker_i)
    integer:: k
    real(p_):: rannum,tmp
    integer:: iseed,next_seed
    ! real(p_),parameter:: eps=1.0d20
    !  real(p_):: initial_delta_f_i !function name
!!$  do k=1,nmarker_i
!!$     w_i(k)=ps_vol_i(k)*eps !initial weight of markers
!!$  enddo  

    do k=1,nmarker_i
       w_i(k)=ps_vol_i(k)*initial_delta_f_i(radcor_i(k),theta_i(k),alpha_i(k),v_i(k)) !initial weight of markers
       !write(*,*) 'w_i(k)=',w_i(k),ps_vol_i(k),initial_delta_f_i(radcor_i(k),theta_i(k),alpha_i(k),v_i(k))
    enddo

!!$  iseed=-(3777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
!!$  ! now generate the random numbers
!!$  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 
!!$  do k=1,nmarker_i
!!$     call sub_random_yj(0,next_seed,rannum) !0 means using last random number as iseed 
!!$     !w_i(k)=ni0*1.0d-6*(rannum-0.5_p_)*2._p_ !random, for testing
!!$w_i(k)=ps_vol_i(k)*initial_delta_f_i(radcor_i(k),theta_i(k),alpha_i(k),v_i(k))*(rannum-0.5_p_)*2._p_ !initial weight of markers
!!$  enddo


!!$if(myid.eq.1) then
!!$   do k=1,nmarker_i
!!$      write(*,*) alpha_i(k),initial_delta_f_i(radcor_i(k),theta_i(k),alpha_i(k),v_i(k)),w_i(k)
!!$   enddo
!!$endif

  end subroutine initialize_weight_i

  function initial_delta_f_i(radcor,theta,alpha,v) result (z) ! for tesing Monte-Carlo integration, v in unit of vn, z in unit 1/(Ln**3*vn**3)
    use constants,only: p_
    use constants,only: two,twopi,kev
    use normalizing,only: vn_i,Ln
    use fk_module,only: mass_i,ti0
    use delta_ne_mod,only: delta_ne !function name
    implicit none
    real(p_):: radcor,theta,alpha,v,z
    real(p_):: v_si

    v_si=v*vn_i
    z=delta_ne(radcor,theta,alpha)*sqrt((mass_i/(twopi*ti0*kev))**3)*exp(-mass_i*v_si**2/(two*ti0*kev))
    z=z*(vn_i**3*Ln**3)
  end function initial_delta_f_i


  subroutine initialize_gk_weight()
    use constants,only: p_, two
    use gk_module,only: nsm, nm_gk, vn_e,ne0, ptcl_num0_e, &
                       mu_e,v_e,vpar_e,radcor_e,theta_e,alpha_e, w_e !output
    use magnetic_coordinates,only:mpol,nflux,radcor_1d_array,theta_1d_array, nsegment
    use control_parameters, only : nh
    use table_in_mc,only: b_mc
    use math, only : sub_random_yj
    use interpolate_module,only: linear_2d_interpolation
    use domain_decomposition,only:myid
    implicit none
    integer:: k,ns
    !real(p_),parameter:: eps=1.0d20
    real(p_):: v,b_val
    !real(p_):: fe2_arbitrary
    real(p_):: rannum,tmp
    integer:: iseed,next_seed


!!$  do k=1,nm_gk
!!$     w_e(k)=ps_vol_e(k)*eps !initial weight of markers
!!$  enddo
  iseed=-(1777+myid*2) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 
do ns=1,nsm
    do k=1,nm_gk(ns)
       !call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,b_mc,theta_e(k),radcor_e(k),b_val)
       !v=sqrt(two*mu_e(k)*b_val*2+vpar_e(k)**2) !in unit of vn_e, wrong! a bug found when tested the deposition procedure
       ! v=sqrt(mu_e(k)*two*b_val+vpar_e(k)**2) !in unit of vn_e
       !write(*,*) v_e(k),v
       !     w_e(k)=ps_vol_e(k)*(mu_e(k)*b_val*2+vpar_e(k)**2)/3.0*fe2_arbitrary(radcor_e(k),theta_e(k),alpha_e(k),v) !initial weight of markers, for testing
       !w_e(k,ns)=ps_vol_e(k,ns)*initial_delta_f_e(radcor_e(k,ns),theta_e(k,ns),alpha_e(k,ns),v_e(k,ns),vpar_e(k,ns),ns) !initial weight of markers
       w_e(k,ns)=ptcl_num0_e(k,ns)*0.001*cos(alpha_e(k,ns)*nsegment*nh)
       !w_e(k,ns)=0.
       !call sub_random_yj(0,next_seed,rannum) !0 means using last random number as iseed        
       !w_e(k,ns)=ps_vol_e(k,ns)*ne0(ns)*10**(-6)*(rannum-0.5)
       !initial weight of markers
       !write(*,*) 'w_e=',w_e(k),ps_vol_e(k),initial_delta_f_e(radcor_e(k),theta_e(k),alpha_e(k),v_e(k))
    enddo
 enddo


!!$ iseed=-(5777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
!!$  ! now generate the random numbers
!!$  call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 
!!$  do k=1,nm_gk
!!$     call sub_random_yj(0,next_seed,rannum) !0 means using last random number as iseed 
!!$     w_e(k)=ne0*1.0d-10*(rannum-0.5_p_)*2._p_ !random, for testing
!!$  enddo

end subroutine initialize_gk_weight

  function initial_delta_f_e(radcor,theta,alpha,v,vpar,ns) result (z) ! v in unit of vn_e, z in unit 1/(Ln**3*vn**3)
    use constants,only: p_
    use constants,only: two,twopi,kev
    use normalizing,only: Ln
    use gk_module,only: mass_e,te0, vn_e
    use delta_ne_mod,only: delta_ne
    implicit none
    integer,intent(in) :: ns
    real(p_),intent(in) :: radcor,theta,alpha,v,vpar
    real(p_) :: z
    real(p_),parameter:: vpar_drift=0. !in unit of vn_e
    real(p_):: vper_sq

    vper_sq=v**2-vpar**2
    z=delta_ne(radcor,theta,alpha)*sqrt((mass_e(ns)/(twopi*te0(ns)*kev))**3) &
         & *exp(-mass_e(ns)*(vper_sq*vn_e(ns)**2+(vpar-vpar_drift)**2*vn_e(ns)**2)/(two*te0(ns)*kev))
    z=z*(vn_e(ns)**3*Ln**3)
  end function initial_delta_f_e


end module initialization_mod
