module load_gk_mod
contains
  subroutine load_gk(ns, nm, radcor_e, theta_e, alpha_e, vpar_e,mu_e,v_e, ptcl_num0_e, touch_bdry_e)
    use constants,only: p_, one,twopi,pi,two,kev,fourpi, zero
    use control_parameters,only: gk_spatial_loading_scheme, & !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
    &  gk_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
    use magnetic_coordinates,only: radcor_low=>radcor_low2,radcor_upp=>radcor_upp2,vol,j_low=>j_low2,j_upp=>j_upp2, & !as input
         & mpol,nflux, theta_1d_array, radcor_1d_array, tor_shift_mc, &
         & jacobian,toroidal_range
    use gk_module,only: total_nm_gk, mass_e, te0, vn_e
    use magnetic_field, only : pfn_func
    use domain_decomposition,only: numprocs,myid
    use interpolate_module
    use math, only : sub_random_yj
    use pputil
    use math,only: shift_to_specified_toroidal_range
    use func_in_mc, only:  br_mc_func, bz_mc_func, bphi_mc_func
    use rnorm_mod, only : rnorm_box_muller_vec
    implicit none
    integer, intent(in) :: ns
    integer,intent(inout) :: nm
    real(p_), intent(out) :: radcor_e(:), theta_e(:), alpha_e(:), vpar_e(:),mu_e(:),v_e(:), ptcl_num0_e(:)
    logical, intent(out) ::  touch_bdry_e(:)
    real(p_) :: ps_vol_e(nm)
    integer :: iseed,next_seed
    integer,parameter :: max_try=10000 !used in rejection method
    real(p_) :: radcor_val,theta_val,rannum1,rannum2, rannum3,tmp
    integer :: i, j, ierr
    real(p_) :: pos1,pos2, pos_max
    real(p_) :: abs_jacobian_func,abs_jacobian_max, jacobian_val
    integer :: np_old, np_new
    real(p_) :: vt,vmin,vmax,v_val, probability_max, maxwellian_max
    real(p_) :: theta_v(nm), phi_v(nm), theta_v0
    real(p_) :: vx(nm), vy(nm), vz(nm), tmp_array(3*nm), vx0,vy0,vz0, v0
    real(p_) :: rp_e(nm), zp_e(nm), phip_e(nm) !particle locations
    real(p_) :: brval, bphival, vr, vphi, angle, bval,bxval,byval,bzval
    real(p_) :: tor_shift_e,   normalizing_factor
    real(p_) :: rg_e(nm), zg_e(nm), phig_e(nm)

    touch_bdry_e=.false. ! initially, all markers are considered as not touching the boundary
    abs_jacobian_max=maxval(abs(jacobian(:,j_low:j_upp)))

    iseed=-(2777+myid*3) !set the iseed in different procs. The random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
    call sub_random_yj(iseed,next_seed,tmp) !just to trigger the use of the iseed, the generated random number is not used, 

    select case(gk_spatial_loading_scheme)
    case(1) !uniform loading in (radcor,theta), my testing indicate this loading give more accurate Monte-Carlo integration than uniform loading in real space
       do i=1,nm     
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          call sub_random_yj(0,next_seed,rannum2) 
          radcor_val = radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
          theta_val = (rannum2-0.5_p_)*twopi
          radcor_e(i) = radcor_val
          theta_e(i) = theta_val
       enddo

    case(2) !uniform in real space
       do i=1,nm
          do j=1,max_try !rejection method to generate nonuniform random numbers
             call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
             call sub_random_yj(0,next_seed,rannum2) 
             radcor_val=radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
             theta_val=(rannum2-0.5_p_)*twopi !scaled to the range [-pi:pi]
             pos1=abs_jacobian_func(theta_val,radcor_val)
             call sub_random_yj(0,next_seed,pos2)
             pos2=pos2*abs_jacobian_max !scaled to the range [0: abs_jacobian_max]
             if(pos1<pos2) then
                cycle
             else
                radcor_e(i)=radcor_val
                theta_e(i)=theta_val
                exit
             endif
          enddo
          if(j.eq.max_try+1) stop "***stop**, rejection method falied"
       enddo
    case default
       stop 'please specify a loading scheme for the spatial distribution of gk markers'
    end select

    do i=1,nm !set toroidal coordinates of markers
       call sub_random_yj(0,next_seed,rannum1)
       alpha_e(i) = toroidal_range*rannum1
    enddo

    do i=1,nm
       call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,tor_shift_mc,theta_e(i),radcor_e(i),tor_shift_e) 
       phip_e(i) = alpha_e(i) + tor_shift_e
    enddo

    vt = sqrt(two*te0(ns)*kev/mass_e(ns)) 

    select case(gk_velocity_loading_scheme)
    case(1) !uniform in velocity using Cartesian coordinates (vx,vy,vz)
       vmin = -3_p_*vt/vn_e(ns) 
       vmax = +3_p_*vt/vn_e(ns)
       do i=1, 3*nm
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          tmp_array(i)=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
       enddo

       do i=1,nm
          vx(i)=tmp_array(i)
          vy(i)=tmp_array(i+nm)
          vz(i)=tmp_array(i+2*nm)
       enddo
       v_e(1:nm)=sqrt(vx(:)*vx(:)+vy(:)*vy(:)+vz(:)*vz(:))
       
    case(2) !maxwellian in veloicty using spherical coordinates
       probability_max = probability_sphere(vt,ns)
       vmax= 10.0_p_*vt/vn_e(ns) !I found I need to make the range wide to enusure that the code does not blowup
       do i=1,nm !set velocity magnitude
          do j=1,max_try !rejection method
             call sub_random_yj(0,next_seed,rannum1) 
             v0 = rannum1*vmax
             pos1 = probability_sphere(v0*vn_e(ns),ns)
             call sub_random_yj(0,next_seed,rannum2)
             pos2 = rannum2*probability_max
             if(pos1>pos2) then !accept this smapling
                v_e(i)= v0
                exit  !move to generating next sampling
             else  
                cycle !reject this sampling, move to the next try
             endif
          enddo
       enddo

       do i=1,nm !set velocity theta
          do j=1,max_try !rejection method to generate nonuniform random numbers
             call sub_random_yj(0,next_seed,rannum1) 
             theta_v0 = rannum1*pi
             pos1 = sin(theta_v0)
             call sub_random_yj(0,next_seed,rannum2)
             pos2 = rannum2*1.0
             if(pos1>pos2) then 
                theta_v(i)= theta_v0
                exit  
             else  
                cycle
             endif
          enddo
       enddo

       do i=1,nm !set velocity phi
          call sub_random_yj(0,next_seed,rannum1)
          phi_v(i)=twopi*rannum1
       enddo

       do i=1,nm !transform veloicty from spherical coordinates to Cartesian coordinates
          vz(i) = v_e(i)*cos(theta_v(i)) 
          vx(i) = v_e(i)*sin(theta_v(i))*cos(phi_v(i))
          vy(i) = v_e(i)*sin(theta_v(i))*sin(phi_v(i))
       end do

    case(3) !maxwellian in veloicty using Cartesian coordinates
       vmin = -3.0_p_*vt/vn_e(ns) 
       vmax = +3.0_p_*vt/vn_e(ns)
       maxwellian_max=1.0
!!$       do i=1,3*nm
!!$          do j=1,max_try !rejection method to generate nonuniform random numbers
!!$             call sub_random_yj(0,next_seed,rannum1)
!!$             v_val = vmin + rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
!!$             pos1=maxwellian_func_electron(v_val*vn_e(ns),ns)
!!$             call sub_random_yj(0,next_seed,pos2) 
!!$             pos2=pos2*maxwellian_max !scaled to [0,maxwellian_max]
!!$             if(pos1<pos2) then
!!$                cycle
!!$             else
!!$                tmp_array(i) = v_val
!!$                exit
!!$             endif
!!$          enddo
!!$       enddo

       tmp_array = rnorm_box_muller_vec(3*nm,zero,sqrt(te0(ns)*kev/mass_e(ns)))/vn_e(ns)

       do i=1,nm
          vx(i)=tmp_array(i)
          vy(i)=tmp_array(i+nm)
          vz(i)=tmp_array(i+2*nm)
       enddo
       v_e(1:nm) = sqrt(vx**2 + vy**2 + vz**2)

    case default
       stop 'please specify a loading scheme for the velocity distribution of electron markers'
    end select

    !if(myid.eq.0) call calculate_possibility_density(vx,nm,100,vmin,vmax)

    do i=1,nm !computing parallel velocity and magnetic moment
       bxval=br_mc_func(theta_e(i),radcor_e(i))*cos(phip_e(i)) +bphi_mc_func(theta_e(i),radcor_e(i))*(-sin(phip_e(i))) !x components in a constant Cartesian coor. system
       byval=br_mc_func(theta_e(i),radcor_e(i))*sin(phip_e(i)) +bphi_mc_func(theta_e(i),radcor_e(i))*cos(phip_e(i)) !y components in a constant Cartesian coor. system
       bzval=bz_mc_func(theta_e(i),radcor_e(i)) !z components in a constant Cartesian coor. system
       bval=sqrt(bxval**2+byval**2+bzval**2) 
       vpar_e(i)=(vx(i)*bxval+vy(i)*byval+vz(i)*bzval)/bval !scalar product
       mu_e(i)=(v_e(i)**2-vpar_e(i)**2)/(two*bval) !normalized magnetic moment
    enddo

    do i=1,nm
       call magnetic_coordinates_to_cylindrical_coordinates(theta_e(i),radcor_e(i),rp_e(i),zp_e(i))
    enddo

    do i=1,nm
       brval=br_mc_func(theta_e(i),radcor_e(i))
       bzval=bz_mc_func(theta_e(i),radcor_e(i))
       bphival=bphi_mc_func(theta_e(i),radcor_e(i))
       angle=atan2(vy(i),vx(i))
       vr = vx(i)*cos(angle)+ vy(i)*sin(angle)
       vphi = -vx(i)*sin(angle)+vy(i)*cos(angle)
       call particle_to_guiding_center_location(rp_e(i),phip_e(i),zp_e(i),vr,vphi,vz(i),brval,bphival,bzval,&
            & rg_e(i),phig_e(i),zg_e(i),ns)
       radcor_e(i) = pfn_func(rg_e(i),zg_e(i)) !calculate the radial coordinate
       call interpolate_from_cylindrical_to_magnetic_coordinates(rg_e(i),zg_e(i),theta_e(i),tor_shift_e)
       alpha_e(i) = phig_e(i) + tor_shift_e
    enddo

    do i=1,nm
       if(radcor_e(i).gt.radcor_upp .or. radcor_e(i) .lt. radcor_low)  touch_bdry_e(i)=.true.
    end do

    !if(myid.eq.0) call calculate_possibility_density(vpar_e,nm,100,minval(vpar_e),maxval(vpar_e))
    !if(myid.eq.0) call calculate_possibility_density(mu_e,nm,100,minval(mu_e),maxval(mu_e))

    if(gk_spatial_loading_scheme.eq.1 .and. gk_velocity_loading_scheme.eq.1) then
       normalizing_factor=total_nm_gk(ns)/(twopi*toroidal_range*(radcor_upp-radcor_low)*(vmax-vmin)**3)
       do i=1,nm
          call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_e(i),radcor_e(i),jacobian_val)
          ps_vol_e(i)=abs(jacobian_val)/(normalizing_factor)
       enddo

    elseif(gk_spatial_loading_scheme.eq.1 .and. (gk_velocity_loading_scheme==2 .or. gk_velocity_loading_scheme==3)) then
       normalizing_factor=total_nm_gk(ns)/(twopi*toroidal_range*(radcor_upp-radcor_low)*(sqrt(pi)*vt/vn_e(ns))**3)
       do i=1,nm
          call linear_2d_interpolation(mpol,nflux,theta_1d_array,radcor_1d_array,jacobian,theta_e(i),radcor_e(i),jacobian_val)
          ps_vol_e(i)=abs(jacobian_val)/(normalizing_factor*maxwellian_func_electron(v_e(i)*vn_e(ns),ns))
       enddo
    elseif(gk_spatial_loading_scheme.eq.2 .and. gk_velocity_loading_scheme.eq.1) then
       normalizing_factor=total_nm_gk(ns)/(vol*(vmax-vmin)**3)
       do i=1,nm
          ps_vol_e(i)=one/(normalizing_factor)
       enddo
    elseif(gk_spatial_loading_scheme.eq.2 .and. (gk_velocity_loading_scheme==2 .or. gk_velocity_loading_scheme==3)) then
       normalizing_factor=total_nm_gk(ns)/(vol*(sqrt(pi)*vt/vn_e(ns))**3) 
       do i=1,nm
          ps_vol_e(i)=one/(normalizing_factor*maxwellian_func_electron(v_e(i)*vn_e(ns),ns))
       enddo
    endif

    !ptcl_num0_e=ne0*sqrt((mass_e/(twopi*te0*kev))**3)*(ln**3*vn_e**3)/normalizing_factor_e !explicitly cancel the exp(-v^2/vt^2) dependence, valid only for the case in which both physical equilibrium distribution and marker distribution are Maxwellian.
    do i=1,nm
       ptcl_num0_e(i) = ps_vol_e(i)*f0(radcor_e(i), v_e(i),ns)
    enddo
    call shift_to_specified_toroidal_range(alpha_e(:)) 
!!$  !assign the loaded particles to the corresponding processors, using the subroutines provided in pputil_yj.f90
    np_old=nm
    call init_pmove(theta_e(:),np_old,twopi,ierr)
    call pmove(theta_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(radcor_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(alpha_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(vpar_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(mu_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(v_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ptcl_num0_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit

    call pmove2(touch_bdry_e(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit

    call end_pmove(ierr)

    nm=np_new

!!$   block
!!$  integer :: file_unit    
!!$  character(5):: filename
!!$  write(filename,'(a1,i4.4)') 'e',myid
!!$  open(newunit=file_unit, file=filename//'.txt')
!!$  do i=1,nm
!!$     write(file_unit,*)  radcor_e(i) ,theta_e(i) ,rg_e(i),zg_e(i),i
!!$  enddo
!!$  close(file_unit)
!!$   end block
!!$  call some_test2(myid,numprocs)
    !        if(myid==0) write(*,*) "w_e(k)=", maxval(ps_vol_e(1:nm)), minval(ps_vol_e(1:nm)), ns, '*****'

  end subroutine load_gk


elemental  function f0(x,v,ns) result (z) ! v in unit of vn, f0 in unit 1/(Ln**3*vn**3)
    use constants,only: p_, two,twopi,kev
    use normalizing,only: Ln
    use gk_module,only: mass_e, vn_e
    use density_temperature_profile_mod,only : te_func, ne_func
    implicit none
    integer, intent(in) :: ns
    real(p_),intent(in) :: x, v
    real(p_) :: z, v_si, te
    v_si=v*vn_e(ns)
    te = te_func(x,ns)
    z = ne_func(x,ns)*sqrt((mass_e(ns)/(twopi*te*kev))**3)*exp(-mass_e(ns)*v_si**2/(two*te*kev))
    z=z*(vn_e(ns)**3*Ln**3)
  end function f0
  
  function maxwellian_func_electron(v,ns) result(z) !v in SI unit
    use constants,only: p_, two,kev
    use gk_module,only: mass_e, te0 !as input
    implicit none
    real(p_):: v,z
    integer, intent(in) :: ns
    z=exp(-mass_e(ns)*v*v/(two*te0(ns)*kev))
  end function maxwellian_func_electron

  function probability_sphere(v,ns) result(z) !v in SI unit
    use constants,only:p_
    use constants,only:two,kev
    use gk_module,only: mass_e,te0 !as input
    implicit none
    real(p_) :: v,z
    integer, intent(in) :: ns
    z = v**2*exp(-mass_e(ns)*v**2/(two*te0(ns)*kev))
  end function probability_sphere

  function probability(vx,vy,vz,ns) result(z) !v in SI unit
    use constants,only:p_
    use constants,only:two,kev
    use gk_module,only: mass_e,te0 !as input
    implicit none
    real(p_):: vx,vy,vz,z
    integer, intent(in) :: ns
    z=exp(-mass_e(ns)*(vx**2+vy**2+vz**2)/(two*te0(ns)*kev))
  end function probability


  subroutine particle_to_guiding_center_location(r,phi,z,vr,vphi,vz,brval,bphival,bzval,rg,phig,zg,ns)
    use constants,only:p_
    use constants,only: one,two,half_pi
    use gk_module,only:mass_e, charge_e
    implicit none
    integer, intent(in) :: ns
    real(p_),intent(in) :: r, phi, z, vr, vphi, vz, brval, bphival, bzval
    real(p_),intent(out) :: rg, phig, zg
    real(p_):: v,bval

    v=sqrt(vr**2+vz**2+vphi**2)
    bval=sqrt(brval**2+bphival**2+bzval**2)

    rg=sqrt((r+mass_e(ns)/(bval**2*charge_e(ns))*(vphi*bzval-vz*bphival))**2+&
         &(mass_e(ns)/(bval**2*charge_e(ns))*(vz*brval-vr*bzval))**2)
    phig=phi+asin(mass_e(ns)/(bval**2*charge_e(ns))*(-vr*bzval+vz*brval)/rg)
    !zg=z+mass/(bval**2*charge)*(-vr*bphival-vphi*brval) !wrong, pointed out by Yingfeng Xu
    zg=z+mass_e(ns)/(bval**2*charge_e(ns))*(vr*bphival-vphi*brval) !corrected.
  end subroutine particle_to_guiding_center_location
  
end module load_gk_mod
