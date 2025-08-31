module load_gk_mod
contains
  subroutine load_gk(ns, nm, xgc, zgc, ygc, vpar_gk,mu_gk,v_gk, ps_vol_gk, ptcl_num0_gk, touch_bdry_gc)
    use constants,only: p_, one,twopi,pi,two,kev,fourpi, zero

    use magnetic_coordinates,only: xlow, xupp,vol, & !as input
         & mpol,nrad, zgrid, xgrid, tor_shift_mc, &
         & jacobian,toroidal_range, toroidal_range
    use gk_module,only: total_nm_gk, mass_gk, tgk0, &
         & vn_gk, gk_spatial_loading_scheme, gk_velocity_loading_scheme 
    use magnetic_field, only : pfn_func
    use gk_profile_funcs,only : gkn_func
    use domain_decomposition,only: numprocs, myid
    use misc, only: magnetic_coordinates_to_cylindrical_coordinates
    use interpolate_module
    use math, only: sub_random_yj
    use misc, only: abs_jacobian_func
    use pputil
    use math,only: shift_toroidal
    use func_in_mc, only:  br_mc_func, bz_mc_func, bphi_mc_func
    use map_to_mc, only: interpolate_from_cylindrical_to_magnetic_coordinates
    implicit none
    integer, intent(in) :: ns
    integer,intent(inout) :: nm
    real(p_), intent(out) :: xgc(:), zgc(:), ygc(:), vpar_gk(:),mu_gk(:),v_gk(:), ps_vol_gk(:), ptcl_num0_gk(:)
    logical, intent(out) ::  touch_bdry_gc(:)

    integer :: iseed,next_seed
    integer, parameter :: max_try=10000 !used in rejection method
    real(p_) :: radcor_val,theta_val,rannum1,rannum2, rannum3, z1, z2, tmp, scalar
    integer :: i, j, ierr, nmax
    real(p_) :: pos1,pos2, pos_max
    real(p_) :: abs_jacobian_max, jacobian_val
    integer :: np_old, np_new
    real(p_) :: vt,vmin,vmax,v_val, probability_max, maxwellian_max
    real(p_) :: theta_v(nm), phi_v(nm), theta_v0
    real(p_) :: vx(nm), vy(nm), vz(nm), tmp_array(3*nm), vx0,vy0,vz0, v0
    real(p_) :: rp_e(nm), zp_e(nm), phip_e(nm) !particle locations
    real(p_) :: brval, bphival, vr, vphi, angle, bval,bxval,byval,bzval
    real(p_) :: tor_shift_e,   normalizing_factor, t(mpol, 1:nrad), tmax
    real(p_) :: rg_e(nm), zg_e(nm), phig_e(nm)

    touch_bdry_gc = .false. ! initially, all markers are considered as not touching the boundary
    abs_jacobian_max = maxval(abs(jacobian(:,1:nrad)))

    do i =1, mpol
       do j = 1, nrad
          t(i,j) = abs(jacobian(i,j))*gkn_func(xgrid(j), ns)
       enddo
    enddo
    tmax = maxval(t)

    !set different seeds in different procs, so that random numbers in different procs are independent from each other.
    iseed = 27777 + myid*300
    call sub_random_yj(iseed,next_seed,tmp) !to trigger the use of the iseed, the generated random number is not used, 

    select case(gk_spatial_loading_scheme)
    case(1) !uniform loading in (radcor,theta) my testing indicate this loading give more accurate Monte-Carlo integration than uniform loading in real space
       do i=1,nm     
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          call sub_random_yj(0,next_seed,rannum2) 
          radcor_val = xlow+(xupp-xlow)*rannum1 !scale the random number to the range [xlow: xupp]
          theta_val = (rannum2-0.5_p_)*twopi
          xgc(i) = radcor_val
          zgc(i) = theta_val
       enddo

    case(2) !uniform in real space
       do i=1,nm
          do j=1,max_try !rejection method to generate nonuniform random numbers
             call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
             call sub_random_yj(0,next_seed,rannum2) 
             radcor_val = xlow +(xupp-xlow)*rannum1 !scaled to the range [xlow: xupp]
             theta_val = (rannum2-0.5_p_)*twopi !scaled to the range [-pi:pi]
             pos1 = abs_jacobian_func(theta_val, radcor_val)
             call sub_random_yj(0, next_seed, pos2)
             pos2=pos2*abs_jacobian_max !scaled to the range [0: abs_jacobian_max]
             if(pos1<pos2) then
                cycle
             else
                xgc(i)=radcor_val
                zgc(i)=theta_val
                exit
             endif
          enddo
          if(j.eq.max_try+1) stop "***stop**, rejection method falied"
       enddo
    case(3) !nonuniform in real space,  proportional to the equilibrium density n0(x)
       do i=1,nm
          do j=1,max_try 
             call sub_random_yj(0,next_seed,rannum1)
             call sub_random_yj(0,next_seed,rannum2) 
             radcor_val = xlow +(xupp-xlow)*rannum1 
             theta_val = (rannum2-0.5_p_)*twopi 
             pos1 = abs_jacobian_func(theta_val, radcor_val)*gkn_func(radcor_val, ns)
             call sub_random_yj(0, next_seed, pos2)
             pos2=pos2*tmax
             if(pos1<pos2) then
                cycle
             else
                xgc(i)=radcor_val
                zgc(i)=theta_val
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
       ygc(i) = toroidal_range*rannum1
    enddo

    do i=1,nm
       call linear_2d_interpolate(mpol, nrad, zgrid, xgrid,tor_shift_mc,zgc(i),xgc(i),tor_shift_e) 
       phip_e(i) = ygc(i) + tor_shift_e
    enddo

    vt = sqrt(two*tgk0(ns)*kev/mass_gk(ns)) 

    select case(gk_velocity_loading_scheme)
    case(1) !uniform in velocity using Cartesian coordinates (vx, vy, vz)
       vmin = -3_p_*vt/vn_gk(ns) 
       vmax = +3_p_*vt/vn_gk(ns)
       do i=1, 3*nm
          call sub_random_yj(0,next_seed,rannum1) !0 means using last random number as iseed 
          tmp_array(i)=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
       enddo

       do i=1,nm
          vx(i) = tmp_array(i)
          vy(i) = tmp_array(i+nm)
          vz(i) = tmp_array(i+2*nm)
       enddo
       v_gk(1:nm)=sqrt(vx(:)*vx(:)+vy(:)*vy(:)+vz(:)*vz(:))

    case(2) !maxwellian in veloicty using spherical coordinates
       probability_max = probability_sphere(vt, ns)
       vmax= 6.0_p_*vt/vn_gk(ns) !I found the range is important in enusuring the code does not blowup
       do i =1, nm !set velocity magnitude
          do j = 1, max_try !rejection method
             call sub_random_yj(0, next_seed, rannum1) 
             v0 = rannum1*vmax
             pos1 = probability_sphere(v0*vn_gk(ns), ns)
             call sub_random_yj(0,next_seed,rannum2)
             pos2 = rannum2*probability_max
             if(pos1>pos2) then !accept this smapling
                v_gk(i) = v0
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
          vz(i) = v_gk(i)*cos(theta_v(i)) 
          vx(i) = v_gk(i)*sin(theta_v(i))*cos(phi_v(i))
          vy(i) = v_gk(i)*sin(theta_v(i))*sin(phi_v(i))
       end do

    case(3) !maxwellian in veloicty using Cartesian coordinates
       vmin = -4.0_p_*vt/vn_gk(ns) 
       vmax = +4.0_p_*vt/vn_gk(ns)
       maxwellian_max = 1.0
!!$       do i = 1, 3*nm
!!$          do j = 1, max_try !rejection method to generate nonuniform random numbers
!!$             call sub_random_yj(0, next_seed, rannum1)
!!$             v_val = vmin + rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
!!$             pos1 = maxwellian0(v_val*vn_gk(ns), ns)
!!$             call sub_random_yj(0, next_seed, pos2) 
!!$             pos2 = pos2*maxwellian_max !scaled to [0,maxwellian_max]
!!$             if(pos1 < pos2) then
!!$                cycle
!!$             else
!!$                tmp_array(i) = v_val
!!$                exit
!!$             endif
!!$          enddo
!!$       enddo
!!$
       scalar = sqrt(tgk0(ns)*kev/mass_gk(ns))/vn_gk(ns)
       do i = 1, 3*nm, 2
          ! Box-Muller method (DOI: 10.1214/aoms/1177706645):
          call sub_random_yj(0, next_seed, rannum1)
          call sub_random_yj(0, next_seed, rannum2)
          z1 = sqrt(-two * log(rannum1)) * cos(twopi * rannum2)
          z2 = sqrt(-two * log(rannum1)) * sin(twopi * rannum2)
          tmp_array(i) = z1*scalar
          if(i+1 > 3*nm) exit
          tmp_array(i+1) = z2*scalar
       enddo

       do i=1,nm
          vx(i) = tmp_array(i)
          vy(i) = tmp_array(i+nm)
          vz(i) = tmp_array(i+2*nm)
       enddo
       v_gk(1:nm) = sqrt(vx**2 + vy**2 + vz**2)

    case default
       stop 'please specify a loading scheme for the velocity distribution of electron markers'
    end select

    !if(myid.eq.0) call calculate_possibility_density(vx,nm,100,vmin,vmax)

    do i=1,nm !computing parallel velocity and magnetic moment
       bxval=br_mc_func(zgc(i),xgc(i))*cos(phip_e(i)) +bphi_mc_func(zgc(i),xgc(i))*(-sin(phip_e(i))) !x components in a constant Cartesian coor. system
       byval=br_mc_func(zgc(i),xgc(i))*sin(phip_e(i)) +bphi_mc_func(zgc(i),xgc(i))*cos(phip_e(i)) !y components in a constant Cartesian coor. system
       bzval=bz_mc_func(zgc(i),xgc(i)) !z components in a constant Cartesian coor. system
       bval=sqrt(bxval**2+byval**2+bzval**2) 
       vpar_gk(i)=(vx(i)*bxval+vy(i)*byval+vz(i)*bzval)/bval ! dot product
       mu_gk(i)=(v_gk(i)**2-vpar_gk(i)**2)/(two*bval) !normalized magnetic moment
    enddo

    do i=1,nm
       call magnetic_coordinates_to_cylindrical_coordinates(zgc(i),xgc(i),rp_e(i),zp_e(i))

    enddo

    do i = 1, nm
       brval = br_mc_func(zgc(i), xgc(i))
       bzval = bz_mc_func(zgc(i), xgc(i))
       bphival = bphi_mc_func(zgc(i), xgc(i))
       angle = atan2(vy(i),vx(i))
       vr = vx(i)*cos(angle)+ vy(i)*sin(angle)
       vphi = -vx(i)*sin(angle)+vy(i)*cos(angle)
       call particle_to_guiding_center_location(ns, rp_e(i), phip_e(i), zp_e(i), vr,vphi,vz(i),brval,bphival,bzval, &
            & rg_e(i),phig_e(i),zg_e(i))

       xgc(i) = pfn_func(rg_e(i), zg_e(i)) !calculate the radial coordinate
       call interpolate_from_cylindrical_to_magnetic_coordinates(rg_e(i), zg_e(i), zgc(i), tor_shift_e)
       !ygc(i) = phig_e(i) + tor_shift_e ! a bug?
       ygc(i) = phig_e(i) - tor_shift_e
    enddo

    do i=1,nm
       if(xgc(i).gt.xupp .or. xgc(i) .lt. xlow)  touch_bdry_gc(i)=.true.
    end do

    !if(myid.eq.0) call calculate_possibility_density(vpar_gk,nm,100,minval(vpar_gk),maxval(vpar_gk))
    !if(myid.eq.0) call calculate_possibility_density(mu_gk,nm,100,minval(mu_gk),maxval(mu_gk))

    if(gk_spatial_loading_scheme==1 .and. gk_velocity_loading_scheme==1) then
       normalizing_factor=total_nm_gk(ns)/(twopi*toroidal_range*(xupp-xlow)*(vmax-vmin)**3)
       do i=1,nm
          call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,zgc(i),xgc(i),jacobian_val)
          ps_vol_gk(i)=abs(jacobian_val)/(normalizing_factor)
       enddo

    elseif(gk_spatial_loading_scheme.eq.1 .and. (gk_velocity_loading_scheme==2 .or. gk_velocity_loading_scheme==3)) then
       normalizing_factor=total_nm_gk(ns)/(twopi*toroidal_range*(xupp-xlow)*(sqrt(pi)*vt/vn_gk(ns))**3)
       do i=1,nm
          call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,zgc(i),xgc(i),jacobian_val)
          ps_vol_gk(i)=abs(jacobian_val)/(normalizing_factor*maxwellian0(v_gk(i)*vn_gk(ns),ns))
       enddo
    elseif(gk_spatial_loading_scheme.eq.2 .and. gk_velocity_loading_scheme.eq.1) then
       normalizing_factor=total_nm_gk(ns)/(vol*(vmax-vmin)**3)
       do i=1,nm
          ps_vol_gk(i)=one/(normalizing_factor)
       enddo
    elseif((gk_spatial_loading_scheme==2) .and. (gk_velocity_loading_scheme==2 .or. gk_velocity_loading_scheme==3)) then
       normalizing_factor = total_nm_gk(ns)/(vol*(sqrt(pi)*vt/vn_gk(ns))**3) 
       do i = 1, nm
          ps_vol_gk(i) = one/(normalizing_factor*maxwellian0(v_gk(i)*vn_gk(ns),ns))
       enddo
    endif

    !ptcl_num0_gk=ngk0*sqrt((mass_gk/(twopi*tgk0*kev))**3)*(ln**3*vn_gk**3)/normalizing_factor_e !explicitly cancel the exp(-v^2/vt^2) dependence, valid only for the case in which both physical equilibrium distribution and marker distribution are Maxwellian.
    do i = 1, nm
       ptcl_num0_gk(i) = ps_vol_gk(i) * f0(xgc(i), v_gk(i), ns)
    enddo
    call shift_toroidal(ygc(:), toroidal_range)

!!$  !assign the loaded particles to the corresponding processors, using the subroutines provided in pputil_yj.f90
    np_old = nm
    call init_pmove(zgc(:), np_old, twopi, ierr)
    call pmove(zgc(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(xgc(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ygc(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(vpar_gk(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(mu_gk(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(v_gk(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ps_vol_gk(:), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit

    call pmove(ptcl_num0_gk(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit

    call pmove2(touch_bdry_gc(:),np_old,np_new,ierr)
    if (ierr.ne.0) call ppexit
    nm = np_new

!!$   block
!!$  integer :: file_unit    
!!$  character(5):: filename
!!$  write(filename,'(a1,i4.4)') 'e',myid
!!$  open(newunit=file_unit, file=filename//'.txt')
!!$  do i=1,nm
!!$     write(file_unit,*)  xgc(i) ,zgc(i) ,rg_e(i),zg_e(i),i
!!$  enddo
!!$  close(file_unit)
!!$   end block

  end subroutine load_gk


  pure real(p_) function maxwellian(x, E, ns) result(f) !E in keV, f in SI unit
    use constants,only: p_, twopi, kev
    use gk_module, only : mass_gk
    use gk_profile_funcs, only : gkn_func, gkt_func
    real(p_), intent(in) :: x, E
    integer, intent(in) :: ns
    real(p_) :: n, T

    n = gkn_func(x,ns)
    T = gkt_func(x,ns)
    f = n*sqrt(mass_gk(ns)/(twopi*T*kev))**3*exp(-E/T)

  end function maxwellian

  pure real(p_) function f0(x,v,ns) result (z) ! v in unit of vn, f0 in unit 1/(Ln**3*vn**3)
    use constants,only: p_, two, twopi, kev
    use normalizing,only: Ln
    use gk_module,only: mass_gk, vn_gk
    use gk_profile_funcs,only : gkt_func, gkn_func
    implicit none
    integer, intent(in) :: ns
    real(p_),intent(in) :: x, v
    real(p_) :: v_si, te

    v_si = v*vn_gk(ns)
    te = gkt_func(x,ns)
    z = gkn_func(x,ns)*sqrt((mass_gk(ns)/(twopi*te*kev))**3)*exp(-mass_gk(ns)*v_si**2/(two*te*kev))
    z=z*(vn_gk(ns)**3*Ln**3)
  end function f0

  pure real(p_)  function maxwellian0(v, ns) result(z) !v in SI unit
    use constants, only: p_, two,kev
    use gk_module, only: mass_gk, tgk0 !as input
    implicit none
    real(p_), intent(in) :: v
    integer, intent(in) :: ns
    z = exp(-mass_gk(ns)*v*v/(two*tgk0(ns)*kev))
  end function maxwellian0

  pure real(p_)  function probability_sphere(v, ns) result(z) !v in SI unit
    use constants,only:p_, two,kev
    use gk_module,only: mass_gk,tgk0 !as input
    implicit none
    real(p_), intent(in) :: v
    integer, intent(in) :: ns
    z = v**2*exp(-mass_gk(ns)*v**2/(two*tgk0(ns)*kev))
  end function probability_sphere

  pure real(p_)  function probability(vx,vy,vz,ns) result(z) !v in SI unit
    use constants,only:p_, two,kev
    use gk_module,only: mass_gk,tgk0 !as input
    implicit none
    real(p_), intent(in) :: vx,vy,vz
    integer, intent(in) :: ns
    z=exp(-mass_gk(ns)*(vx**2+vy**2+vz**2)/(two*tgk0(ns)*kev))
  end function probability


  pure  subroutine particle_to_guiding_center_location(ns, r,phi,z,vr,vphi,vz,brval,bphival,bzval,rg,phig,zg)
    use constants,only: p_, one,two,half_pi
    use gk_module,only: mass_gk, charge_gk
    implicit none
    integer, intent(in) :: ns
    real(p_),intent(in) :: r, phi, z, vr, vphi, vz, brval, bphival, bzval
    real(p_),intent(out) :: rg, phig, zg
    real(p_):: v,bval

    v = sqrt(vr**2+vz**2+vphi**2)
    bval = sqrt(brval**2+bphival**2+bzval**2)

    rg=sqrt((r+mass_gk(ns)/(bval**2*charge_gk(ns))*(vphi*bzval-vz*bphival))**2+&
         &(mass_gk(ns)/(bval**2*charge_gk(ns))*(vz*brval-vr*bzval))**2)
    phig=phi+asin(mass_gk(ns)/(bval**2*charge_gk(ns))*(-vr*bzval+vz*brval)/rg)
    !zg=z+mass_gk(ns)/(bval**2*charge_gk(ns))*(-vr*bphival-vphi*brval) !wrong, pointed out by Yingfeng Xu
    zg=z+mass_gk(ns)/(bval**2*charge_gk(ns))*(vr*bphival-vphi*brval) !corrected.
  end subroutine particle_to_guiding_center_location

end module load_gk_mod
