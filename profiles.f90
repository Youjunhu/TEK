module radial_profile_class
  use constants,only:p_

  type radial_profile
     real(p_),dimension(:),allocatable::  data 
   contains
     procedure :: set_profile, func
  end type radial_profile

contains
  pure real(p_) function func(this, pfn) result(z) !radial profile
    use constants,only:one
    use radial_module, only :npsi, pfn_npsi
    use interpolate_module, only: linear_1d_interpolate
    implicit none
    class(radial_profile), intent(in) :: this
    real(p_), intent(in) :: pfn !pfn is the normalized poloidal magnetic flux

    if(pfn>1) then
       z = 0.1*this%data(npsi)
    else
       call linear_1d_interpolate(npsi, pfn_npsi, this%data, pfn, z)
    endif
  end function func

  subroutine set_profile(this, filename, unit_of_data, my_unit, radial_coordinate_type)
    !set interpolating table by reading a file
    use constants,only: one, two, kev
    use radial_module,only: npsi, pfn_npsi, tfn_npsi !as input
    use magnetic_field, only : tfn_func_pfn
    use interpolate_module, only : linear_1d_interpolate_nonuniform
    implicit none
    class(radial_profile), intent(inout) :: this
    real(p_),intent(in) :: unit_of_data, my_unit
    character(*),intent(in):: filename,radial_coordinate_type
    integer,parameter:: max_num=3000
    real(p_):: radial_coordinate(max_num),tmp_ti(max_num)
    real(p_),dimension(:),allocatable::  tfn_sqrt_ndata
    real(p_)::tmp_y2b(npsi), tmp
    integer:: j, u, ndata
    real(p_),dimension(:),allocatable::  pfn_ndata, ti_ndata

    open(newunit=u,file=filename, status='old')
    do j=1,max_num
       read(u,*,end=111) radial_coordinate(j), tmp_ti(j) 
    enddo
111 close(u)
    ndata=j-1
    !write(*,*) 'number of data of the density radial profile=',ndata
    if(ndata.le.1) stop 'profile data are missing'

    allocate(ti_ndata(ndata))
    allocate(pfn_ndata(ndata))
    allocate(tfn_sqrt_ndata(ndata))
    allocate(this%data(npsi))

    ti_ndata(:) = tmp_ti(1:ndata)*unit_of_data/my_unit

    if (trim(radial_coordinate_type).eq.'toroidal-flux-sqrt') then
       tfn_sqrt_ndata(1:ndata)=radial_coordinate(1:ndata) 
!!$         call spline(sqrt(tfn_npsi),pfn_npsi,npsi,2.d30,2.d30,tmp_y2b) !prepare the second order derivative needed in the cubic spline interpolation
!!$         do j=1,ndata    
!!$            tfn_sqrt_ndata(j)=radial_coordinate(j)
!!$         enddo
!!$         do j=1,ndata !interpolating to get the corresponding pfn_ndata
!!$            call splint(sqrt(tfn_npsi),pfn_npsi,tmp_y2b,npsi,tfn_sqrt_ndata(j),pfn_ndata(j))
!!$         enddo
       do j=1,npsi
          tmp=tfn_func_pfn(pfn_npsi(j))
          call linear_1d_interpolate_nonuniform(ndata,tfn_sqrt_ndata,ti_ndata,sqrt(tmp),this%data(j))  
       enddo

    else if(trim(radial_coordinate_type).eq.'poloidal-flux') then   
       do j=1,npsi
          call linear_1d_interpolate_nonuniform(ndata,radial_coordinate,ti_ndata,pfn_npsi(j),this%data(j))  
       enddo

    else if(trim(radial_coordinate_type).eq.'poloidal-flux-sqrt') then   
       pfn_ndata(1:ndata)=radial_coordinate(1:ndata)**2
       do j=1,npsi
          call linear_1d_interpolate_nonuniform(ndata,pfn_ndata,ti_ndata,pfn_npsi(j),this%data(j))  
       enddo

    else 
       stop 'please specify the type of the radial grids used in the profile file'
    endif

  end subroutine set_profile
end module radial_profile_class


module adiabatic_e_profiles
  use radial_profile_class, only: radial_profile
  type(radial_profile) :: ne_object, te_object
contains

  subroutine initialize_adiabatic_electron()
    use constants, only: p_, kev
    implicit none
    character(100) :: ne_file, ne_radcor, te_file, te_radcor
    real(p_) :: ne_unit, te_unit
    namelist/adiabatic_electron_nmlt/ne_file, ne_unit, ne_radcor, te_file, te_unit, te_radcor
    integer :: u
    
    open(newunit=u,file='input.nmlt')
    read(u, adiabatic_electron_nmlt)
    close(u)

    call ne_object%set_profile(ne_file, ne_unit, 1.0d0, ne_radcor) 
    call te_object%set_profile(te_file, te_unit, kev, te_radcor) 

  end subroutine initialize_adiabatic_electron

end module adiabatic_e_profiles


module gk_radial_profiles
  use radial_profile_class, only: radial_profile
  type(radial_profile), allocatable :: density_object(:)
  type(radial_profile), allocatable :: temperature_object(:)
  type(radial_profile) :: alpha_ecrit, alpha_normc
  type(radial_profile) :: nalpha_object
contains
  subroutine initialize_gk_radial_profiles(nsm, density_file, density_unit, density_radcor, &
       &   temperature_file, temperature_unit, temperature_radcor)
    ! after this function is called, profile functions are ready to be used
    use constants, only : p_, kev, Mev, one, one
    implicit none
    integer, intent(in) :: nsm
    character(*), intent(in) :: density_file(nsm), temperature_file(nsm)
    character(*), intent(in) :: density_radcor(nsm), temperature_radcor(nsm)
    real(p_), intent(in) :: density_unit(nsm), temperature_unit(nsm)
    integer :: i, k, u

    allocate(density_object(nsm))
    allocate(temperature_object(nsm))

    do i = 1, nsm
       call density_object(i)%set_profile(density_file(i), density_unit(i), 1.0d0, density_radcor(i)) 
       call temperature_object(i)%set_profile(temperature_file(i), temperature_unit(i), kev, temperature_radcor(i)) 
    enddo !after this, the profile functions are ready to be used.
    
    nalpha_object = density_object(nsm)
    call alpha_ecrit%set_profile("cfetr/alpha_ecrit.txt", Mev, kev, "poloidal-flux-sqrt")
    call alpha_normc%set_profile("cfetr/alpha_normc.txt", one, one, "poloidal-flux-sqrt")

  end subroutine initialize_gk_radial_profiles
end module gk_radial_profiles


module gk_profile_funcs
  use constants, only: p_
  use gk_radial_profiles, only : temperature_object, density_object
  implicit none

contains
  
  pure real(p_) function gkn_func(x, ns) result (z) !unit: m^(-3)
    real(p_), intent(in) :: x
    integer, intent(in) :: ns

    z = density_object(ns)%func(x)

  end function gkn_func
  
  pure real(p_) function gkt_func(x, ns) result (z) !unit: kev
    real(p_), intent(in) :: x
    integer, intent(in) :: ns

    z = temperature_object(ns)%func(x)

  end function gkt_func

  pure real(p_) function gkdndx_func(x, ns) result (z)
    real(p_), intent(in) :: x
    integer, intent(in) :: ns
    real(p_), parameter :: dx = 1d-3
    
    z = (density_object(ns)%func(x+dx) -density_object(ns)%func(x-dx))/(2*dx)
    
  end function gkdndx_func
 
  
  pure real(p_) function gkdtdx_func(x, ns) result (z)
    real(p_), intent(in) :: x
    integer, intent(in) :: ns
    real(p_), parameter :: dx = 2d-3
    
    z = (temperature_object(ns)%func(x+dx) -temperature_object(ns)%func(x-dx))/(2*dx)
    
  end function gkdtdx_func


end module gk_profile_funcs
