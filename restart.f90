module restart
contains
  subroutine write_data_for_restarting(kend)
    !use fk_module
    use gk_module, only: nm_gk, w_gk, ps_vol_gk, lost_gc, v_gk, xgc, zgc, ygc, mu_gk, vpar_gk
    use domain_decomposition, only: myid
    use perturbation_field, only: potential, phix, phiy, phiz, apara, ax, ay, az 
    implicit none
    character(len=64) :: filename
    integer :: u, kend

    filename = 'myidxxxxx.pd'
    write(filename(5:9),'(i5.5)') myid
    open(newunit=u, file=filename, form='unformatted') 
    write(u) kend, nm_gk, w_gk, ps_vol_gk, lost_gc, v_gk, xgc, zgc, ygc, mu_gk, vpar_gk, &
         & potential, phix, phiy, phiz, apara, ax, ay, az !,&
    !& vt_i,vmin_i,vmax_i, normalizing_factor,&
    !& nmarker_i,active_i, touch_bdry_i, w_i,ps_vol_i, r_i,z_i,phi_i, vr_i,vz_i,vphi_i,&
    !& vr_i_mid, vz_i_mid, vphi_i_mid
    close(u)

  end subroutine write_data_for_restarting


  subroutine read_data_for_restarting(kstart)
    !use fk_module
    use gk_module, only: nm_gk, w_gk, ps_vol_gk, lost_gc, v_gk, xgc, zgc, ygc, mu_gk, vpar_gk
    use domain_decomposition, only: myid
    use perturbation_field, only: potential, phix, phiy, phiz, apara, ax, ay, az, apara_s, apara_h, ahx, ahy, ahz
    use mpi
    implicit none
    character(len=64) :: filename
    integer, intent(in) :: kstart
    integer :: u, kend, ierr

    filename = 'myidxxxxx.pd'
    write(filename(5:9),'(i5.5)') myid
    open(newunit=u, file=filename, form='unformatted', status='old')
      read(u) kend, nm_gk, w_gk, ps_vol_gk, lost_gc, v_gk, xgc, zgc, ygc, mu_gk, vpar_gk, &
           &  potential, phix, phiy, phiz, apara, ax, ay, az !,&
    !& vt_i, vmin_i, vmax_i, normalizing_factor, &
    !& nmarker_i, active_i, touch_bdry_i, w_i, ps_vol_i, r_i,z_i,phi_i, vr_i,vz_i,vphi_i, &
    !& vr_i_mid, vz_i_mid, vphi_i_mid
    close(u)

    apara_s = apara
    apara_h = 0.
    ahx = 0.
    ahy = 0.
    ahz = 0.
    
    if(kstart .ne. kend) then
       call MPI_FINALIZE(ierr)
       write(*,"(a,i5,a,i5,a,i5)") '***error****kend in the restarting file=', kend, ' kstart=', kstart, ' myid=', myid
       stop "***error***the restarting file is not for this time-step"
    endif

  end subroutine read_data_for_restarting
end module restart
