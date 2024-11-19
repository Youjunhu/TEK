module restart_mod
contains
subroutine write_data_for_restarting(kend)
  use fk_module
  use gk_module,only: nm_gk,w_e, ptcl_num0_e,touch_bdry_e,radcor_e,theta_e,alpha_e,mu_e,vpar_e
  use domain_decomposition,only:myid
!  use perturbation_field_matrix,only: mf_par_left,mf_x_left,mf_y_left !field at the left-boundary (smaller theta) of the present cell 
  implicit none
  character(len=64)::cfile
  integer:: ufile !file nunit
  integer:: kend

  cfile = 'myidxxxxx.pd'
  write(cfile(5:9),'(i5.5)') myid
  open(newunit=ufile,file=cfile,form='unformatted') !newunit is a keyword argument (as output) of open()
  write(ufile) kend,  vt_i,vmin_i,vmax_i, normalizing_factor,&
       & nmarker_i,active_i, touch_bdry_i, w_i,ps_vol_i, r_i,z_i,phi_i, vr_i,vz_i,vphi_i,&
       & vr_i_mid, vz_i_mid, vphi_i_mid,&
       & nm_gk,w_e,ptcl_num0_e, touch_bdry_e,radcor_e,theta_e,alpha_e,mu_e,vpar_e
  close(ufile)

end subroutine write_data_for_restarting


subroutine read_data_for_restarting(kstart)
  use fk_module
  use gk_module,only: nm_gk,w_e,ptcl_num0_e,touch_bdry_e,radcor_e,theta_e,alpha_e,mu_e,vpar_e
  use domain_decomposition,only:myid

  use communication_connection
  use mpi
  implicit none
  character(len=64)::cfile
  integer:: ufile !file nunit
  integer:: kstart, kend, ierr
  cfile = 'myidxxxxx.pd'
  write(cfile(5:9),'(i5.5)') myid
  open(newunit=ufile,file=cfile,form='unformatted') !newunit is a keyword argument (as output) of open()

  read(ufile) kend, vt_i,vmin_i,vmax_i,normalizing_factor,&
       & nmarker_i,active_i, touch_bdry_i, w_i,ps_vol_i, r_i,z_i,phi_i, vr_i,vz_i,vphi_i, &
       & vr_i_mid, vz_i_mid, vphi_i_mid,&
       & nm_gk,w_e,ptcl_num0_e, touch_bdry_e,radcor_e,theta_e,alpha_e,mu_e,vpar_e
  close(ufile)

  if(kstart.ne.kend+1) then
     call MPI_FINALIZE(ierr)
     write(*,"(a,i5,a,i5,a,i5)") '***error****kend in the restarting file=',kend, ' kstart=',kstart, ' myid=',myid
     stop "***error***the restarting file is not for this time-step"
!call MPI_Abort(MPI_COMM_WORLD, errcode, ierr) !In a catastrophic error condition, I may use this, other cases I prefer to use mpi_finalize and then use fortran "stop"
  endif
!  call update_efield_at_right_boundary_of_present_cell() !Electric field value at right-boundary of the present cell is needed when pushing particle weight, smoothing and taking z derivatives
!  call update_bfield_at_right_boundary_of_present_cell() !magnetic field value at right-boundary of the present cell is needed when pushing particle weight, smoothing and taking z derivatives
!  call update_efield_at_second_left_boundary() !electric field value at the second left-boundary of the present cell is needed when smoothing and taking z derivatives
!  call update_bfield_at_second_left_boundary() !magnetic field value at the second left-boundary of the present cell is needed when smoothing and taking z derivatives

end subroutine read_data_for_restarting
end module restart_mod
