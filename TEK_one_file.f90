module constants
  implicit none
  save
  integer,parameter:: p_=kind(1.0d0) !precision
  real(p_),parameter :: coulomb_log=15._p_ !assumed to be a constant
  real(p_),parameter :: kev=1.6022d-16   !unit J
  real(p_),parameter :: Mev=1.6022d-13   !unit J
  real(p_),parameter :: elementary_charge=1.6022d-19   !in unit of Coulomb
  real(p_),parameter :: epsilon0=8.8542d-12 !Permittivity of free space 
  real(p_),parameter :: pi=3.1415926_p_
  real(p_),parameter :: twopi=pi*2.0_p_
  real(p_),parameter :: fourpi=pi*4.0_p_
  real(p_),parameter :: half_pi=pi*0.5_p_
  real(p_),parameter :: mu0=fourpi*1.0d-7 !permeability in SI unit
  real(p_),parameter :: zero=0.0_p_
  real(p_),parameter :: one=1.0_p_
  real(p_),parameter :: two=2.0_p_
  real(p_),parameter :: three=3.0_p_
  real(p_),parameter :: four=4.0_p_
  real(p_),parameter :: five=5.0_p_
  real(p_),parameter :: six=6.0_p_
  real(p_),parameter :: seven=7.0_p_
  real(p_),parameter :: eight=8.0_p_
  real(p_),parameter :: nine=9.0_p_
  real(p_),parameter :: ten=10.0_p_
  real(p_),parameter :: eleven=11.0_p_
  real(p_),parameter :: twelve=12.0_p_
  real(p_),parameter :: thirteen=13.0_p_
  real(p_),parameter :: fourteen=14.0_p_
  real(p_),parameter :: fifteen=15.0_p_
  real(p_),parameter :: one_half=0.5_p_
  real(p_),parameter :: one_third=one/three
  real(p_),parameter :: one_fifth=0.2_p_
  real(p_),parameter :: three_halfs=1.5_p_
  real(p_),parameter :: c=299792458d0 !speed of light in vaccuum
  real(p_),parameter:: atom_mass_unit=1.660539066d-27  !kg
  complex(p_), parameter :: ii=(0.0_p_, 1._p_) 
end module constants

module normalizing
  use constants, only : p_
  implicit none
  save
  real(p_), parameter :: Ln=1.0d0 !length unit (in unit of meter), do not change the value because I have assumed this value (1.0) in some parts of the code
  real(p_), parameter :: bn=1.0d0 !nmagnetic field strength (in unit of Tesla), do not change the value
  real(p_) :: tu  !temperature unit in keV
  real(p_) :: nu !SI unit (m^-3)
  real(p_) :: qu   !elementary_charge  in SI unit (Coulomb)
  real(p_) :: vu !in SI unit
  real(p_) :: dtao_main
end module normalizing

module  poloidal_flux_2d !poloidal flux and its partial derivatives in (R,Z) plane
  use constants,only:p_
  implicit none
  save
  integer:: nx, nz !number of gridpoints in R and Z directions (specified in G-file)
  real(p_), dimension(:), allocatable :: xarray, zarray ! R and Z array
  real(p_), dimension(:,:), allocatable :: psi !psi is related to poloidal magnetic field by Bp=\nabla{psi}\times\nabla{phi}
  real(p_),dimension(:,:),allocatable :: psi_x, psi_z  !partial derivatives
  real(p_),dimension(:,:),allocatable :: psi_xx, psi_zz, psi_xz, psi_zx !partial derivatives
  real(p_),dimension(:,:),allocatable :: psi_gradient, y2a_gradient
  !strength of the gradient of the poloidal flux, y2a_gradient is the second derivative array used in 2D spline interpolation
  real(p_), dimension(:,:), allocatable :: y2a_psi, y2a_psi_x,y2a_psi_z,y2a_psi_xx,y2a_psi_zz,y2a_psi_xz,y2a_psi_zx
  ! y2a_* is the second derivative array used in 2D spline interpolation
end module poloidal_flux_2d


module radial_module
  use constants, only: p_
  implicit none
  save
  integer :: npsi
  real(p_) :: r_axis, z_axis, baxis, psi_axis, psi_lcfs, minor_a, sign_bphi
  integer :: j_fixed
  real(p_) :: radcor_fixed !the radcor of the center of computational region
  real(p_), dimension(:), allocatable :: fpsi, ffprime, fprime, qpsi
  real(p_), dimension(:), allocatable :: pressure, pprime
  real(p_), dimension(:), allocatable :: psi_1d, pfn_npsi, tfn_npsi
  real(p_), dimension(:), allocatable :: q_with_sign
  real(p_), dimension(:), allocatable :: qrad
end module radial_module

module boundary
  use constants, only: p_
  save
  integer :: nlim, np_lcfs
  real(p_), dimension(:), allocatable :: rlim, zlim, x_lcfs, z_lcfs
end module boundary


module mapping_module !from cylindrical coordinates to magnetic coordinates
  use constants, only: p_
  implicit none
  save
  integer, parameter :: nx_mapping=100, nz_mapping=100
  real(p_):: r_cyl(nx_mapping), z_cyl(nz_mapping)
  real(p_):: radcor(nx_mapping,nz_mapping)
  real(p_):: theta_a(nx_mapping,nz_mapping), theta_b(nx_mapping,nz_mapping)
  real(p_):: tor_shift_a(nx_mapping,nz_mapping), tor_shift_b(nx_mapping,nz_mapping)
  real(p_):: dr,dz
  integer:: i0, j0 !index of the point at magnetic axis
  real(p_):: dtheta_dr(nx_mapping,nz_mapping), dtheta_dz(nx_mapping,nz_mapping)
  real(p_):: ddelta_dr_a(nx_mapping,nz_mapping), ddelta_dz_a(nx_mapping,nz_mapping)
  real(p_):: ddelta_dr_b(nx_mapping,nz_mapping), ddelta_dz_b(nx_mapping,nz_mapping)
  real(p_):: dradial_dr(nx_mapping,nz_mapping), dradial_dz(nx_mapping,nz_mapping)
end module mapping_module


module domain_decomposition
  use constants, only: p_
  implicit none
  integer :: numprocs, myid, nvp
  integer :: GRID_COMM, TUBE_COMM
  integer :: GCLR, TCLR, ntube, my_left, my_right, my_left2, my_right2
  integer :: multi_eq_cells
  !1D domain decomposion along theta coordinates, a processor treats [theta_start, theta_start+dtheta2]
  real(p_) :: dtheta2, theta_start 
  integer :: ipol_eq !index of theta_start in the equilibrium zgrid
  real(p_), allocatable :: dvol(:)
end module domain_decomposition

module control_parameters
  use constants,only:p_, elementary_charge,kev
  implicit none
  save
  integer  :: kstart,kend
  real(p_) :: dt_omega_i_axis
  character(100):: poloidal_angle_type
  logical :: adiabatic_electrons
  logical :: store_restart_data
  integer  :: iplot_mode_structure
  logical  :: filter_radial, diagnosis !if ture, output additional testing information
  integer  :: order_faraday_scheme
  integer :: space_charge_switch, fk_switch
  integer :: ismooth
  integer :: nh_min !non-negative integer
  integer :: nh_max !non-negative integer, nh_max musch be larger than or equal to nh_min
  !toroidal harmonics included are exp(i*nsegment*y) with i = nh_min, ..., nh_max and i = -nh_max, ..., -nh_min
end module control_parameters


module perturbation_field
  use constants, only: p_
  implicit none
  save
  real(p_), dimension(:,:,:), allocatable :: potential, phix, phiy, phiz !electrostatic potential and its derivatives
  real(p_), dimension(:,:,:), allocatable :: apara, ax, ay, az !parallel component of the vector potnetial and its derivative
  real(p_), dimension(:,:,:), allocatable :: apara_h, apara_s, apara_s_old !used in the mixed-variable pullback method
  real(p_), dimension(:,:,:), allocatable :: ahx, ahy, ahz

  real(p_), dimension(:,:), allocatable :: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
  real(p_), dimension(:,:), allocatable :: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right

contains

  subroutine allocate_field_matrix(m, n)
    implicit none
    integer, intent(in) :: m, n !toroidal, radial

    allocate(potential(m, n, 2)) ! (:,:,1) is at the left-boundary of the present z-cell, (:,:,2) is at the right-boundary
    allocate(phix(m, n, 2)) !dphi/dx
    allocate(phiy(m, n, 2)) !dphi/dy
    allocate(phiz(m, n, 2)) !dphi/dz

    allocate(apara  (m,n,2))
    allocate(apara_h(m,n,2))
    allocate(apara_s(m,n,2))
    allocate(apara_s_old(m,n,2))
    allocate(ax(m,n,2)) !dApar/dx
    allocate(ay(m,n,2)) !dApar/dy
    allocate(az(m,n,2)) !dApar/dz
    allocate(ahx(m,n,2)) !dApar_h/dx
    allocate(ahy(m,n,2)) !dApar_h/dy
    allocate(ahz(m,n,2)) !dApar_h/dz

    allocate(ef_cyl_r_left(m+1,n))
    allocate(ef_cyl_z_left(m+1,n))
    allocate(ef_cyl_phi_left(m+1,n))
    allocate(ef_cyl_r_right(m+1,n))
    allocate(ef_cyl_z_right(m+1,n))
    allocate(ef_cyl_phi_right(m+1,n))

  end subroutine allocate_field_matrix

end module perturbation_field

MODULE pputil
  use constants, only: p_
  use domain_decomposition, only: myid, nvp, GCLR, TCLR, TUBE_COMM
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ppexit, init_pmove, pmove, pmove2
  INTEGER, SAVE :: pmove_tag=0
  REAL(P_), DIMENSION(:), ALLOCATABLE, SAVE :: s_buf, r_buf
  logical, DIMENSION(:), ALLOCATABLE, SAVE :: s_buf2, r_buf2 !yj added
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: s_counts, s_displ
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: r_counts, r_displ
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ipsend, iphole

CONTAINS
  SUBROUTINE init_pmove(xp, np, lz, ierr)
!INCLUDE 'mpif.h'
    use mpi
    REAL(P_), DIMENSION(:), INTENT(in) :: xp
    INTEGER, INTENT(in) :: np
    REAL(P_), INTENT(in) :: lz
    INTEGER, INTENT(out) :: ierr
    INTEGER :: nsize, ksize
    INTEGER :: i, ip, iz, ih, iwork
    REAL(P_) :: dzz
    INTEGER, DIMENSION(0:nvp-1) :: isb

    IF( .not. ALLOCATED(s_counts) ) ALLOCATE(s_counts(0:nvp-1))
    IF( .not. ALLOCATED(s_displ) ) ALLOCATE(s_displ(0:nvp-1))
    IF( .not. ALLOCATED(r_counts) ) ALLOCATE(r_counts(0:nvp-1))
    IF( .not. ALLOCATED(r_displ) ) ALLOCATE(r_displ(0:nvp-1))

    !----------------------------------------------------------------------
    !              1.  Construct send buffer

    dzz = lz / nvp
    s_counts = 0
    DO ip = 1, np
       iz = INT((xp(ip)+ lz/2)/dzz)
       IF( iz .ne. GCLR )  s_counts(iz) = s_counts(iz) + 1
    END DO
    s_displ(0) = 0
    DO i=1,nvp-1
       s_displ(i) = s_displ(i-1) + s_counts(i-1)
    END DO

    nsize = sum(s_counts)
    IF( .not. ALLOCATED(s_buf) ) THEN
       ksize = 2*nsize         ! To prevent too much futur reallocations
       ALLOCATE(s_buf(1:ksize))
       ALLOCATE(s_buf2(1:ksize))
       ALLOCATE(ipsend(1:ksize))
       ALLOCATE(iphole(1:ksize))
    ELSE IF ( SIZE(s_buf) .LT. nsize ) THEN
       DEALLOCATE(s_buf)
       DEALLOCATE(s_buf2)
       DEALLOCATE(ipsend)
       DEALLOCATE(iphole)
       ALLOCATE(s_buf(1:nsize))
       ALLOCATE(s_buf2(1:nsize))
       ALLOCATE(ipsend(1:nsize))
       ALLOCATE(iphole(1:nsize))
    END IF


    ! Construct (sorted) pointers to holes
    isb(0:nvp-1) = s_displ(0:nvp-1)
    ih = 0
    DO ip = 1, np
       iz = INT((xp(ip)+ lz/2)/dzz)
       IF( iz .ne. GCLR ) THEN
          isb(iz) = isb(iz) + 1
          ipsend(isb(iz)) = ip
          ih = ih+1
          iphole(ih) = ip
       END IF
    END DO

    ! Construct receive buffer
    CALL MPI_ALLTOALL(s_counts, 1, MPI_INTEGER, &
         & r_counts, 1, MPI_INTEGER, TUBE_COMM, ierr)
    r_displ(0) = 0
    DO i=1,nvp-1
       r_displ(i) = r_displ(i-1) + r_counts(i-1)
    END DO

    nsize = sum(r_counts)
    IF( .not. ALLOCATED(r_buf) ) THEN
       ksize=2*nsize         ! To prevent too much futur reallocations
       ALLOCATE(r_buf(1:ksize))
       ALLOCATE(r_buf2(1:ksize))
    ELSE IF ( SIZE(r_buf) .LT. nsize ) THEN
       DEALLOCATE(r_buf)
       DEALLOCATE(r_buf2)
       ALLOCATE(r_buf(1:nsize))
       ALLOCATE(r_buf2(1:nsize))
    END IF

    ! Check for particle array overflow
    ierr = 0
    nsize = np - sum(s_counts) + sum(r_counts)
    if( nsize .gt. size(xp) ) then
       write(*,*) 'myid=', myid, 'Particle array overflow', nsize, '>', size(xp)
       ierr = 1
    end if

    pmove_tag = 101

  END SUBROUTINE init_pmove

  !===========================================================================
  SUBROUTINE pmove(xp, np_old, np_new, ierr)
    !
    !INCLUDE 'mpif.h'
    use mpi
    !
    REAL(P_), DIMENSION(:), INTENT(inout) :: xp
    INTEGER, INTENT(in) :: np_old
    INTEGER, INTENT(out) :: np_new
    INTEGER, INTENT(out) :: ierr
    !
    !  Local vars
    INTEGER :: nsize
    INTEGER :: i, ip, iz, ih, isrc
    INTEGER :: nhole, mhole, nrrecv, nrsend, nptot_old, nptot
    INTEGER :: ind, count, tot_count, iwork
    !
    !  Local arrays
    INTEGER, DIMENSION(1:nvp) :: s_requ, r_requ, id_source
    INTEGER :: isrt, iend
    INTEGER :: status(MPI_STATUS_SIZE), arr_status(MPI_STATUS_SIZE,nvp)
    !
    !----------------------------------------------------------------------
    !              1.  Fill send buffer
    !
    DO i=0,nvp-1
       IF( s_counts(i) .GT. 0 ) THEN
          isrt = s_displ(i)+1
          iend = s_displ(i)+s_counts(i)
          s_buf(isrt:iend) = xp(ipsend(isrt:iend))
       END IF
    END DO
    !----------------------------------------------------------------------
    !              2.   Initiate non-blocking send/receive
    !
    pmove_tag = pmove_tag+1
    nrrecv=0             !......... Start non-blocking receive
    DO i=0,nvp-1
       IF(r_counts(i) .GT. 0 ) THEN
          nrrecv=nrrecv+1
          id_source(nrrecv) = i
          isrt = r_displ(i)+1
          CALL MPI_IRECV(r_buf(isrt), r_counts(i), MPI_REAL8,&
               & i, pmove_tag, TUBE_COMM, r_requ(nrrecv), ierr)
       END IF
    END DO
    nrsend=0             !......... Start non-blocking SYNCHRONOUS send
    DO i=0,nvp-1
       IF(s_counts(i) .GT. 0 ) THEN
          nrsend=nrsend+1
          isrt = s_displ(i)+1
          CALL MPI_ISSEND(s_buf(isrt), s_counts(i), MPI_REAL8,&
               & i, pmove_tag, TUBE_COMM, s_requ(nrsend), ierr)
       END IF
    END DO
    !
    !-------------------------------------------------------------
    !3.   Remove holes and compress particle arrays

    nhole = sum(s_counts)
    ip = np_old
    DO ih = nhole, 1, -1
       xp(iphole(ih)) = xp(ip)
       ip = ip-1
    END DO
    np_new = ip
    !
    !-------------------------------------------------------------
    !4.   Store incoming particle to the partitcle arrays

    tot_count = 0
    DO i=1,nrrecv
       CALL MPI_WAITANY(nrrecv, r_requ, ind, status, ierr)
       isrc = id_source(ind)
       CALL MPI_GET_COUNT(status, MPI_REAL8, count, ierr)
       IF( count .ne. r_counts(isrc) ) THEN
          WRITE(*,*) 'PE',myid, '  Counts mismatched from PE',isrc,&
               & count,  r_counts(isrc)
       END IF
       tot_count =  tot_count+count
    END DO
    !
    IF( tot_count .GT. 0 ) THEN
       isrt = np_new + 1
       iend = np_new + tot_count
       xp(isrt:iend) = r_buf(1:tot_count)
       np_new = iend
    END IF

    !---------------------------------------------------
    ! 5.   Epilogue
    !... Wait for any non-blocking comm. requests
    IF( nrsend.gt.0) CALL MPI_WAITALL(nrsend, s_requ, arr_status, ierr)
    !
    !... Check consistency
    ierr = 0
    CALL MPI_ALLREDUCE(np_old, nptot_old, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(np_new, nptot, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    IF( nptot .ne. nptot_old ) THEN
       IF(myid.eq.0) WRITE(*,*) 'PMOVE: mismatch in total numbers:', &
            & nptot_old, nptot
       ierr = 1
    END IF

  END SUBROUTINE pmove


  SUBROUTINE pmove2(xp, np_old, np_new, ierr) !differ from pmove() in that it deal with logical array, instead of real array, by yjhu
    !
    !INCLUDE 'mpif.h'
    use mpi
    !
    logical, DIMENSION(:), INTENT(inout) :: xp
    INTEGER, INTENT(in) :: np_old
    INTEGER, INTENT(out) :: np_new
    INTEGER, INTENT(out) :: ierr
    !
    !  Local vars
    INTEGER :: nsize
    INTEGER :: i, ip, iz, ih, isrc
    INTEGER :: nhole, mhole, nrrecv, nrsend, nptot_old, nptot
    INTEGER :: ind, count, tot_count, iwork
    !
    !  Local arrays
    INTEGER, DIMENSION(1:nvp) :: s_requ, r_requ, id_source
    INTEGER :: isrt, iend
    INTEGER :: status(MPI_STATUS_SIZE), arr_status(MPI_STATUS_SIZE,nvp)
    !
    !----------------------------------------------------------------------
    !              1.  Fill send buffer
    !
    DO i=0,nvp-1
       IF( s_counts(i) .GT. 0 ) THEN
          isrt = s_displ(i)+1
          iend = s_displ(i)+s_counts(i)
          s_buf2(isrt:iend) = xp(ipsend(isrt:iend))
       END IF
    END DO
    !----------------------------------------------------------------------
    !              2.   Initiate non-blocking send/receive
    !
    pmove_tag = pmove_tag+1
    nrrecv=0             !......... Start non-blocking receive
    DO i=0,nvp-1
       IF(r_counts(i) .GT. 0 ) THEN
          nrrecv=nrrecv+1
          id_source(nrrecv) = i
          isrt = r_displ(i)+1
          CALL MPI_IRECV(r_buf2(isrt), r_counts(i), MPI_logical,&
               & i, pmove_tag, TUBE_COMM, r_requ(nrrecv), ierr)
       END IF
    END DO
    nrsend=0             !......... Start non-blocking SYNCHRONOUS send
    DO i=0,nvp-1
       IF(s_counts(i) .GT. 0 ) THEN
          nrsend=nrsend+1
          isrt = s_displ(i)+1
          CALL MPI_ISSEND(s_buf2(isrt), s_counts(i), MPI_logical,&
               & i, pmove_tag, TUBE_COMM, s_requ(nrsend), ierr)
       END IF
    END DO
    !
    !----------------------------------------------------------------------
    !              3.   Remove holes and compress part. arrays
    !
    nhole = sum(s_counts)
    ip = np_old
    DO ih = nhole, 1, -1
       xp(iphole(ih)) = xp(ip)
       ip = ip-1
    END DO
    np_new = ip
    !
    !----------------------------------------------------------------------
    !              4.   Store incoming part. to the part. arrays
    !
    tot_count = 0
    DO i=1,nrrecv
       CALL MPI_WAITANY(nrrecv, r_requ, ind, status, ierr)
       isrc = id_source(ind)
       CALL MPI_GET_COUNT(status, MPI_logical, count, ierr)
       IF( count .ne. r_counts(isrc) ) THEN
          WRITE(*,*) 'PE',myid, '  Counts mismatched from PE',isrc,&
               & count,  r_counts(isrc)
       END IF
       tot_count =  tot_count+count
    END DO
    !
    IF( tot_count .GT. 0 ) THEN
       isrt = np_new + 1
       iend = np_new + tot_count
       xp(isrt:iend) = r_buf2(1:tot_count)
       np_new = iend
    END IF
    !
    !----------------------------------------------------------------------
    !              5.   Epilogue
    !
    !... Wait for any non-blocking comm. requests
    IF( nrsend.gt.0) CALL MPI_WAITALL(nrsend, s_requ, arr_status, ierr)
    !
    !... Check consistency
    ierr = 0
    CALL MPI_ALLREDUCE(np_old, nptot_old, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(np_new, nptot, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    IF( nptot.ne.nptot_old ) THEN
       IF(myid.eq.0) WRITE(*,*) 'PMOVE: mismatch in total numbers:',&
            & nptot_old, nptot
       ierr = 1
    END IF
    !    call ppsum(ierr)
    !
    !----------------------------------------------------------------------!
  END SUBROUTINE pmove2


  SUBROUTINE ppexit
    INTEGER :: ierr
    CALL MPI_FINALIZE(ierr)
    STOP
  END SUBROUTINE ppexit

END MODULE pputil
module lapack
implicit none

! This is the precision that LAPACK "d" routines were compiled with (typically
! double precision, unless a special compiler option was used while compiling
! LAPACK). This "dp" is only used in lapack.f90
! The "d" routines data type is defined as "double precision", so
! we make "dp" the same kind as 0.d0 ("double precision"), so
! as long as LAPACK and this file were compiled with the same compiler options,
! it will be consistent. (If for example all double precision is promoted to
! quadruple precision, it will be promoted both in LAPACK and here.)
integer, parameter :: dp=kind(0.d0)

interface

    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                       EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, &
                       IWORK, INFO )
    import :: dp
    CHARACTER          EQUED, FACT, TRANS
    INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
    REAL(dp)           RCOND
    INTEGER            IPIV( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), &
                       C( * ), FERR( * ), R( * ), WORK( * ), X( LDX, * )
    END SUBROUTINE

    SUBROUTINE ZGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE ZGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, &
                       EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, &
                       WORK, RWORK, INFO )
    import :: dp
    CHARACTER          EQUED, FACT, TRANS
    INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
    REAL(dp)           RCOND
    INTEGER            IPIV( * )
    REAL(dp)           BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )
    COMPLEX(dp)        A( LDA, * ), AF( LDAF, * ), B( LDB, * ), WORK( * ), &
                       X( LDX, * )
    END SUBROUTINE

    SUBROUTINE DGBSV( N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO )
    import :: dp
    INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           AB( LDAB, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DSYSV( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYTRS( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B, &
                       LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK, &
                       IWORK, INFO )
    import :: dp
    CHARACTER          FACT, UPLO
    INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
    REAL(dp)           RCOND
    INTEGER            IPIV( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), AF( LDAF, * ), B( LDB, * ), &
                       BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
    END SUBROUTINE

    SUBROUTINE DSYEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, &
                       LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), W( * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, &
                       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                       LWORK, IWORK, IFAIL, INFO )
    import :: dp
    CHARACTER          JOBZ, RANGE, UPLO
    INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
    REAL(dp)           ABSTOL, VL, VU
    INTEGER            IFAIL( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * ), &
                       Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, &
                       VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, &
                       LWORK, IWORK, IFAIL, INFO )
    import :: dp
    CHARACTER          JOBZ, RANGE, UPLO
    INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
    REAL(dp)           ABSTOL, VL, VU
    INTEGER            IFAIL( * ), IWORK( * )
    REAL(dp)           A( LDA, * ), W( * ), WORK( * ), &
                       Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, &
                      BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          JOBVL, JOBVR
    INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
    REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), &
                       B( LDB, * ), BETA( * ), VL( LDVL, * ), &
                       VR( LDVR, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB, &
                       ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, ILO, IHI, &
                       LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV, WORK, &
                       LWORK, IWORK, BWORK, INFO )
    import :: dp
    CHARACTER          BALANC, JOBVL, JOBVR, SENSE
    INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
    REAL(dp)           ABNRM, BBNRM
    LOGICAL            BWORK( * )
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), ALPHAI( * ), ALPHAR( * ), B( LDB, * ), &
                       BETA( * ), LSCALE( * ), RCONDE( * ), RCONDV( * ), &
                       RSCALE( * ), VL( LDVL, * ), VR( LDVR, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                      LDVR, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          JOBVL, JOBVR
    INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), &
                       WORK( * ), WR( * )
    END SUBROUTINE

    SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI, &
                       VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, &
                       RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
    import :: dp
    CHARACTER          BALANC, JOBVL, JOBVR, SENSE
    INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           ABNRM
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), RCONDE( * ), RCONDV( * ), &
                       SCALE( * ), VL( LDVL, * ), VR( LDVR, * ), &
                       WI( * ), WORK( * ), WR( * )
    END SUBROUTINE

    SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, &
                      WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          JOBVL, JOBVR
    INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           RWORK( * )
    COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
                       WORK( * )
    END SUBROUTINE

    SUBROUTINE ZGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, &
                       LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, &
                       RCONDV, WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          BALANC, JOBVL, JOBVR, SENSE
    INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
    REAL(dp)           ABNRM
    REAL(dp)           RCONDE( * ), RCONDV( * ), RWORK( * ), SCALE( * )
    COMPLEX(dp)        A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), W( * ), &
                       WORK( * )
    END SUBROUTINE

    SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, &
                       LWORK, IWORK, LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
    END SUBROUTINE

    REAL(dp) FUNCTION DLAMCH( CMACH )
    import :: dp
    CHARACTER          CMACH
    END FUNCTION

    INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
    CHARACTER*( * )    NAME, OPTS
    INTEGER            ISPEC, N1, N2, N3, N4
    END FUNCTION

    SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
    import :: dp
    INTEGER            INFO, LDA, M, N
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * )
    END SUBROUTINE

    SUBROUTINE ZGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
    import :: dp
    CHARACTER          TRANS
    INTEGER            INFO, LDA, LDB, N, NRHS
    INTEGER            IPIV( * )
    COMPLEX(dp)         A( LDA, * ), B( LDB, * )
    END SUBROUTINE

    SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    COMPLEX(dp)        A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
    import :: dp
    INTEGER            INFO, LDA, M, N
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * )
    END SUBROUTINE

    SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LWORK, N
    INTEGER            IPIV( * )
    REAL(dp)           A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LWORK, N
    REAL(dp)           RWORK( * ), W( * )
    COMPLEX(dp)        A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZHEEVD( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, &
                       LRWORK, IWORK, LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, LDA, LIWORK, LRWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           RWORK( * ), W( * )
    COMPLEX(dp)        A( LDA, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZHEGVD( ITYPE,  JOBZ,  UPLO,  N,  A,  LDA,  B, LDB, W, &
                       WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, &
                       INFO )
    import :: dp
    CHARACTER          JOBZ, UPLO
    INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           RWORK( * ), W( * )
    COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
                       WORK, LWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
    REAL(dp)           RCOND
    INTEGER            JPVT( * )
    REAL(dp)           A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE ZGELSY( M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, &
                       WORK, LWORK, RWORK, INFO )
    import :: dp
    INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
    REAL(dp)           RCOND
    INTEGER            JPVT( * )
    REAL(dp)           RWORK( * )
    COMPLEX(dp)        A( LDA, * ), B( LDB, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, &
                       LDVT, WORK, LWORK, INFO )
    import :: dp
    CHARACTER          JOBU, JOBVT
    INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
    REAL(dp)           A( LDA, * ), S( * ),  U( LDU,  * ), VT( LDVT, * ), &
                       WORK( * )
    END SUBROUTINE

    SUBROUTINE ZGESVD( JOBU, JOBVT,  M,  N,  A,  LDA, S, U, LDU, VT, LDVT, &
                       WORK, LWORK, RWORK, INFO )
    import :: dp
    CHARACTER          JOBU, JOBVT
    INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
    REAL(dp)           RWORK( * ), S( * )
    COMPLEX(dp)        A( LDA, * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
    END SUBROUTINE

    SUBROUTINE DSTEVD( JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, &
                       LIWORK, INFO )
    import :: dp
    CHARACTER          JOBZ
    INTEGER            INFO, LDZ, LIWORK, LWORK, N
    INTEGER            IWORK( * )
    REAL(dp)           D( * ), E( * ), WORK( * ), Z( LDZ, * )
    END SUBROUTINE

    SUBROUTINE XERBLA( SRNAME, INFO )
    CHARACTER*(*)      SRNAME
    INTEGER            INFO
    END SUBROUTINE

! BLAS

    SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
    import :: dp
    INTEGER INCX,INCY,N
    COMPLEX(dp) ZX(*),ZY(*)
    END SUBROUTINE

    SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
    import :: dp
    integer :: INCX, INCY, N
    real(dp) :: DA, DX(*), DY(*)
    END SUBROUTINE

    SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    import :: dp
    DOUBLE PRECISION ALPHA,BETA
    INTEGER K,LDA,LDB,LDC,M,N
    CHARACTER TRANSA,TRANSB
    REAL(dp) A(LDA,*),B(LDB,*),C(LDC,*)
    END SUBROUTINE

    real(dp) FUNCTION DNRM2(N,X,INCX)
    import :: dp
    integer :: INCX, N
    real(dp) :: X(*)
    END FUNCTION

    SUBROUTINE DSCAL(N,DA,DX,INCX)
    import :: dp
    real(dp) :: DA, DX(*)
    integer :: INCX, N
    END SUBROUTINE

    SUBROUTINE DSYMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    import :: dp
    REAL(dp) ALPHA,BETA
    INTEGER LDA,LDB,LDC,M,N
    CHARACTER SIDE,UPLO
    REAL(dp) A(LDA,*),B(LDB,*),C(LDC,*)
    END SUBROUTINE

    SUBROUTINE DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
    import :: dp
    INTEGER  INFO, LDA, LWORK, M, N
    REAL(dp) A(LDA, *), TAU(*), WORK(*)
    END SUBROUTINE

    SUBROUTINE DORGQR(M, N, K, A, LDA, TAU, WORK, LWORK, INFO)
    import :: dp
    INTEGER  INFO, K, LDA, LWORK, M, N
    REAL(dp) A(LDA,*), TAU(*), WORK(*)
    END SUBROUTINE

    SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO )
    import :: dp
    CHARACTER          UPLO
    INTEGER            INFO, LDA, N
    REAL(dp)           A( LDA, * )
    END SUBROUTINE

    SUBROUTINE DTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, LDA, B, LDB, INFO )
    import :: dp
    CHARACTER          DIAG, TRANS, UPLO
    INTEGER            INFO, LDA, LDB, N, NRHS
    REAL(dp)           A( LDA, * ), B( LDB, * )
    END SUBROUTINE

end interface

contains

end module
module utils

! Various general utilities.
! Based on a code by John E. Pask, LLNL.

use constants, only: dp=>p_
implicit none
private
public upcase, lowcase, whitechar, blank, numstrings, getstring, &
    stop_error, arange, loadtxt, savetxt, newunit, assert, str

interface str
    module procedure str_int, str_real, str_real_n
end interface

contains

function upcase(s) result(t)
! Returns string 's' in uppercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z')) then
        ! if lowercase, make uppercase
        t(i:i) = char(ichar(t(i:i)) + diff)
    end if
end do
end function

function lowcase(s) result(t)
! Returns string 's' in lowercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
        ! if uppercase, make lowercase
        t(i:i) = char(ichar(t(i:i)) - diff)
    end if
end do
end function

logical function whitechar(char) ! white character
! returns .true. if char is space (32) or tab (9), .false. otherwise
character, intent(in) :: char
if (iachar(char) == 32 .or. iachar(char) == 9) then
    whitechar = .true.
else
    whitechar = .false.
end if
end function

logical function blank(string)
! Returns true if string contains only white characters
character(*), intent(in) :: string
integer :: i
do i = 1, len(string)
    if (.not. whitechar(string(i:i))) exit
end do
blank = (i>len(string))
end function

integer function numstrings(s) result(n)
! Returns number of substrings contained in input string 's' delimited
! by white space.
character(*), intent(in) :: s    ! input string
character(len(s)+2) :: t         ! temporary string to facilitate analysis
integer :: i
t = " " // s // " "
n = 0
do i = 1, len(t)-1
    if (whitechar(t(i:i)) .and. .not. whitechar(t(i+1:i+1))) n = n + 1
end do
end function

!--------------------------------------------------------------------------------------------------!

subroutine getstring(s,is,ss)
! Returns first substring ss in string s, delimited by white space, starting at
! index is in s. If ss is found, is is set to (index of last character of ss in
! s) + 1; else is is set to 0. If is is out of range on input, routine
! terminates with is = -1.
character(*), intent(in) :: s   ! input string
integer, intent(inout) :: is    ! on input: starting index for search for ss in
                                ! s on output: (index of last character of ss in
                                ! s) + 1
character(*), intent(out) :: ss ! first substring in s, starting from index is
character(len(s)+1) :: t        ! temporary string to facilitate search
integer i, i1, i2
logical prevwhite, curwhite
if (is <= 0 .or. is > len(s)) then
    ss = ""; is = -1; return
end if
t = s // " "
if (is == 1) then
    prevwhite = .true.
else
    prevwhite = whitechar(t(is-1:is-1))
end if
i1 = 0; i2 = 0
do i = is, len(t)
    curwhite = whitechar(t(i:i))
    if (prevwhite .and. .not. curwhite) i1 = i   ! beginning of substring
    if (i1>0 .and. curwhite) then                ! end of substring
        i2 = i-1; exit
    end if
    prevwhite=curwhite
end do
if (i2 > 0) then
    ss = t(i1:i2); is = i2+1
else
    ss = ""; is = 0
end if
end subroutine

integer function newunit(unit) result(n)
! Returns lowest i/o unit number not in use (to be used in older compilers).
!
! Starting at 10 to avoid lower numbers which are sometimes reserved.
! Note: largest valid unit number may be system-dependent.
!
! Arguments
! ---------
!
! If present, the new unit will be returned into it
integer, intent(out), optional :: unit
!
! Example
! -------
!
! integer :: u
! open(newunit(u), file="log.txt", status="old")
! read(u, *) a, b
! close(u)
!
! In new compilers, just use the "newunit" keyword argument:
!
! integer :: u
! open(newunit=u, file="log.txt", status="old")
! read(u, *) a, b
! close(u)

logical inuse
integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
integer, parameter :: nmax=999  ! may be system-dependent
do n = nmin, nmax
    inquire(unit=n, opened=inuse)
    if (.not. inuse) then
        if (present(unit)) unit=n
        return
    end if
end do
call stop_error("newunit ERROR: available unit not found.")
end function

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine

subroutine loadtxt(filename, d)
! Loads a 2D array from a text file.
!
! Arguments
! ---------
!
! Filename to load the array from
character(len=*), intent(in) :: filename
! The array 'd' will be automatically allocated with the correct dimensions
real(dp), allocatable, intent(out) :: d(:, :)
!
! Example
! -------
!
! real(dp), allocatable :: data(:, :)
! call loadtxt("log.txt", data)  ! 'data' will be automatically allocated
!
! Where 'log.txt' contains for example::
!
!     1 2 3
!     2 4 6
!     8 9 10
!     11 12 13
!     ...
!
character :: c
integer :: s, ncol, nrow, ios, i
logical :: lastwhite
real(dp) :: r

open(newunit=s, file=filename, status="old")

! determine number of columns
ncol = 0
lastwhite = .true.
do
   read(s, '(a)', advance='no', iostat=ios) c
   if (ios /= 0) exit
   if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
   lastwhite = whitechar(c)
end do

rewind(s)

! determine number or rows
nrow = 0
do
   read(s, *, iostat=ios) r
   if (ios /= 0) exit
   nrow = nrow + 1
end do

rewind(s)

allocate(d(nrow, ncol))
do i = 1, nrow
    read(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine savetxt(filename, d)
! Saves a 2D array into a textfile.
!
! Arguments
! ---------
!
character(len=*), intent(in) :: filename  ! File to save the array to
real(dp), intent(in) :: d(:, :)           ! The 2D array to save
!
! Example
! -------
!
! real(dp) :: data(3, 2)
! call savetxt("log.txt", data)

integer :: s, i
open(newunit=s, file=filename, status="replace")
do i = 1, size(d, 1)
    write(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine arange(a, b, dx, u)
! Returns an array u = [a, a+dx, a+2*dx, ..., b-dx]
!
! Arguments
! ---------
!
real(dp), intent(in) :: a, b, dx
real(dp), allocatable, intent(out) :: u(:)
!
! Example
! -------
!
! real(dp), allocatable :: u(:)
! call arange(1, 5, 1, u)   ! u = [1, 2, 3, 4]
integer :: n, i
n = int((b-a) / dx)
allocate(u(n))
do i = 1, n
    u(i) = a + (i-1)*dx
end do
end subroutine

subroutine assert(condition)
! If condition == .false., it aborts the program.
!
! Arguments
! ---------
!
logical, intent(in) :: condition
!
! Example
! -------
!
! call assert(a == 5)

if (.not. condition) call stop_error("Assert failed.")
end subroutine

pure integer function str_int_len(i) result(sz)
! Returns the length of the string representation of 'i'
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, '(i0)') i
sz = len_trim(s)
end function

pure function str_int(i) result(s)
! Converts integer "i" to string
integer, intent(in) :: i
character(len=str_int_len(i)) :: s
write(s, '(i0)') i
end function

pure integer function str_real_len(r, fmt) result(sz)
! Returns the length of the string representation of 'i'
real(dp), intent(in) :: r
character(len=*), intent(in) :: fmt
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fmt) r
sz = len_trim(s)
end function

pure function str_real(r) result(s)
! Converts the real number "r" to string with 7 decimal digits.
real(dp), intent(in) :: r
character(len=*), parameter :: fmt="(f0.6)"
character(len=str_real_len(r, fmt)) :: s
write(s, fmt) r
end function

pure function str_real_n(r, n) result(s)
! Converts the real number "r" to string with 'n' decimal digits.
real(dp), intent(in) :: r
integer, intent(in) :: n
character(len=str_real_len(r, "(f0." // str_int(n) // ")")) :: s
write(s, "(f0." // str_int(n) // ")") r
end function

end module utils



module splines

! Splines are fully specified by the interpolation points, except that
! at the ends, we have the freedom to prescribe the second derivatives.
! If we know a derivative at an end (exactly), then best is to impose that.
! Otherwise, it is better to use the "consistent" end conditions: the second
  ! derivative is determined such that it is smooth at the first and last interior knots
  !(i.e., third derivatives are continuous at those two points) jargon: "not-a-knot" condition
!
! High level API: spline3, spline3ders.
! Low level API: the rest of public soubroutines.
!
! Use the high level API to obtain cubic spline fit with consistent boundary
! conditions and optionally the derivatives. Use the low level API if more fine
! grained control is needed.
!
! This module is based on a code written by John E. Pask, LLNL.
! ref: https://github.com/certik/fortran-utils/blob/master/src/splines.f90
use constants, only: dp=>p_
use lapack, only: dgesv, dgbsv
use utils, only: stop_error
implicit none
private
public spline3pars, spline3valder, iix, iixmin, iixun, iixexp, poly3, dpoly3, &
    d2poly3, spline3, spline3ders

contains

function spline3(x, y, xnew) result(ynew)
! Takes the function values 'y' on the grid 'x' and returns new values 'ynew'
! at the given grid 'xnew' using cubic splines interpolation with such
! boundary conditions so that the 2nd derivative is consistent with the
! interpolating cubic.
real(dp), intent(in) :: x(:), y(:), xnew(:)
real(dp) :: ynew(size(xnew))
real(dp) :: c(0:4, size(x)-1)
integer :: i, ip

call spline3pars(x, y, [2, 2], [0._dp, 0._dp], c) ! get spline parameters

ip = 0
do i = 1, size(xnew)
    ip = iixmin(xnew(i), x, ip)
    ynew(i) = poly3(xnew(i), c(:, ip))
end do
end function


subroutine spline3ders(x, y, xnew, ynew, dynew, d2ynew)
  ! Just like 'spline3', but also calculate 1st and 2nd derivatives
  real(dp), intent(in) :: x(:), y(:), xnew(:)
  real(dp), intent(out), optional :: ynew(:), dynew(:), d2ynew(:)
  real(dp) :: c(0:4, size(x)-1)
  integer :: i, ip

  !call spline3pars(x, y, [2, 2], [0._dp, 0._dp], c) ! get spline parameters
  call spline3pars(x, y, [1, 1], [0._dp, 0._dp], c) ! get spline parameters
  ip = 0
  do i = 1, size(xnew)
     ip = iixmin(xnew(i), x, ip)
     if (present(  ynew))   ynew(i) =   poly3(xnew(i), c(:, ip))
     if (present( dynew))  dynew(i) =  dpoly3(xnew(i), c(:, ip))
     if (present(d2ynew)) d2ynew(i) = d2poly3(xnew(i), c(:, ip))
  end do
end subroutine spline3ders


subroutine splint(x, y, c, xnew, ynew, dynew, d2ynew)
  ! Just like 'spline3ders', but (1) for scaler xnew; (2) coefficients are assumed ready; (3) assume uniform grid
  real(dp), intent(in) :: x(:), y(:), c(0:4, size(x)-1)
  real(dp),intent(in) :: xnew
  real(dp), intent(out), optional :: ynew, dynew, d2ynew
  integer :: ip

  ip = int((xnew-x(1))/(x(2)-x(1))) + 1 !assuming uniform grid
  if (present(  ynew))   ynew =   poly3(xnew, c(:, ip))
  if (present( dynew))  dynew =  dpoly3(xnew, c(:, ip))
  if (present(d2ynew)) d2ynew = d2poly3(xnew, c(:, ip))

end subroutine splint

!the following naive implementation of cubic spline in 2D is computationally expensive, which makes it useless for large scale simulations.
subroutine spline_2d(x1a,x2a,ya,m,n,c2d)
  use constants,only:p_
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: x1a(m), x2a(n), ya(m,n)
  real(p_), intent(out) :: c2d(m, 0:4, n-1)
  integer ::  i

  do  i=1,m
     call spline3pars(x2a, ya(i,:), [2, 2], [0._p_, 0._p_], c2d(i,:,:))
  enddo

end subroutine spline_2d

subroutine splint_2d(x1a,x2a,ya,c2d,x1,x2,y, dy, d2y)
  use constants,only:p_
  implicit none
  real(p_), intent(in) :: x1a(:), x2a(:), ya(:,:), c2d(:,:,:), x1, x2
  real(p_), intent(out), optional ::  y, dy, d2y
  real(p_) :: ytmp(size(x2a)), c(0:4, size(x1a)-1)
  integer :: i

  do  i=1,size(x1a)
     call splint(x2a, ya(i,:), c2d(i,:,:), x2, ytmp(i))
  enddo

  call spline3pars(x1a, ytmp, [2, 2], [0._dp, 0._dp], c) ! get spline parameters (output in c)
  call splint     (x1a, ytmp, c, x1, y, dy, d2y)
end subroutine splint_2d


subroutine spline_2d_x1x2(x1a,x2a,ya,m,n,c2d)
  use constants,only:p_
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: x1a(m), x2a(n), ya(m,n)
  real(p_), intent(out) :: c2d(n, 0:4, m-1)
  integer ::  j

  do  j=1,n
     call spline3pars(x1a, ya(:,j), [2, 2], [0._p_, 0._p_], c2d(j,:,:))
  enddo

end subroutine spline_2d_x1x2


subroutine splint_2d_x1x2(x1a,x2a,ya,c2d,x1,x2,y, dy, d2y)
  use constants,only:p_
  implicit none
  real(p_), intent(in) :: x1a(:), x2a(:), ya(:,:), c2d(:,:,:), x1, x2
  real(p_), intent(out), optional ::  y, dy, d2y
  real(p_) :: ytmp(size(x2a)), c(0:4, size(x2a)-1)
  integer :: j

  do  j=1,size(x2a)
     call splint(x1a, ya(:,j), c2d(j,:,:), x1, ytmp(j))
  enddo

  call spline3pars(x2a, ytmp, [2, 2], [0._dp, 0._dp], c) ! get spline parameters (output in c)
  call splint     (x2a, ytmp, c, x2, y, dy, d2y)
end subroutine splint_2d_x1x2



subroutine spline3pars(xi,yi,bctype,bcval,c)
! Returns parameters c defining cubic spline interpolating x-y data xi, yi, with
! boundary conditions specified by bcytpe, bcvals
real(dp), intent(in):: xi(:)        ! x values of data
real(dp), intent(in):: yi(:)        ! y values of data
integer, intent(in):: bctype(2)     ! type of boundary condition at each end:
   ! bctype(1) = type at left end, bctype(2) = type at right end.
   ! 1 = specified 2nd derivative, 2 = 2nd derivative consistent with interpolating cubic.
real(dp), intent(in):: bcval(2)     ! boundary condition values at each end:
   ! bcval(1) = value at left end, bcval(2) = value at right end
real(dp), intent(out):: c(0:,:)     ! parameters defining spline: c(i,j) = ith parameter of jth
   ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
   ! dimensions: c(0:4,1:n-1)
real(dp) As(5,2*size(c,2))             ! spline eq. matrix -- LAPACK band form
real(dp) bs(2*size(c,2))               ! spline eq. rhs vector
real(dp) cs(2*size(c,2))               ! spline eq. solution vector
real(dp) hi(size(c,2))                 ! spline intervals
real(dp) Ae(4,4)                       ! end-cubic eq. matrix
real(dp) be(4)                         ! end-cubic eq. rhs vector
real(dp) ce(4)                         ! end-cubic eq. solution vector
real(dp) xe(4),ye(4)                   ! x,y values at ends
real(dp) d2p1,d2pn                     ! 2nd derivatives at ends
real(dp) x0                            ! expansion center
real(dp) c1,c2,c3,c4                   ! expansion coefficients
integer n                              ! number of data points
integer i,j,i2
! lapack variables
integer ipiv(4),ipiv2(2*size(c,2))
real(dp) bemat(4,1),bmat(2*size(c,2),1)
integer info

! check input parameters
if (bctype(1) < 1 .or. bctype(1) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
if (bctype(2) < 1 .or. bctype(2) > 2) call stop_error("spline3pars error: bctype /= 1 or 2.")
if (size(c,1) /= 5) call stop_error("spline3pars error: size(c,1) /= 5.")
if (size(c,2) /= size(xi)-1) call stop_error("spline3pars error: size(c,2) /= size(xi)-1.")
if (size(xi) /= size(yi)) call stop_error("spline3pars error: size(xi) /= size(yi)")

! To get rid of compiler warnings:
d2p1 = 0
d2pn = 0

! initializations
n=size(xi)
do i=1,n-1
   hi(i)=xi(i+1)-xi(i)
end do

! compute interpolating-cubic 2nd derivs at ends, if required
   ! left end
if(bctype(1)==2) then
   if (n < 4) call stop_error("spline3pars error: n < 4")
   xe=xi(1:4)
   ye=yi(1:4)
   x0=xe(1) ! center at end
   do i=1,4
      do j=1,4
         Ae(i,j) = (xe(i)-x0)**(j-1)
      end do
   end do
   Ae(:,1) = 1    ! set 0^0 = 1
   be=ye; bemat(:,1)=be
   call dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
   if (info /= 0) call stop_error("spline3pars error: dgesv error.")
   ce=bemat(:,1)
   d2p1=2*ce(3)
end if
   ! right end
if(bctype(2)==2) then
   if (n < 4) call stop_error("spline3pars error: n < 4")
   xe=xi(n-3:n)
   ye=yi(n-3:n)
   x0=xe(4) ! center at end
   do i=1,4
      do j=1,4
         Ae(i,j) = (xe(i)-x0)**(j-1)
      end do
   end do
   Ae(:,1) = 1    ! set 0^0 = 1
   be=ye; bemat(:,1)=be
   call dgesv(4, 1, Ae, 4, ipiv, bemat, 4, info)
   if (info /= 0) call stop_error("spline3pars error: dgesv error.")
   ce=bemat(:,1)
   d2pn=2*ce(3)
end if

! set 2nd derivs at ends
if(bctype(1)==1) d2p1=bcval(1)
if(bctype(2)==1) d2pn=bcval(2)
!write(*,*) d2p1,d2pn

! construct spline equations -- LAPACK band form
! basis: phi1 = -(x-x_i)/h_i, phi2 = (x-x_{i+1})/h_i, phi3 = phi1^3-phi1, phi4 = phi2^3-phi2
! on interval [x_i,x_{i+1}] of length h_i = x_{i+1}-x_i
!A=0  ! full matrix
As=0
   ! left end condition
!A(1,1)=6/hi(1)**2   ! full matrix
As(4,1)=6/hi(1)**2
bs(1)=d2p1
   ! internal knot conditions
do i=2,n-1
   i2=2*(i-1)
!   A(i2,i2-1) = 1/hi(i-1)    ! full matrix ...
!   A(i2,i2)   = 2/hi(i-1)
!   A(i2,i2+1) = 2/hi(i)
!   A(i2,i2+2) = 1/hi(i)
!   A(i2+1,i2) = 1/hi(i-1)**2
!   A(i2+1,i2+1) = -1/hi(i)**2
   As(5,i2-1) = 1/hi(i-1)
   As(4,i2)   = 2/hi(i-1)
   As(3,i2+1) = 2/hi(i)
   As(2,i2+2) = 1/hi(i)
   As(5,i2)   = 1/hi(i-1)**2
   As(4,i2+1) = -1/hi(i)**2
   bs(i2) = (yi(i+1) - yi(i))/hi(i) - (yi(i) - yi(i-1))/hi(i-1)
   bs(i2+1) = 0
end do
   ! right end condition   
!A(2*(n-1),2*(n-1))=6/hi(n-1)**2 ! full matrix
As(4,2*(n-1))=6/hi(n-1)**2
bs(2*(n-1))=d2pn

! solve spline equations -- full matrix
!bmat(:,1)=bs
!call dgesv(2*(n-1), 1, A, 2*(n-1), ipiv2, bmat, 2*(n-1), info)
!if (info /= 0) call stop_error("spline3pars error: dgesv error.")
!cs=bmat(:,1)

! solve spline equations -- LAPACK band form
bmat(:,1)=bs
call dgbsv(2*(n-1), 1, 2, 1, As, 5, ipiv2, bmat, 2*(n-1), info)
if (info /= 0) call stop_error("spline3pars error: dgbsv error.")
cs=bmat(:,1)
!write(*,*) cs(1:6)
!write(*,*) cs(2*(n-1)-5:2*(n-1))

! transform to (x-x0)^(i-1) basis and return
do i=1,n-1
   ! coefficients in spline basis:
   c1=yi(i)
   c2=yi(i+1)
   c3=cs(2*i-1)
   c4=cs(2*i)
   ! coefficients in (x-x0)^(i-1) basis
   c(0,i)=xi(i)
   c(1,i)=c1
   c(2,i)=-(c1-c2+2*c3+c4)/hi(i)
   c(3,i)=3*c3/hi(i)**2
   c(4,i)=(-c3+c4)/hi(i)**3
end do
end subroutine spline3pars

!--------------------------------------------------------------------------------------------------!

subroutine spline3valder(x,xi,c,val,der)
! Returns value and 1st derivative of spline defined by knots xi and parameters c
! returned by spline3pars
real(dp), intent(in):: x            ! point at which to evaluate spline
real(dp), intent(in):: xi(:)        ! spline knots (x values of data)
real(dp), intent(in):: c(0:,:)      ! spline parameters: c(i,j) = ith parameter of jth
   ! spline polynomial, p_j = sum_{i=1}^4 c(i,j) (x-c(0,j))^(i-1), j = 1..n-1, n = # of data pts.
   ! dimensions: c(0:4,1:n-1)
real(dp), intent(out):: val         ! value of spline at x
real(dp), intent(out):: der         ! 1st derivative of spline at x
integer n                           ! number of knots
integer i1

! initialize, check input parameters
n=size(xi)
if (size(c,1) /= 5) call stop_error("spline3 error: size(c,1) /= 5.")
if (size(c,2) /= size(xi)-1) call stop_error("spline3 error: size(c,2) /= size(xi)-1.")
! find interval containing x
i1=iix(x,xi)
! return value and derivative
val=poly3(x,c(:,i1))
der=dpoly3(x,c(:,i1))
end subroutine

!--------------------------------------------------------------------------------------------------!

integer function iix(x, xi) result(i1)
! Returns index i of interval [xi(i),xi(i+1)] containing x in mesh xi,
! with intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
! Uses bisection, except if "x" lies in the first or second elements (which is
! often the case)
real(dp), intent(in) :: x            ! target value
real(dp), intent(in) :: xi(:)        ! mesh, xi(i) < xi(i+1)
integer n                            ! number of mesh points
integer i2, ic

n = size(xi)
i1 = 1
if (n < 2) then
    call stop_error("error in iix: n < 2")
elseif (n == 2) then
    i1 = 1
elseif (n == 3) then
    if (x <= xi(2)) then ! first element
        i1 = 1
    else
        i1 = 2
    end if
elseif (x <= xi(1)) then ! left end
    i1 = 1
elseif (x <= xi(2)) then ! first element
    i1 = 1
elseif (x <= xi(3)) then ! second element
    i1 = 2
elseif (x >= xi(n)) then  ! right end
    i1 = n-1
else
    ! bisection: xi(i1) <= x < xi(i2)
    i1 = 3; i2 = n
    do
        if (i2 - i1 == 1) exit
        ic = i1 + (i2 - i1)/2
        if (x >= xi(ic)) then
            i1 = ic
        else
            i2 = ic
        endif
    end do
end if
end function

integer function iixmin(x, xi, i_min) result(ip)
  ! Just like iix, but assumes that x >= xi(i_min)
  real(dp), intent(in) :: x, xi(:)
  integer, intent(in) :: i_min
  if (i_min >= 1 .and. i_min <= size(xi)-1) then
     ip = iix(x, xi(i_min:)) + i_min - 1
  else
     ip = iix(x, xi)
  end if
end function iixmin

!--------------------------------------------------------------------------------------------------!

function iixun(x,n,x1,xn)
! Returns index i of interval [x(i),x(i+1)] containing x in uniform mesh defined by
!   x(i) = x1 + (i-1)/(n-1)*(xn-x1), i = 1 .. n,
! with intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
integer iixun                       ! index i of interval [x(i),x(i+1)] containing x
real(dp), intent(in):: x            ! target value
integer, intent(in):: n             ! number of mesh points
real(dp), intent(in):: x1           ! initial point of mesh
real(dp), intent(in):: xn           ! final point of mesh
integer i

! compute index
i=int((x-x1)/(xn-x1)*(n-1))+1
! reset if ouside 1..n
if (i<1) i=1
if (i>n-1) i=n-1
iixun=i
end function

!--------------------------------------------------------------------------------------------------!

function iixexp(x,n,x1,alpha,beta)
! Returns index i of interval [x(i),x(i+1)] containing x in exponential mesh defined by
!   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
! and intervals indexed by left-most points.
! N.B.: x outside [x1,xn] are indexed to nearest end.
integer iixexp                      ! index i of interval [x(i),x(i+1)] containing x
real(dp), intent(in):: x            ! target value
integer, intent(in):: n             ! number of mesh points
real(dp), intent(in):: x1           ! initial point of mesh
real(dp), intent(in):: alpha        ! mesh parameter:
!   x(i) = x1 + alpha [ exp(beta(i-1)) - 1 ], i = 1 .. n,
! where alpha = (x(n) - x(1))/[ exp(beta(n-1)) - 1 ],
! beta = log(r)/(n-2), r = (x(n)-x(n-1))/(x(2)-x(1)) = ratio of last to first interval,
real(dp), intent(in):: beta         ! mesh parameter
integer i

! compute index
i=int(log((x-x1)/alpha + 1)/beta) + 1
! reset if outside 1..n
if (i<1) i=1
if (i>n-1) i=n-1
iixexp=i
end function

!--------------------------------------------------------------------------------------------------!

function poly3(x,c)
! returns value of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) poly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
poly3=c(1)+c(2)*dx+c(3)*dx**2+c(4)*dx**3
end function

!--------------------------------------------------------------------------------------------------!

function dpoly3(x,c)
! returns 1st derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dpoly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
dpoly3=c(2)+2*c(3)*dx+3*c(4)*dx**2
end function

!--------------------------------------------------------------------------------------------------!

function d2poly3(x,c)
! returns 2nd derivative of polynomial \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) d2poly3
real(dp), intent(in):: x      ! point at which to evaluate polynomial
real(dp), intent(in):: c(0:)  ! coefficients: poly = \sum_{i=1}^4 c(i) (x-c(0))^(i-1)
real(dp) dx
dx=x-c(0)
d2poly3=2*c(3)+6*c(4)*dx
end function

end module splines



!subroutines from numerical reciple

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      use constants,only:p_
      implicit none
      INTEGER n,NMAX
      REAL(p_) yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=800)
      INTEGER i,k
      REAL(p_) p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+&
     &1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
     &u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      use constants,only:p_
      implicit none
      INTEGER n
      REAL(p_) x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL(p_) a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h** &
     &2)/6.
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


            SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      use constants,only:p_
      implicit none
      INTEGER m,n,NN
      REAL(p_) x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=200)
!CU    USES spline
      INTEGER j,k
      REAL(p_) y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.d30,1.d30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      use constants,only:p_
      implicit none
      INTEGER m,n,NN
      REAL(p_) x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=200)
!CU    USES spline,splint
      INTEGER j,k
      REAL(p_) y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.d30,1.d30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END
!C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
module interpolate_module
  use constants,only:p_, one
  implicit none
  private
  public linear_2d_interpolate_kernel, linear_2d_interpolate, linear_2d_interpolate0, linear_1d_interpolate, &
       & linear_1d_interpolate_nonuniform, locate
contains


  
  pure subroutine locate(m,x,dx,xval,i) !uniform grid is assumed
    use constants,only:p_
    implicit none
    integer, intent(in) :: m
    real(p_),intent(in) :: x(m),  dx,  xval
    integer,intent(out) :: i

    i = floor(one+(xval-x(1))/dx) 
    if(i==0) i=1
    if(i==m) i=m-1

  end subroutine locate

 pure subroutine linear_2d_interpolate(nx,nz,xarray,zarray,psi,x,z,psival)  !uniform xarray and zarray are assumed
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in):: xarray(nx),zarray(nz),psi(nx,nz)
    real(p_),intent(in):: x,z
    real(p_),intent(out)::psival
    real(p_):: dx,dz,t1,t2,slope
    integer:: i,j ,ii,jj

    dx = xarray(2)-xarray(1)
    i = floor(one+(x-xarray(1))/dx)

    dz = zarray(2)-zarray(1)
    j = floor(one+(z-zarray(1))/dz)

    if(i.ge.nx) i=nx-1
    if(j.ge.nz) j=nz-1
    if(i.le.1) i=1 
    if(j.le.1) j=1 
    slope = (psi(i+1,j)-psi(i,j))/dx
    t1 = psi(i,j)+slope*(x-xarray(i))
    slope = (psi(i+1,j+1)-psi(i,j+1))/dx
    t2 = psi(i,j+1)+slope*(x-xarray(i))
    slope = (t2-t1)/dz
    psival = t1+slope*(z-zarray(j))
  end subroutine linear_2d_interpolate
  
 pure subroutine linear_2d_interpolate0(nx,nz,xarray,zarray,dx,dz,psi,x,z,i,j,psival)  !uniform xarray and zarray are assumed
    use constants,only: p_, one
    implicit none
    integer,intent(in) :: nx,nz
    real(p_),intent(in) :: xarray(nx), zarray(nz), dx, dz, psi(nx,nz), x, z
    integer,intent(in) :: i,j
    real(p_),intent(out) :: psival
    real(p_) :: t1,t2,slope

    slope = (psi(i+1,j)-psi(i,j))/dx
    t1 = psi(i,j)+slope*(x-xarray(i))
    slope = (psi(i+1,j+1)-psi(i,j+1))/dx
    t2 = psi(i,j+1)+slope*(x-xarray(i))
    slope = (t2-t1)/dz
    psival = t1+slope*(z-zarray(j))
  end subroutine linear_2d_interpolate0
  
  pure  subroutine linear_2d_interpolate_kernel(x1a,x2a,ya,x1,x2,y)
    real(p_),intent(in) :: x1a(2), x2a(2), ya(2,2), x1, x2
    real(p_),intent(out) :: y
    real(p_) :: ytmp(2),slope
    integer :: j

    do j=1,2
       slope=(ya(2,j)-ya(1,j))/(x1a(2)-x1a(1))
       ytmp(j)=ya(1,j)+slope*(x1-x1a(1))
    enddo
    slope=(ytmp(2)-ytmp(1))/(x2a(2)-x2a(1))
    y=ytmp(1)+slope*(x2-x2a(1))

  end subroutine linear_2d_interpolate_kernel


  pure subroutine linear_1d_interpolate(n, x, y, xval, yval)
    use constants, only: p_, one
    implicit none
    integer, intent(in) :: n
    real(p_), intent(in) :: x(n), y(n), xval
    real(p_), intent(out) :: yval
    real(p_) :: dx, slope
    integer :: i

    dx = x(2)-x(1)
    i = floor(one+(xval-x(1))/dx) !uniform xarray is assumed

    if(i<=1) then
       yval = y(1)
       return
    endif

    if(i>=n) then
       yval = y(n)
       return
    endif

    slope = (y(i+1)-y(i))/(x(i+1)-x(i))
    yval = y(i)+slope*(xval-x(i))
  end subroutine linear_1d_interpolate



 pure subroutine linear_1d_interpolate_nonuniform(n,x,y,xval,yval)  !non-uniform x array
    use constants,only:p_, one
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(in):: xval
    real(p_),intent(out):: yval
    real(p_):: slope
    integer:: i

    !dx=x(2)-x(1)
    !i=floor(one+(xval-x(1))/dx) !this for uniform x, otherwise we need to call location() subroutine to locate xval
    if(xval.ge.x(n) ) then
       i=n-1
    elseif(xval.le.x(1)) then
       i=1
    else
       call location(n,x,xval,i)
    endif

    slope=(y(i+1)-y(i))/(x(i+1)-x(i))
    yval=y(i)+slope*(xval-x(i))

  end subroutine linear_1d_interpolate_nonuniform


pure  subroutine location(n,x,xval,k) !use bisection method to locate xval in an array
    !return k (xval is located between x(k) and x(k+1)
    use constants,only:p_
    implicit none
    integer,intent(in):: n
    real(p_),intent(in):: x(n),xval
    integer,intent(out)::k
    integer:: kl,ku,km

!!$    if(xval.gt.x(n) ) then
!!$       write(*,*) "***warning****, x provided is not in the range"
!!$       k=n-1
!!$       return
!!$    elseif(xval.lt.x(1)) then
!!$       write(*,*) "***warning****, x provided is not in the range"
!!$       k=1
!!$       return
!!$    endif
    kl=1
    ku=n
30  if(ku-kl .gt. 1) then  !use bisection method to search location of theta
       km=(ku+kl)/2
       if((x(n).ge.x(1)).eqv.(xval.ge.x(km))) then
          kl=km
       else
          ku=km
       endif
       goto 30
    endif
    k=kl
  end subroutine location



end module interpolate_module

module math
contains

real(p_)  function random_yj(seed) result (z) !return a random number uniform distributed in [0:1)
  !linear congruental method to genereate random number
  !This corresponds to Park and Miller's choice implemented using Schrange's algorithm to aviod integer overflow
  !refer to R. Fitzpatrick's book "Computational physics, an introduction course" for the details
  use constants, only: p_
  implicit none
  integer, intent(in) :: seed
  integer, parameter :: a = 16807, m=2147483647 !m=2^31-1,c=0, this choice is called Park and Miller method
  integer, parameter :: q = 127773 !q=m/a
  integer, parameter :: r = 2836 !r=mod(m,a)
  real(p_), parameter :: RANDMAX = 2147483646._p_
  integer, save :: next = 1

  if (seed .ne.0) next=seed
  next = a*mod(next,q)-r*(next/q)
  if(next<0) next=next+m
  z=next/RANDMAX
end function random_yj



  subroutine arc_length(x_contour,z_contour,np,dl)
    !calculate the poloidal arc length between neighbour points on a contour line.
    use constants,only: p_
    implicit none
    integer,intent(in) :: np
    real(p_),intent(in) :: x_contour(np), z_contour(np) 
    real(p_),intent(out) :: dl(np-1)
    integer :: i, i_plus_one, i_plus_two, i_minus_one
    real(p_) :: x(4), z(4)

    do i=1,np-1
       i_plus_one=i+1  !i_plus_one indicates the right point
       i_minus_one=i-1 !i_minus_one indicates the left point
       i_plus_two=i+2
       if (i == 1)  i_minus_one=np-1 !deal with boundary points
       if (i == np-1) i_plus_two=2 !deal with boundary points

       x(1)=x_contour(i_minus_one)
       x(2)=x_contour(i)
       x(3)=x_contour(i_plus_one)
       x(4)=x_contour(i_plus_two)

       z(1)=z_contour(i_minus_one)
       z(2)=z_contour(i)
       z(3)=z_contour(i_plus_one)
       z(4)=z_contour(i_plus_two)

       !call arc_between_two_points(x,z,dl(i))
       !use simple formula to calculate the arc length:
       dl(i) = sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2) 
    enddo

  end subroutine arc_length

  subroutine   arc_between_two_points(x,z,dl)
    !calculate the arc length between point (x(2),z(2)) and point (x(3),z(3))
    use constants,only:p_
    use constants,only:one,two
    implicit none
    real(p_),intent(in):: x(4),z(4)
    real(p_),intent(out):: dl
    real(p_):: ds,a_length,b_length
    real(p_):: dot_a_and_ds,dot_b_and_ds, cos_tha,cos_thb,m1,m2

    !ds is the length of straight-line segment passing through (x(2),z(2)) and (x(3),z(3))
    ds=sqrt((x(3)-x(2))**2+(z(3)-z(2))**2) 

    a_length=sqrt((x(3)-x(1))**2+(z(3)-z(1))**2)
    b_length=sqrt((x(4)-x(2))**2+(z(4)-z(2))**2)
    dot_a_and_ds=(x(3)-x(1))*(x(3)-x(2)) &
         +(z(3)-z(1))*(z(3)-z(2))
    dot_b_and_ds=(x(4)-x(2))*(x(3)-x(2)) &
         +(z(4)-z(2))*(z(3)-z(2))
    cos_tha=dot_a_and_ds/(a_length*ds)
    cos_thb=dot_b_and_ds/(b_length*ds)

    m1=sqrt(one-cos_tha**2)/cos_tha 
    m2=sqrt(one-cos_thb**2)/cos_thb
    !the value of m1 and m2 should be positive for most cases
    dl=ds*(one+(two*m1**2+two*m2**2+m1*m2)/30._p_) !calculate arc length using Eq. (5.38) in  S. Jardin's book. Here I assume that the dot product of the slope is negative, so the -m1*m2 term is replaced with +abs(slope1)*abs(slope2), need checking the correctness for general case
    !dl=ds*(1._p_+0.) !use linear function to approximate the arc length
  end subroutine arc_between_two_points

  subroutine ZGETRS_wrapper(kn, matrix, ipiv, rhs_dft, solution_dft)
    !solve the field equation for a single toroidal harmonic
    use constants, only: p_
    implicit none
    integer, intent(in) :: kn !index of the toroidal harmonic
    complex(p_), intent(in) :: matrix(:,:,0:) !LU factorization of the radial coefficient matrix
    integer, intent(in) :: ipiv(:,0:) 
    complex(p_), intent(in):: rhs_dft(:)
    complex(p_), intent(out):: solution_dft(:)
    complex(p_) :: s(size(rhs_dft),1)
    integer:: info, n

    n = size(rhs_dft)
    s(:,1) = rhs_dft(:)
    call ZGETRS("N", n, 1, matrix(:,:,kn), n, IPIV(:,kn), s, n, INFO) !solve using the LU factorization computed by ZGETRF.
    !call svd_back_substitution(s(:,1),kn) !solve using the SVD computed by ZGESDD, the results agree with that given by ZGETRS.
    solution_dft=s(:,1)
  end subroutine ZGETRS_wrapper


  subroutine srcbes(biz,gam0,gam1)
    use constants,only:p_
    implicit none
    REAL(p_) :: t1,t2,biz,gam0,gam1
    !.....Calculates gamma nought and gamma 1. (Abramowitz and Stegun).
446 if (biz.gt.3.75d0) go to 148
    t1=(biz/3.75d0)**2
    t2=exp(-biz)
    gam0=t2*((((((.0045813d0*t1+.0360768d0)*t1+.2659732d0)*t1+ &
         1.2067492d0)*t1+3.0899424d0)*t1+3.5156229d0)*t1+1.d0)
    gam1=t2*biz*((((((.00032411d0*t1+.00301532d0)*t1+.02658733d0) &
         *t1+.15084934d0)*t1+.51498869d0)*t1+.87890594d0)*t1+.5d0)
    go to 149
148 t2=1.d0/sqrt(biz)
    t1=3.75d0/biz
    gam0=t2*((((((((.00392377d0*t1-.01647633d0)*t1+.02635537d0) &
         *t1-.02057706d0)*t1+.00916281d0)*t1-.00157565d0)*t1+ &
         .00225319d0)*t1+.01328592d0)*t1+.39894228d0)
    gam1=t2*((((((((-.00420059d0*t1+.01787654d0)*t1-.02895312d0) &
         *t1+.02282967d0)*t1-.01031555d0)*t1+.00163801d0)*t1- &
         .00362018d0)*t1-.03988024d0)*t1+.39894228d0)
149 continue
    return
  end subroutine srcbes

  !I have verified that gam0 obtained by calling srcbes(x,gam0,gam1) is equal to bessi0(x)*exp(-x)
  FUNCTION bessi0(x) ! the 0th modified Bessel function of the first kind.
    use constants,only:p_
    implicit none
    REAL(p_)::  bessi0,x
    REAL(p_):: ax
    real(p_):: p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0, &
         & 1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
    DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,&
         & 0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1, &
         & -0.1647633d-1,0.392377d-2/
    if (abs(x).lt.3.75) then
       y=(x/3.75)**2
       bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
    else
       ax=abs(x)
       y=3.75/ax
       bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y* &
            &(q7+y*(q8+y*q9))))))))
    endif
  END FUNCTION bessi0
  !  (C) Copr. 1986-92 Numerical Recipes Software .

  subroutine cross_product_in_cartesian(ax,ay,az,bx,by,bz,cx,cy,cz)
    use constants,only:p_
    implicit none

    real(p_),intent(in):: ax,ay,az,bx,by,bz
    real(p_),intent(out)::cx,cy,cz

    cx=ay*bz-az*by
    cy=az*bx-ax*bz
    cz=ax*by-ay*bx
  end subroutine cross_product_in_cartesian


    elemental subroutine shift_toroidal(a, range) !shift "a" into the range [0:range]
    use constants,only:p_
    implicit none
    real(p_), intent(in) :: range !is positive
    real(p_), intent(inout) :: a

!!$ a=a-floor(a/range)*range

    ! a = mod(a,range)
    ! if(a<0) a = a + range

    a = modulo(a, range)
  end subroutine shift_toroidal


  subroutine shift_to_minus_pi_positive_pi_range(a) !shift "a" into the range [-pi:pi]
    use constants,only:p_
    use constants,only: twopi,pi
    implicit none
    real(p_),intent(inout):: a
    integer:: ishift

    ishift=floor(a/twopi)
    a=a-ishift*twopi
    if(a>pi) a=a-twopi
  end subroutine shift_to_minus_pi_positive_pi_range


  subroutine one_dimensional_derivative(n,x,y,dydx)
    use constants,only:p_
    use constants,only:zero,one,two,twopi,one_half
    implicit none

    integer,intent(in):: n
    real(p_),intent(in):: x(n),y(n)
    real(p_),intent(out):: dydx(n)
    real(p_):: tmp0,dx
    integer:: j

    dx=x(2)-x(1) !uniform interval is assumed

    do j=2,n-1 !use center difference scheme for inner points
       dydx(j)=(y(j+1)-y(j-1))/(two*dx)
    enddo

    !use linear interpolation to get the value  j=n
    tmp0=(y(n)-y(n-1))/dx
    dydx(n)=two*tmp0-dydx(n-1)
    !use linear interpolation to get the value j=1
    tmp0=(y(2)-y(1))/dx
    dydx(1)=two*tmp0-dydx(2)
  end subroutine one_dimensional_derivative


  subroutine calculate_determinant_complex(a,n,d)
    use constants ,only: p_
    implicit none
    integer:: n
    complex(p_)::a(n,n),d  
    integer::info,j,ipiv(n)

    ! Compute LU Factorization
    if ( p_ == kind(1.0e1) ) then   !cgetrf and zgetrf are programs in Lapack
       call cgetrf(n,n,a,n,ipiv,info)
    else
       call zgetrf(n,n,a,n,ipiv,info)
    endif

    !  compute determinant
    if ( info .ge. 0 ) then
       d = (1.0,0.0)
       do j=1,n
          d = d*a(j,j)
          if ( ipiv(j) /= j ) d=-d
       enddo
    else
       print *," *** Error in computing determinant"
       print *," info =",info
    endif
  end subroutine calculate_determinant_complex

  !C PNP1                                                                
  !C       http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html                                                                
  !C     ..................................................................
  !C                                                                       
  !C        SUBROUTINE PNPOLY                                              
  !C                                                                       
  !C        PURPOSE                                                        
  !C           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
  !C                                                                       
  !C        USAGE                                                          
  !C           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
  !C                                                                       
  !C        DESCRIPTION OF THE PARAMETERS                                  
  !C           PX      - X-COORDINATE OF POINT IN QUESTION.                
  !C           PY      - Y-COORDINATE OF POINT IN QUESTION.                
  !C           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
  !C                     VERTICES OF POLYGON.                              
  !C           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
  !C                     VERTICES OF POLYGON.                              
  !C           N       - NUMBER OF VERTICES IN THE POLYGON.                
  !C           INOUT   - THE SIGNAL RETURNED:                              
  !C                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
  !C                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
  !C                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
  !C                                                                       
  !C        REMARKS                                                        
  !C           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
  !C           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
  !C           OPTIONALLY BE INCREASED BY 1.                               
  !C           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
  !C           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
  !C           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
  !C           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
  !C           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
  !C           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
  !C           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
  !C                                                                       
  !C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
  !C           NONE                                                        
  !C                                                                       
  !C        METHOD                                                         
  !C           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
  !C           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
  !C           POINT IS INSIDE OF THE POLYGON.                             
  !C                                                                       
  !C     ..................................................................
  !C                                                                       

  pure  SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)
    use constants,only:p_
    implicit none
    integer,intent(in):: n
    integer,intent(inout):: inout
    real(p_),intent(in):: px,py
    REAL(p_),intent(in):: XX(N),YY(N)                                    
    integer:: i,j,maxdim
    REAL(p_):: X(n),Y(n), tmp
    LOGICAL MX,MY,NX,NY                                               
    !      INTEGER O                                                         
    !      OUTPUT UNIT FOR PRINTED MESSAGES                                 
    !      DATA O/6/                                                         

    !yj    MAXDIM=200
    MAXDIM=n
    IF(N.LE.MAXDIM)GO TO 6                                            
    !      WRITE(O,7)                                                        
7   FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. RESULTS INVALID')     
    RETURN                                                            
6   DO  I=1,N                                                        
       X(I)=XX(I)-PX                                                     
       Y(I)=YY(I)-PY
    enddo
    INOUT=-1                                                          
    DO I=1,N                                                        
       J=1+MOD(I,N)                                                      
       MX=X(I).GE.0.0                                                    
       NX=X(J).GE.0.0                                                    
       MY=Y(I).GE.0.0                                                    
       NY=Y(J).GE.0.0                                                    
       IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) cycle
       IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
       INOUT=-INOUT                                                      
       cycle
!!$3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
!!$4     INOUT=0                                                           
!!$      RETURN                                                            
!!$5     INOUT=-INOUT                                                      
!!$2     CONTINUE                                                          

3      tmp=(Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))
       IF(tmp<0) then
          cycle
       elseif(tmp==0) then
          goto 4
       else
          goto 5
       endif
4      INOUT=0                                                           
       RETURN                                                            
5      INOUT=-INOUT                                                      
    enddo
  END SUBROUTINE PNPOLY


  subroutine laplace_cylindrical2d(psi, x, z, nx, nz, jphi)
    !  use math, only : partial_derivatives_2d
    use constants, only : p_, mu0
    implicit none
    integer, intent(in) :: nx, nz
    real(p_), intent(in) :: psi(nx,nz),  x(nx), z(nz)
    real(p_), intent(out) :: jphi(nx,nz)
    real(p_) :: psi_x(nx,nz), psi_z(nx,nz), psi_xx(nx,nz), psi_xz(nx,nz), psi_zx(nx,nz),psi_zz(nx,nz)
    real(p_) :: dx, dz
    integer :: i, j  

    dx = x(2) - x(1)
    dz = z(2) - z(1)
    call partial_derivatives_2d(nx,nz,x,z, psi,   psi_x,  psi_z)
    call partial_derivatives_2d(nx,nz,x,z, psi_x, psi_xx, psi_xz)
    call partial_derivatives_2d(nx,nz,x,z, psi_z, psi_zx, psi_zz)

    do i= 1, nx
       do j=1,nz
          jphi(i,j) = psi_zz(i,j) + psi_xx(i,j)  - 1/x(i)*psi_x(i,j)
          jphi(i,j) = -jphi(i,j)/(mu0*x(i)) !to toroidal current density
       enddo
    enddo

  end subroutine laplace_cylindrical2d

  subroutine partial_derivatives_2d(nx,nz,rarray,zarray,b,b_x,b_z)
    use constants,only:p_
    implicit none
    integer,intent(in):: nx,nz
    real(p_),intent(in):: rarray(nx),zarray(nz),b(nx,nz)
    real(p_),intent(out):: b_x(nx,nz),b_z(nx,nz)
    integer:: i,j,i1,j1,i2,j2

    do i=1,nx
       do j=1,nz
          i2=i+1
          i1=i-1
          j2=j+1
          j1=j-1
          if(i.eq.1) i1=i
          if(j.eq.1) j1=j
          if(i.eq.nx) i2=i
          if(j.eq.nz) j2=j
          b_x(i,j)=(b(i2,j)-b(i1,j))/(rarray(i2)-rarray(i1))
          b_z(i,j)=(b(i,j2)-b(i,j1))/(zarray(j2)-zarray(j1))
       enddo
    end do
  end subroutine partial_derivatives_2d


  subroutine arrange_lcfs(x_lcfs, z_lcfs, np, x_axis, z_axis)
    !arrange the arrays x_lcfs_new and z_lcfs_new so that the starting point (x_lcfs_new(1),z_lcfs_new(1)) is on the high-field-side of the midplane
    use constants,only:p_
    use domain_decomposition, only : myid
    implicit none
    integer, intent(in) :: np
    real(p_), intent(inout) :: x_lcfs(np), z_lcfs(np)
    real(p_), intent(in) :: x_axis, z_axis
    real(p_) :: x_lcfs_new(np), z_lcfs_new(np), dz(np)
    real(p_) :: Rout
    real(p_) :: r_major, r_minor, eps, direction
    integer :: kk(1), k, i

    !kk=minloc(x_lcfs) ! to determine the high-field side of the midplane
    !k=kk(1)
    dz = abs((z_lcfs - z_axis))
    do i=1,np
       if(x_lcfs(i) > x_axis) dz(i) = 10d20 !filtering
    enddo
    kk = minloc(dz)
    k = kk(1)

    !call major_radius_on_midplane(np,x_lcfs,z_lcfs,x_axis,z_axis,Rout)
    ! write(*,*) 'inverse aspect ratio of LCFS is (definied at low-field side) ', (Rout-x_axis)/x_axis
    !  write(*,*) 'r_axis=',x_axis, 'r_major=', r_major, 'r_minor=',r_minor
    !  eps= r_minor/r_major !standard definition
    ! write(*,*) 'inverse aspect ratio of LCFS (i.e., r_minor/r_major) is ', eps
    ! write(*,*) 'ellipticity (elongation) of LCFS is ', (maxval(z_lcfs)-minval(z_lcfs))/2._p_/r_minor
    !write(*,*) 'upper triangularity of LCFS is ', (r_major-x_lcfs(maxloc(z_lcfs)))/r_minor, &
    !          & 'lower triangularity of LCFS is ', (r_major-x_lcfs(minloc(z_lcfs)))/r_minor

    !replace one point of LCFS with the new point
    !x_lcfs(k)=Rout
    !z_lcfs(k)=z_axis

    !arrange the arrays so that (x_lcfs_new(1), z_lcfs_new(1)) corresponds to (x_lcfs(k), z_lcfs(k))
    do i = 1, np
       if(k+i-1.le.np) then
          x_lcfs_new(i)=x_lcfs(k+i-1)
          z_lcfs_new(i)=z_lcfs(k+i-1)
       else
          x_lcfs_new(i)=x_lcfs(k+i-np)
          z_lcfs_new(i)=z_lcfs(k+i-np)
       endif
    enddo

    !use x_lcfs and z_lcfs to store the new data
    x_lcfs = x_lcfs_new
    z_lcfs = z_lcfs_new

    !check wheter the direction of the sequecne (r(i),z(i)) with i increasing is clockwise or anticlockwise
    !when viewed along grad_phi direction, if clockwise, switch it to anticlockwise
    !This is achieved by using the determination of the direction matrix (a well known method in graphic theory).
    !Because the contours of Psi considered here are always convex (rather than concave) polygons,
    !we can select any vertex on the curve to calculate the direction matrix. (refer to wikipedia about the direction matrix)
    direction=(x_lcfs(2)-x_lcfs(1))*(z_lcfs(3)-z_lcfs(1))-(x_lcfs(3)-x_lcfs(1))*(z_lcfs(2)-z_lcfs(1))
    if(direction .lt. 0.) then
       if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) &
            & with i increasing is clockwise, switch it to anticlockwise'

       do i=1,np !switch it to anticlockwise
          x_lcfs(i)=x_lcfs_new(np+1-i)
          z_lcfs(i)=z_lcfs_new(np+1-i)
       enddo
    else if (direction .gt. 0.) then
       if(myid.eq.0) write(*,*) 'the direction of sequency (x_lcfs(i,j),z_lcfs(i,j)) with i increasing is anticlockwise'
    else
       stop 'the three vertex (points) used in calculating the direction matrix is collinear'
    endif

  end subroutine arrange_lcfs



  subroutine major_radius_on_midplane(mpol,rs,zs,r_axis,z_axis,Rout)
    !given a flux surface, this subroutine determines the major radius of the point on the low/high-field side of the middle plane
    use constants,only:p_
    use interpolate_module,only: linear_1d_interpolate_nonuniform

    implicit none
    integer,intent(in):: mpol
    real(p_),intent(in):: rs(mpol), zs(mpol),r_axis,z_axis
    real(p_),intent(out):: Rout !the major radius of the point on the low/high-field-side of the mid-plane
    integer:: i,j,k1,k2,n
    real(p_):: r_select(mpol),z_select(mpol)
    !real(p_),dimension(:),allocatable:: x,z,tmp_y2
    real(p_):: tmp

    n=1
    do i=1,mpol-1 
       !if(rs(i).gt.r_axis) then !select the low-field side (i.e., the part with r larger than r_axis) of a flux surface
       if(rs(i).lt.r_axis) then !select the high-field side (i.e., the part with r less than r_axis) of a flux surface
          r_select(n)=rs(i)
          z_select(n)=zs(i)
          n=n+1
       endif
    enddo
    !  write(*,*) 'n-1= ', n-1
    !order the array according to the value of z_select
    do i=1,n-1
       do j=i+1,n-1
          if(z_select(j).le.z_select(i)) then
             !exchange the z value 
             tmp=z_select(i)
             z_select(i)=z_select(j)
             z_select(j)=tmp
             !also exchange the r value (I forgot this step in the older version, which cause a serious error)
             tmp=r_select(i)
             r_select(i)=r_select(j)
             r_select(j)=tmp
          endif
       enddo
    end do

    call linear_1d_interpolate_nonuniform(n-1,z_select,r_select,z_axis,Rout)  

  end subroutine major_radius_on_midplane


  subroutine solver_toroidal_mode_number_parallel(matrix, IPIV, source_dft, potential_dft)
    use constants, only : p_
    use control_parameters, only : nh_min, nh_max
    use domain_decomposition, only : TCLR, ntube, grid_comm
    use mpi
    implicit none
    complex(p_), intent(in) :: matrix(:,:,0:) !LU factorization of the radial coefficient matrix
    integer, intent(in) :: ipiv(:,0:) 
    complex(p_), intent(in) :: source_dft(0:,:)
    complex(p_), intent(inout) :: potential_dft(0:,:)
    logical, save :: is_first=.true.
    integer, save :: my_kn_start, my_kn_end, sendcount, n, nht
    integer, save, allocatable :: recvcounts(:), displace(:)
    complex(p_), allocatable  :: my_field(:,:), field(:,:)
    integer :: i, ierr, kn, ntask

    if (is_first.eqv..true.) then !needs to be genreated only for the first time and used multiple times
       is_first=.false.
       allocate(recvcounts(0:ntube-1))
       allocate(displace(0:ntube-1))
       n=size(source_dft,2) !radial grid number
       nht = nh_max - nh_min + 1 !total number of toroidal harmonics
       if(nht <= ntube) then
          my_kn_start = TCLR + nh_min
          my_kn_end = my_kn_start
          if(TCLR>nht-1) then
             my_kn_start = 0
             my_kn_end = -1
          endif
          recvcounts(0:nht-1)=1*n
          recvcounts(nht:)= 0
          do i=0,nht-1
             displace(i)=i*1*n
          enddo
       else
          ntask = nht/ntube
          my_kn_start = TCLR*ntask + nh_min
          my_kn_end = my_kn_start + ntask-1
          if(TCLR .eq. (ntube-1)) my_kn_end = nht-1 !the last process handles all the remainder part
          recvcounts(:)=ntask*n
          recvcounts(ntube-1)= (nht-ntask*(ntube-1))*n  !last process contains additional elements
          do i=0,ntube-1
             displace(i)=i*ntask*n
          enddo
       endif
       sendcount= (my_kn_end - my_kn_start +1)*n
    endif

    allocate(my_field(n, my_kn_start:my_kn_end))
    allocate(field(n, nh_min:nh_max))

    do kn = my_kn_start, my_kn_end !for non-negative toroidal mode number
       call ZGETRS_wrapper(kn, matrix, IPIV, source_dft(kn,:), my_field(:,kn)) 
    enddo

    call MPI_Allgatherv(my_field, sendcount, MPI_Double_Complex, &
         &      field, recvcounts, displace, MPI_Double_Complex, Grid_comm, ierr)

    do kn = nh_min, nh_max !re-arange the gathered data
       potential_dft(kn,:) = field(:,kn)
    enddo

  end subroutine solver_toroidal_mode_number_parallel
  
end module math





subroutine read_and_process_equilibrium()
  use constants,only: p_, zero,one,two,three,four,five,twopi 
  use poloidal_flux_2d,only:xarray,zarray,nx,nz,psi,psi_gradient,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx
  use radial_module,only:psi_axis,psi_lcfs,npsi,psi_1d,fpsi,qpsi,pfn_npsi,tfn_npsi,baxis,r_axis,z_axis 
  use domain_decomposition,only:myid,numprocs
  implicit none
  integer :: i, u

  call read_gfile()
  
  allocate(psi_x(nx,nz))
  allocate(psi_z(nx,nz))
  allocate(psi_xx(nx,nz))
  allocate(psi_zz(nx,nz))
  allocate(psi_xz(nx,nz))
  allocate(psi_zx(nx,nz))
  allocate(psi_gradient(nx,nz))
  
  !call calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  call calculate_poloidal_flux_partial_derivatives2(nx,nz,xarray,zarray,psi,psi_x,psi_z,psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)

  npsi=nx !nx in g-file is also used to define the number of radial grid points.
  allocate(psi_1d(npsi))
  do i=1,npsi
     psi_1d(i)=psi_axis+(psi_lcfs-psi_axis)/(npsi-1)*(i-1) 
  enddo

  allocate(tfn_npsi(npsi))
  allocate(pfn_npsi(npsi))

  call calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
  pfn_npsi=(psi_1d-psi_1d(1))/(psi_1d(npsi)-psi_1d(1))

  if(myid==0) then
     open(newunit=u,file='q1.txt')
     do i=1,npsi
        write(u,'(10ES18.4)') psi_1d(i), pfn_npsi(i), tfn_npsi(i), qpsi(i)
     enddo
     close(u)
  endif
  
end subroutine read_and_process_equilibrium


subroutine read_gfile()
  use constants, only : p_, two
  use poloidal_flux_2d, only: xarray, zarray, nx, nz, psi !as output
  use radial_module,only: psi_axis,psi_lcfs ,fpsi,ffprime,fprime,qpsi,baxis,r_axis,z_axis, pressure, pprime, sign_bphi
  use radial_module, only : 
  use boundary,only: nlim,rlim,zlim,np_lcfs, x_lcfs,z_lcfs !as output
  use math, only : arrange_lcfs
  use domain_decomposition, only : myid
  implicit none
  character(100) :: gfile_name
  logical :: reverse_tf,reverse_ip
  character(len=100) :: format1, format2, format3
  character(len=8) :: ntitle(5), vid
  integer:: neq,ipestg
  real(p_)::xdim, zdim, rleft, zmid
  real(p_)::  btorus, current, major_R
  real(p_):: dumaraya4(4),dumarayb5(5)
  integer:: i,j, u

  namelist /magnetic_configuration/gfile_name, reverse_tf, reverse_ip
  open(newunit=u, file='input.nmlt')
  read(u, magnetic_configuration)
  close(u)
  if(myid==0)  write(*,magnetic_configuration)

  format1='(5e16.9)'
  format2='(6a8, 3i4)'
  format3='(2i5)'
  !open and read in eqdsk file (refer to G EQDSK.pdf (or weqdsku.f in onetwo) for the gfile format)
  open(newunit=neq,file=gfile_name,status='old')
  ipestg = 4
  read (neq, format2) (ntitle(i), i=1,5), vid, ipestg, nx, nz
  allocate(psi(nx,nz))
  allocate(qpsi(nx))
  allocate(fpsi(nx))
  allocate(ffprime(nx))
  allocate(pressure(nx))
  allocate(pprime(nx))
  read (neq, format1) xdim, zdim, major_R, rleft, zmid
  read (neq, format1) r_axis, z_axis, psi_axis, psi_lcfs, btorus
  read (neq, format1) current, dumaraya4
  read (neq, format1) dumarayb5
  read (neq ,format1) (fpsi(i), i=1,nx)
  read (neq ,format1) (pressure(i), i=1,nx)
  read (neq ,format1) (ffprime(i), i=1,nx)
  read (neq ,format1) (pprime(i), i=1,nx)
  read (neq ,format1) ((psi(i,j), i=1,nx), j=1,nz)
  read (neq ,format1) (qpsi(i), i=1,nx)
  read (neq ,format3) np_lcfs, nlim
  allocate(x_lcfs(np_lcfs))
  allocate(z_lcfs(np_lcfs))
  !    nlim=56
  allocate(rlim(nlim))
  allocate(zlim(nlim))
  read (neq ,format1) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
  read (neq ,format1) (rlim(i), zlim(i), i=1,nlim)
  close(neq)


!!$    open(newunit=u,file='hl/first_wall.txt')
!!$    do i=1,nlim
!!$       read(u,*) rlim(i), zlim(i)
!!$    enddo
!!$    close(u)

  !to verify that I have read the eqdsk file correctly, I write the data read to a new file called 'tmp.gfile'.
  !After the program finished, I compare the file 'tmp.gfile' with the original file using diff command
  !the output of diff command indicates that the two files are idential, which shows I have read the eqdsk file correctly
  if(myid==0) then
!!$         open(newunit=neq,file='tmp.gfile')
!!$         write (neq, format2) (ntitle(i), i=1,5), vid, ipestg, nx, nz
!!$         write (neq, format1) xdim, zdim, major_R, rleft, zmid
!!$         write (neq, format1) r_axis, z_axis, psi_axis, psi_lcfs, btorus
!!$         write (neq, format1) current, dumaraya4
!!$         write (neq, format1) dumarayb5
!!$         write (neq ,format1) (fpsi(i), i=1,nx)
!!$         write (neq ,format1) (pressure(i), i=1,nx)
!!$         write (neq ,format1) (ffprime(i), i=1,nx)
!!$         write (neq ,format1) (pprime(i), i=1,nx)
!!$         write (neq ,format1) ((psi(i,j), i=1,nx), j=1,nz)
!!$         write (neq ,format1) (qpsi(i), i=1,nx)
!!$         write (neq ,format3) np_lcfs, nlim
!!$         write (neq ,format1) (x_lcfs(i), z_lcfs(i), i=1,np_lcfs)
!!$         write (neq ,format1) (rlim(i), zlim(i), i=1,nlim)
!!$         close(neq)
  endif


  !Somtimes, I alter some quantities (e.g. reverse the toroidal magnetic field, or increase the pressure by a constant),in this case, the out gfile is different from the original one
!!$  do i=1,nx
!!$     pressure(i)=pressure(i)+0.5*pressure(1) !increase the presure
!!$  enddo

  if(reverse_tf.eqv..true.)  then
     fpsi=-fpsi !revert the toroidal magnetic field
  endif
  baxis=fpsi(1)/r_axis

  if(reverse_ip.eqv..true.) then
     psi=-psi !revert direction of the torodial current
     psi_axis=-psi_axis
     psi_lcfs=-psi_lcfs
     ffprime=-ffprime
     pprime=-pprime
  endif

  if(fpsi(1)>0) then 
     sign_bphi=1
     if(myid==0) write(*,*) 'bphi>0'
  else
     sign_bphi=-1
     if(myid==0) write(*,*) 'bphi<0'
  endif

  if(myid == 0) then
     open(newunit=u,file='lcfs.txt')
     do i = 1, np_lcfs
        write(u,*) x_lcfs(i), z_lcfs(i)
     enddo
     close(u)
  endif

  allocate(xarray(nx))
  allocate(zarray(nz))

  do i=1,nx !construct the X array
     xarray(i)=rleft+xdim/(nx-1)*(i-1)
  enddo

  do j=1,nz !construct the Z array
     zarray(j)=(zmid-zdim/two)+zdim/(nz-1)*(j-1)
  enddo

  allocate(fprime(nx))
  fprime=ffprime/fpsi

  call arrange_lcfs(x_lcfs, z_lcfs, np_lcfs, r_axis, z_axis) !change starting point and direction (counter-clockwise)

  if(myid==0) then
     open(newunit=u,file='lcfs2.txt')
     do i = 1, np_lcfs
        write(u,*) x_lcfs(i), z_lcfs(i)
     enddo
     close(u)

     write(*,*) 'nR, nZ =', nx, nz
     write(*,*) 'rleft=',rleft, 'rright=', rleft+xdim,'zlow=',zmid-zdim/2._p_,'zupp=',zmid+zdim/2._p_
     write(*,*) 'magnetic location (meter) (r,z)=', r_axis, z_axis
     write(*,*) 'baxis (Tesla)=', baxis
     write(*,*) 'rcenter=',major_R, 'vacuum magnetic field at rcenter=',btorus
     !write(*,*)  'total current in all TF coils (MegaAmpere)=', btorus*major_R/(2d-7)/10**6 !Ampere's circuital law
     write(*,*) 'psi_axis=',psi_axis,'psi_lcfs=',psi_lcfs
     write(*,*) 'np_lcfs=',np_lcfs
     !write(*,*) 'x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)=', x_lcfs(1),z_lcfs(1),x_lcfs(np_lcfs),z_lcfs(np_lcfs)
     write(*,*) 'Cyclotron angular frequency of Deuterium ion at magnetic axis (10^6 rad/s)=', &
          & fpsi(1)/r_axis*1.6022d-19/3.3452d-27/1.d6

     if(psi_lcfs>psi_axis) then
        write(*,*) 'Iphi<0'
     else
        write(*,*) 'Iphi>0'
     endif
     block !find out the direction of the troidal current
       use math, only : laplace_cylindrical2d
       real(p_), allocatable:: jphi(:,:)
       real(p_) :: dx, dz, s

       allocate(jphi(nx,nz))
       call laplace_cylindrical2d(psi(1:nx,1:nz), xarray, zarray, nx, nz, jphi)
       if (sum(jphi(nx/5:nx/2,nz/5:nz/2))>0) then
          write(*,*) 'Jphi>0'
       else
          write(*,*) 'Jphi<0'
       endif
       dx = xarray(2) -xarray(1)
       dz = zarray(2) -zarray(1)
       if(myid==0) then
          open(newunit=u,file='psi_and_jphi.txt')
          s = 0
          do i=1,nx
             do j =1,nz
                write(u,'(4e16.5)') xarray(i), zarray(j), psi(i,j), jphi(i,j)
                s = s + jphi(i,j)*dx*dz
             enddo
             write(u,*)
          enddo
          write(*,*) 'total plasma current calculated from the poloidal magnetic flux (kA)=', s/1000.
          write(*,*) "total plasma current given in g-file (kA)= ", current/1000.
          close(u)
       endif

     end block
     open(newunit=u,file='limiter.txt')
     do i=1,nlim
        write(u,*) rlim(i), zlim(i)
     enddo
     close(u)
     write(*,*) 'wall r_min=',minval(rlim), 'wall r_max=',maxval(rlim), 'wall z_min=',minval(zlim), 'wall z_max=',maxval(zlim)
  endif

end subroutine read_gfile


subroutine calculate_poloidal_flux_partial_derivatives(nx,nz,xarray,zarray,psi,psi_x,psi_z,&
     & psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  !subroutine calculate_poloidal_flux_gradient(nx,nz,xarray,zarray,psi,psi_gradient)
  use constants,only:p_
  use constants,only:one,two,twopi
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: psi(nx,nz)
  real(p_),intent(in):: xarray(nx),zarray(nz)
  real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz),psi_xx(nx,nz),psi_zz(nx,nz),psi_xz(nx,nz),psi_zx(nx,nz)
  real(p_),intent(out):: psi_gradient(nx,nz)

  integer:: i,j,i1,i2,j1,j2

  !first-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_x(i,j)=(psi(i2,j)-psi(i1,j))/(xarray(i2)-xarray(i1))
        psi_z(i,j)=(psi(i,j2)-psi(i,j1))/(zarray(j2)-zarray(j1))
     enddo
  enddo

  !second-order partial derivatives
  do i=1,nx
     do j=1,nz
        i2=i+1
        i1=i-1
        j2=j+1
        j1=j-1
        if(i.eq.1) i1=i
        if(j.eq.1) j1=j
        if(i.eq.nx) i2=i
        if(j.eq.nz) j2=j
        psi_xx(i,j)=(psi_x(i2,j)-psi_x(i1,j))/(xarray(i2)-xarray(i1))
        psi_zz(i,j)=(psi_z(i,j2)-psi_z(i,j1))/(zarray(j2)-zarray(j1))
        psi_xz(i,j)=(psi_x(i,j2)-psi_x(i,j1))/(zarray(j2)-zarray(j1))
        psi_zx(i,j)=(psi_z(i2,j)-psi_z(i1,j))/(xarray(i2)-xarray(i1))
     enddo
  enddo


  psi_gradient=sqrt(psi_x**2+psi_z**2)


end subroutine calculate_poloidal_flux_partial_derivatives


subroutine calculate_poloidal_flux_partial_derivatives2(nx,nz,xarray,zarray,psi,psi_x,psi_z,&
     & psi_xx,psi_zz,psi_xz,psi_zx,psi_gradient)
  use constants,only:p_, one,two,twopi
  use splines, only: spline3ders
  implicit none
  integer,intent(in)::nx,nz
  real(p_),intent(in):: xarray(nx),zarray(nz), psi(nx,nz)
  real(p_),intent(out):: psi_x(nx,nz),psi_z(nx,nz),psi_xx(nx,nz),psi_zz(nx,nz),psi_xz(nx,nz),psi_zx(nx,nz)
  real(p_),intent(out):: psi_gradient(nx,nz)
  real(p_) :: tmp(nx,nz)
  integer :: i, j
  
  do i = 1, nx
     call spline3ders(zarray, psi(i,:), zarray, dynew=psi_z(i,:), d2ynew=psi_zz(i,:))
  enddo

  do j = 1, nz
     call spline3ders(xarray, psi(:,j), xarray, dynew=psi_x(:,j), d2ynew=psi_xx(:,j))
  enddo

  do i = 1, nx
     call spline3ders(zarray, psi_x(i,:), zarray, dynew=psi_xz(i,:))
  enddo

!!$  psi_zx = psi_xz
  do j = 1, nz
     call spline3ders(xarray, psi_z(:,j), xarray, dynew=psi_zx(:,j))
  enddo

  
 psi_gradient=sqrt(psi_x**2+psi_z**2)

end subroutine calculate_poloidal_flux_partial_derivatives2

subroutine calculate_tfn(npsi,psi_1d,qpsi,tfn_npsi)
  !calculate the toroidal magnetic flux
  use constants,only:p_
  use constants,only: two,twopi
  use domain_decomposition,only:myid
  implicit none
  integer,intent(in):: npsi
  real(p_),intent(in):: psi_1d(npsi),qpsi(npsi)
  real(p_),intent(out):: tfn_npsi(npsi)
  real(p_):: dpsi,tf_npsi(npsi),total_tf
  integer:: j

  dpsi=(psi_1d(npsi)-psi_1d(1))/(npsi-1)
  tf_npsi(1)=0._p_
  do j=2,npsi 
     tf_npsi(j)=tf_npsi(j-1)+(qpsi(j)+qpsi(j-1))/two*twopi*dpsi  !using the formula dtf=q*dpf=q*d(pf_gs)*twopi
  enddo

  if(myid.eq.0) then
     open(213,file='pf_tf.txt')
     do j=1,npsi
        write(213,*) psi_1d(j),tf_npsi(j)
     enddo
     close(213)
  endif
  total_tf=tf_npsi(npsi)
  tfn_npsi= tf_npsi/total_tf !normalized toroidal magnetic flux
end subroutine calculate_tfn


subroutine draw_rect_region(nx,nz,r_1d,z_1d)
  use constants,only:p_
  implicit none
  integer,intent(in):: nx,nz
  real(p_),intent(in):: r_1d(nx),z_1d(nz)
  integer:: i,j

  open(12,file='rectangular')
  do j=1,nz
     write(12,*) r_1d(1),z_1d(j)
  enddo
  do i=1,nx
     write(12,*) r_1d(i),z_1d(nz)
  enddo
  do j=1,nz
     write(12,*) r_1d(nx),z_1d(nz-j+1)
  enddo

  do i=1,nx
     write(12,*) r_1d(nx-i+1),z_1d(1)
  enddo
  close(12)
end subroutine draw_rect_region


module magnetic_field
contains
  !all independent and dependent variables are in SI units
  pure real(p_) function psi_func(x,z) result(f)!SI units, poloidal magnetic flux function
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi,x,z,f)  
  end function psi_func

  pure real(p_) function pfn_func(x,z) result(f)!normalized poloidal magnetic flux
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi
    use radial_module,only:psi_axis,psi_lcfs
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z

    call linear_2d_interpolate(nx,nz,xarray,zarray,psi,x,z,f)  
    f=(f-psi_axis)/(psi_lcfs-psi_axis)
  end function pfn_func

  pure real(p_) function psi_r_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_x
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_x,x,z,psi_r_func)  
  end function psi_r_func


  pure real(p_) function psi_z_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_z !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in):: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_z,x,z,psi_z_func)  
  end function psi_z_func


  pure real(p_) function psi_rr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xx !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_xx,x,z,psi_rr_func)  
  end function psi_rr_func

  pure real(p_) function psi_zz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zz !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_zz,x,z,psi_zz_func)  
  end function psi_zz_func


  pure real(p_) function psi_rz_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_xz 
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_xz,x,z,psi_rz_func)  
  end function psi_rz_func


  pure real(p_) function psi_zr_func(x,z)
    use constants,only:p_
    use poloidal_flux_2d,only: nx,nz,xarray,zarray,psi_zx 
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: x,z
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_zx,x,z,psi_zr_func)  
  end function psi_zr_func

  pure real(p_) function psi_gradient_func(xval,zval) result (f)
    use constants,only:p_
    use poloidal_flux_2d,only: xarray,zarray,psi_gradient,y2a_gradient,nx,nz
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: xval,zval
    call linear_2d_interpolate(nx,nz,xarray,zarray,psi_gradient,xval,zval,f)
  end function psi_gradient_func


  pure real(p_) function br(r,z) !R component of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in) :: r,z
    br=-psi_z_func(r,z)/r
  end function br

  pure real(p_) function bz(r,z) !Z component of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in) :: r,z
    bz=psi_r_func(r,z)/r
  end function bz

  pure real(p_) function bphi(r,z) !phi component of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    bphi=g_func(psi_func(r,z))/r
  end function bphi


  pure real(p_) function b(r,z) !strength of magnetic field
    use constants,only:p_
    implicit none
    real(p_), intent(in) :: r,z
    b=sqrt(br(r,z)**2+bz(r,z)**2+bphi(r,z)**2)
  end function b

  pure real(p_) function br_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    br_r=psi_z_func(r,z)/r**2-psi_rz_func(r,z)/r
  end function br_r


  pure real(p_) function br_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    br_z=-psi_zz_func(r,z)/r
  end function br_z

  pure real(p_) function bz_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_) , intent(in) :: r,z
    bz_r=-psi_r_func(r,z)/r**2+psi_rr_func(r,z)/r
  end function bz_r


  pure real(p_) function bz_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    bz_z=psi_rz_func(r,z)/r
  end function bz_z


  pure real(p_) function bphi_r(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    bphi_r=gprime(psi_func(r,z))*psi_r_func(r,z)/r-g_func(psi_func(r,z))/r**2
  end function bphi_r

  pure real(p_) function bphi_z(r,z) !partial derivative of component of magnetic field
    use constants,only:p_
    use constants,only:one,zero
    implicit none
    real(p_), intent(in):: r,z
    bphi_z=gprime(psi_func(r,z))*psi_z_func(r,z)/r
  end function bphi_z

  pure real(p_) function b_r(r,z) !partial derivative of  magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in) :: r,z
    b_r=-b(r,z)/r+one/(b(r,z)*r*r)*(psi_r_func(r,z)*psi_rr_func(r,z) &
         & +psi_z_func(r,z)*psi_rz_func(r,z)+g_func(psi_func(r,z))*g_r(r,z))
  end function b_r

  pure real(p_) function b_z(r,z) !partial derivative of magnetic field
    use constants,only:p_
    use constants,only:one
    implicit none
    real(p_), intent(in):: r,z
    b_z=one/(b(r,z)*r*r)*(psi_r_func(r,z)*psi_rz_func(r,z)+&
         & psi_z_func(r,z)*psi_zz_func(r,z)+g_func(psi_func(r,z))*g_z(r,z))
  end function b_z

  pure real(p_) function b_phi(r,z) result(f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    f=0._p_
  end function b_phi

  pure real(p_) function unitbr_r(r,z)  result(f)!partial derivative of the component of magnetic unit vector
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    bval=b(r,z)
    f=(br_r(r,z)*bval-b_r(r,z)*br(r,z))/bval**2
  end function unitbr_r

  pure real(p_) function unitbr_z(r,z) result(f)
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval

    bval=b(r,z)
    f=(br_z(r,z)*bval-b_z(r,z)*br(r,z))/bval**2
  end function unitbr_z


  pure real(p_) function unitbr_phi(r,z) result(f)
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval

    !bval=b(r,z)
    !f=(0*bval-b_phi(r,z)*br(r,z))/bval**2
    f=0.
  end function unitbr_phi


  pure real(p_) function unitbz_r(r,z) result(f) 
    use constants,only:p_
    use normalizing,only:Ln
    implicit none
    real(p_), intent(in):: r,z
    real(p_)::bval
    bval=b(r,z)
    f=(bz_r(r,z)*bval-b_r(r,z)*bz(r,z))/bval**2
  end function unitbz_r

  pure real(p_) function unitbz_z(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval

    bval=b(r,z)
    f=(bz_z(r,z)*bval-b_z(r,z)*bz(r,z))/bval**2
  end function unitbz_z

  pure real(p_) function unitbz_phi(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    !bval=b(r,z)
    !f=(0*bval-b_phi(r,z)*bz(r,z))/bval**2
    f=0.
  end function unitbz_phi


  pure real(p_) function unitbphi_r(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    bval=b(r,z)
    f=(bphi_r(r,z)*bval-b_r(r,z)*bphi(r,z))/bval**2
  end function unitbphi_r


  pure real(p_) function unitbphi_z(r,z) result (f)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: bval
    bval=b(r,z)
    f=(bphi_z(r,z)*bval-b_z(r,z)*bphi(r,z))/bval**2
  end function unitbphi_z


  pure real(p_) function g_r(r,z)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: psival
    psival=psi_func(r,z)
    g_r=gprime(psival)*psi_r_func(r,z)
  end function g_r

  pure real(p_) function g_z(r,z)
    use constants,only:p_
    implicit none
    real(p_), intent(in):: r,z
    real(p_):: psival
    psival=psi_func(r,z)
    g_z=gprime(psival)*psi_z_func(r,z)
  end function g_z

  pure real(p_) function g_func(psival) result(f) !all quantities are in S.I units.
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fpsi
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_), intent(in) :: psival
    
    call linear_1d_interpolate(npsi,psi_1d,fpsi,psival,f)  
  end function g_func

  pure real(p_) function gprime(psival) result(f)
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,fprime
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_), intent(in):: psival
    call linear_1d_interpolate(npsi,psi_1d,fprime,psival,f)  
  end function gprime

  pure  real(p_) function qfunc(pfn) result(f) !safety factor, with correct sign
    use constants, only: p_
    use radial_module, only: npsi, pfn_npsi, q_with_sign
    use interpolate_module, only: linear_1d_interpolate
    implicit none
    real(p_), intent(in) :: pfn
    
    call linear_1d_interpolate(npsi, pfn_npsi, q_with_sign, pfn, f)  
  end function qfunc

  pure real(p_) function qfunc0(psival) result(f) !safety factor
    !not used in computing the orbits, only used to do analytical estimation of some quantities, such as bounce frequency
    use constants,only:p_
    use radial_module,only:npsi,psi_1d,qpsi
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_), intent(in):: psival
    call linear_1d_interpolate(npsi,psi_1d,qpsi,psival,f)  
  end function qfunc0


  pure real(p_) function radcor_as_func_of_pfn(pfn) result (z)
    use constants,only: p_
    implicit none
    real(p_), intent(in):: pfn
    z=pfn
  end function radcor_as_func_of_pfn

  pure real(p_) function tfn_func_pfn(pfn0) result(z)
    use constants, only: p_
    use radial_module, only: npsi, pfn_npsi, tfn_npsi
    use interpolate_module, only: linear_1d_interpolate
    implicit none
    real(p_), intent(in) :: pfn0

    call linear_1d_interpolate(npsi, pfn_npsi, tfn_npsi, pfn0, z)  
  end function tfn_func_pfn


end module magnetic_field
module magnetic_coordinates
  use constants, only: p_
  implicit none
  save
  integer :: mpol, nrad, mpol2, mtor 
  !(mpol, mtor, nrad) is number of gridpoints along the (poloidal, toroidal, radial) direction, respectively.
  !mpol2 is number of poloidal gridpoints for perturbed field
  real(p_), dimension(:), allocatable :: xgrid, ygrid, zgrid  !(radial, toroidal, poloidal) grid
  real(p_) :: pfn_inner, pfn_bdry, radial_width
  real(p_) :: dtheta, dradcor, dtor
  real(p_) :: xlow, xupp
  integer :: nsegment !computational region is 1/nsegment of the full torus
  real(p_) :: toroidal_range !toroidal_range=twopi/nsegment

  integer :: itheta0 !poloidal index corresponding to theta=0
  real(p_) :: vol, GSpsi_prime !dGSpsi/dx, x is the normalized poloidal magnetic flux
  real(p_) :: sign_of_jacobian, sign_of_GSpsi_prime
  real(p_), dimension(:,:), allocatable :: r_mc, z_mc !SI units
  real(p_), dimension(:,:), allocatable :: tor_shift_mc, jacobian, abs_jacobian, qhat
  real(p_), dimension(:), allocatable :: tor_shift_left_bdry_minus_one 
  real(p_), dimension(:), allocatable :: jacobian_av
  !GSpsi_array is Grad-Shafranov poloidal flux in SI units, i.e. poiloidal_magnetic_flux/twopi
  real(p_), dimension(:), allocatable :: GSpsi_array, pfn, tfn, minor_r_array, minor_r_prime_array

  real(p_), dimension(:,:),allocatable :: dl_mc, &
       & rth, zth, rpsi, zpsi, &
       & grad_psi_r, grad_psi_z, &
       & grad_theta_r, grad_theta_z, &
       & grad_alpha_r, grad_alpha_z, grad_alpha_phi, &
       & grad_psi, grad_alpha, grad_theta,&
       & grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, grad_alpha_dot_grad_theta
  real(p_), dimension(:), allocatable :: grad_alpha_r_left_bdry_minus_one, grad_alpha_z_left_bdry_minus_one

end  module magnetic_coordinates

module contour_mod
contains
  subroutine find_contour(psival, x_contour, z_contour)
    ! Given a value of the poloidal flux, psival,
    ! this subroutine find the magnetic surface corresponding to psival
    use constants, only: p_
    use boundary, only: x_lcfs, z_lcfs, np_lcfs
    use radial_module, only: r_axis, z_axis
    implicit none
    real(p_), intent(in) :: psival
    real(p_), intent(out) :: x_contour(np_lcfs), z_contour(np_lcfs)
    real(p_), parameter:: xacc = 1.0d-6 !tolerance used in bi-section root-finder
    real(p_), parameter:: huge = 1d30
    real(p_) :: x1, x2, z1, z2, slope(np_lcfs), slope2(np_lcfs)
    integer :: i

    do i=1,np_lcfs
       if(x_lcfs(i)-r_axis .ne. 0._p_) then 
          slope(i)= (z_lcfs(i)-z_axis)/(x_lcfs(i)-r_axis) !the slope for function Z=Z(X)
       else
          slope(i) = huge !I use compiler option that catches all erroneous arithmetic operation, I need to avoid dividing by zero
       endif
       if(z_lcfs(i)-z_axis .ne. 0._p_) then
          slope2(i)=(x_lcfs(i)-r_axis)/(z_lcfs(i)-z_axis) !the slope for function X=X(Z)
       else
          slope2(i) = huge
       endif
    enddo

    do i = 1, np_lcfs-1 
       if(abs(slope(i)).le.1.0_p_) then !use Z=Z(X) function, the reason that I switch between using function X=X(Z) and Z=Z(X) is to aviod large slope.
          x1 = r_axis
          x2 = x_lcfs(i) !+0.01 !shift left a little to gurrantee that the range is enough for a root to lie in
          x_contour(i)=rtbis(one_dim_psi_func,x1,x2,xacc,r_axis,z_axis,slope(i),psival)
          z_contour(i)=zfunc(r_axis,z_axis,slope(i),x_contour(i))
       else !switch to using X=X(Z) function
          z1=z_axis
          z2=z_lcfs(i)
          z_contour(i)=rtbis(one_dim_psi_func2,z1,z2,xacc,r_axis,z_axis,slope2(i),psival)
          x_contour(i)=xfunc(r_axis,z_axis,slope2(i),z_contour(i)) 
       endif
    enddo

    x_contour(np_lcfs) = x_contour(1) !i=1 and i=np_lcfs are identical
    z_contour(np_lcfs) = z_contour(1) 

  end subroutine find_contour


  function one_dim_psi_func(r_axis,z_axis,slope,psival,x) 
    !poloidal flux as a function of x on a straight line with slope "slope" in poloidal plane
    use constants, only: p_
    use magnetic_field, only: psi_func
    implicit none
    real(p_) :: one_dim_psi_func,x,r_axis,z_axis,slope,psival
    one_dim_psi_func = psi_func(x,zfunc(r_axis,z_axis,slope,x))-psival
  end function one_dim_psi_func


  function zfunc(r_axis,z_axis,slope,x) !straight line Z=Z(x) with slope "slope" in poloidal plane starting from the location of magnetic axis
    use constants,only:p_
    implicit none
    real(p_):: zfunc,x,r_axis,z_axis,slope
    zfunc=z_axis+slope*(x-r_axis)
  end function zfunc


  function one_dim_psi_func2(r_axis,z_axis,slope,psival,z) result(fun_val)
    !poloidal flux as a function of z on a straight line with slope "slope" in poloidal plane
    use constants,only:p_
    use magnetic_field, only : psi_func
    implicit none
    real(p_):: fun_val,z
    real(p_):: r_axis,z_axis,slope,psival
    fun_val = psi_func(xfunc(r_axis,z_axis,slope,z),z)-psival
  end function one_dim_psi_func2


  function xfunc(r_axis,z_axis,slope,z) !straight line X=X(Z) with slope "slope" in poloidal plane starting from the location of magnetic axis
    use constants,only:p_
    implicit none
    real(p_):: xfunc, z
    real(p_)::r_axis,z_axis,slope
    xfunc=r_axis+slope*(z-z_axis)
  end function xfunc


  FUNCTION rtbis(func,x1,x2,xacc,xmaxis,zmaxis,slope,psival)
    !find a root of func by using the bisection method
    use constants,only: p_
    implicit none
    INTEGER, parameter :: JMAX=40
    REAL(p_) rtbis, x1, x2, func, xacc, xmaxis, zmaxis, slope, psival
    external :: func
    INTEGER j
    REAL(p_) dx,f,fmid,xmid
    fmid=func(xmaxis,zmaxis,slope,psival,x2)
    f=   func(xmaxis,zmaxis,slope,psival,x1)
    !      write(*,*) 'f1=', f, 'f2=',fmid
    if(f*fmid.ge.0.) stop 'root must be bracketed in rtbis'

    if(f.lt.0.)then
       rtbis=x1
       dx=x2-x1
    else
       rtbis=x2
       dx=x1-x2
    endif
    do  j=1,JMAX
       dx=dx*.5
       xmid=rtbis+dx
       fmid=func(xmaxis,zmaxis,slope,psival,xmid)
       if(fmid.le.0.)rtbis=xmid
       if(abs(dx).lt.xacc .or. fmid.eq.0.) return
    enddo
    stop 'too many bisections in rtbis'
  end FUNCTION rtbis

end module contour_mod


subroutine construct_magnetic_coordinates()
  use constants, only : p_, zero, twopi,pi
  use control_parameters, only: diagnosis
  use boundary, only: x_lcfs, z_lcfs, np_lcfs
  use radial_module,only: r_axis,z_axis, minor_a
  use magnetic_coordinates, only: mpol, mtor, nsegment, nrad, pfn, tfn, GSpsi_array, &
       & r_mc,z_mc, & !as output
       & zgrid, dtheta, & !as output
       & minor_r_array,minor_r_prime_array, itheta0, dl_mc, & !output
       & toroidal_range, dtor, ygrid, vol !output
  use domain_decomposition, only: myid
  use contour_mod, only: find_contour
  use splines, only: spline3ders
  use math, only:  one_dimensional_derivative, arc_length
  implicit none
  real(p_), allocatable :: r_mag_surf0(:,:), z_mag_surf0(:,:)
  integer :: i, j

  call choose_radial_grids()

  allocate(r_mag_surf0(np_lcfs, nrad))
  allocate(z_mag_surf0(np_lcfs, nrad))
  
  do j =1, nrad
     call find_contour(GSpsi_array(j),r_mag_surf0(:,j),z_mag_surf0(:,j))
  enddo

  if(myid.eq.0) call diagnostic1()
  minor_a = (maxval(x_lcfs) - minval(x_lcfs))/2
  if(myid==0) write(*,*) 'minor_a=', minor_a
  allocate(minor_r_array(nrad))
  allocate(minor_r_prime_array(nrad))
  do j=1,nrad
     !minor_r_array(j) = r_axis - r_mag_surf0(1,j) !defined on the high-field-side midplane
     minor_r_array(j) = (maxval(r_mag_surf0(:,j))-minval(r_mag_surf0(:,j)))/2
  enddo
  call one_dimensional_derivative(nrad, pfn, minor_r_array, minor_r_prime_array) !to be used as interpolating table
  !call spline3ders(pfn(:), minor_r_array(:), pfn(:), dynew=minor_r_prime_array(:))

  dtheta = twopi/(mpol-1) !poloidal grid spacing for equilibrium
  allocate(zgrid(mpol))
  zgrid = [ (-pi+dtheta*(i-1), i = 1, mpol) ] !equilibrium theta grid
  itheta0 = (mpol+1)/2 !poloidal index corresponding to theta=0 (mpol is assumed odd)

  allocate(r_mc(mpol,nrad))
  allocate(z_mc(mpol,nrad))

  do j = 1, nrad
     call construct_poloidal_coordinate(r_mag_surf0(:,j),z_mag_surf0(:,j),np_lcfs, &
          & mpol,zgrid, r_mc(:,j), z_mc(:,j))
  enddo

  allocate(dl_mc(mpol-1,nrad))
  do j=1,nrad
     call arc_length(r_mc(:,j), z_mc(:,j), mpol, dl_mc(:,j))
  enddo

  call calculate_metric() 

  toroidal_range = twopi/nsegment
  dtor = toroidal_range/mtor
  allocate(ygrid(mtor+1)) 
  ygrid = [ (zero + dtor*(i-1), i = 1, mtor+1) ]

  call plasma_volume_of_computational_region(vol)
  if ((myid==0) .and. (diagnosis .eqv. .true.)) call diagnostic2()
  if((myid==0) .and. (diagnosis.eqv..true.)) call diagnostic3()
  deallocate(r_mag_surf0, z_mag_surf0)

contains

  subroutine diagnostic1()
    integer:: u
    open(newunit=u,file='mag_surf_shape0.txt')
    do j=1,nrad
       do i=1,np_lcfs
          write(u,*) r_mag_surf0(i,j), z_mag_surf0(i,j)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)
  end subroutine diagnostic1

  subroutine diagnostic2()
    integer:: u
    open(newunit=u,file='theta_line.txt')
    do i=1,mpol
       do j=1,nrad
          write(u,*) r_mc(i,j),z_mc(i,j),zgrid(i)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)

    open(newunit=u,file='mag_surf_shape.txt')
    do j=1,nrad
       do i=1,mpol
          write(u,*) r_mc(i,j),z_mc(i,j)
       enddo
       write(u,*);  write(u,*)
       write(u,*);  write(u,*)
    enddo
    close(u)
  end subroutine diagnostic2

  subroutine diagnostic3()
    integer:: u
    open(newunit=u,file='minor_r.txt')
    do j=1,nrad
       write(u,*) pfn(j), minor_r_array(j),  minor_r_prime_array(j)
    enddo
    close(u)
  end subroutine diagnostic3

end subroutine construct_magnetic_coordinates



subroutine choose_radial_grids()
  use constants, only: p_
  use radial_module, only: psi_axis, psi_lcfs, j_fixed, radcor_fixed
  use magnetic_field, only: radcor_as_func_of_pfn, tfn_func_pfn
  use magnetic_coordinates, only: nrad, pfn_inner, pfn_bdry, &
       & pfn, tfn, GSpsi_array,GSpsi_prime, xgrid,dradcor, & !output
       & xlow, xupp, & !as output
       & radial_width !as output
  implicit none
  real(p_) :: dpfn
  integer :: j

  allocate(pfn(nrad))
  allocate(tfn(nrad))  
  allocate(GSpsi_array(nrad))
  allocate(xgrid(nrad))
  dpfn=(pfn_bdry-pfn_inner)/real(nrad-1)
  do j = 1, nrad !select some flux surfaces (labeld by GSpsi_array)
     pfn(j) = pfn_inner +dpfn*(j-1)
     GSpsi_array(j)=psi_axis+pfn(j)*(psi_lcfs-psi_axis) !GSpsi=Aphi*R, this array is used in finding magnetic surfaces
     xgrid(j)=radcor_as_func_of_pfn(pfn(j))
  enddo

  xlow = xgrid(1)
  xupp = xgrid(nrad)
  
  radial_width = xgrid(nrad) - xgrid(1)
  dradcor = xgrid(2) - xgrid(1) !radial grid interval
  GSpsi_prime = psi_lcfs - psi_axis !dGSpsi/dx, x is the normalized poloidal magnetic flux

  j_fixed = nrad/2
  radcor_fixed = xgrid(j_fixed) !the radcor of the center of computational region, used in flux tube model

  do j = 1, nrad
     tfn(j) = tfn_func_pfn(pfn(j))
  enddo

  
end subroutine choose_radial_grids


subroutine construct_poloidal_coordinate(r_old,z_old,mpol_old, mpol, theta_new, r_new, z_new) !on a magnetic surface
  use constants,only: p_, two,pi,twopi
  use control_parameters,only: poloidal_angle_type
  use magnetic_field, only : psi_gradient_func, b
  use math, only: arc_length
  implicit none
  integer, intent(in) :: mpol_old, mpol
  real(p_), intent(in) :: r_old(mpol_old), z_old(mpol_old)
  real(p_), intent(in) :: theta_new(mpol)
  real(p_), intent(out) :: r_new(mpol), z_new(mpol)
  real(p_) :: theta_old(mpol_old), dl(mpol_old-1)
  real(p_) :: rmid, zmid, y2(mpol_old)
  integer :: i

  call arc_length(r_old, z_old, mpol_old, dl)
  
  theta_old(1)=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,mpol_old
        theta_old(i) = theta_old(i-1) + dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'straight-field-line') then
     do i=2,mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)/(rmid*psi_gradient_func(rmid,zmid)) !straight-field-line poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'Boozer') then
     do i=2,mpol_old
        rmid=0.5_p_*(r_old(i-1)+r_old(i))
        zmid=0.5_p_*(z_old(i-1)+z_old(i))
        theta_old(i)=theta_old(i-1)+dl(i-1)*b(rmid,zmid)**2*rmid/psi_gradient_func(rmid,zmid)
     enddo
  else
     stop 'please choose poloidal angle type among equal-arc/equal-volume/Boozer/straight-field-line'
  endif

  !normalize and shift to the range [-pi:pi]
  theta_old = theta_old*twopi/theta_old(mpol_old) - pi 

  !interpolate R  to theta_new gridpoints
  call spline(theta_old, r_old, mpol_old, 2.d30, 2.d30, y2) !prepare the second order derivative 
  do i = 2, mpol-1
     call splint(theta_old,r_old,y2,mpol_old,theta_new(i), r_new(i))
  enddo

  !interpolate Z  to theta_new gridpoints
  call spline(theta_old, z_old, mpol_old, 2.d30, 2.d30, y2) !prepare the second order derivative 
  do i = 2, mpol-1
     call splint(theta_old,z_old,y2,mpol_old,theta_new(i), z_new(i))
  enddo

  r_new(1) = r_old(1) !ending points are not included in the above interpolation
  z_new(1) = z_old(1)
  r_new(mpol) = r_new(1)
  z_new(mpol) = z_new(1)

end subroutine construct_poloidal_coordinate


subroutine plasma_volume_of_computational_region(vol)
  use constants,only: p_, two, twopi
  use domain_decomposition,only: myid
  use magnetic_coordinates,only: mpol,nrad, dradcor, dtheta, toroidal_range, jacobian
  implicit none
  real(p_), intent(out) :: vol
  real(p_) :: dv(mpol, nrad)
  integer :: i, j

  do i=1,mpol
     do j=1,nrad
        dv(i,j)=abs(jacobian(i,j))*dradcor*dtheta*toroidal_range
     enddo
  enddo

  vol=0._p_
  do i=1,mpol-1
     do j=1,nrad
        vol = vol + dv(i,j)
     enddo
  enddo

  if(myid==0) write(*,*) 'Volume of the toroidal wedge (m^3):',  vol
  if(myid==0) write(*,*) 'Volume of full torus of the computational domain (m^3):', vol*twopi/toroidal_range

end subroutine plasma_volume_of_computational_region


module misc

contains

subroutine magnetic_coordinates_to_cylindrical_coordinates(theta,radcor,r,z) !given (theta,radcor), return (R,Z)
  use constants, only: p_
  use magnetic_coordinates, only: mpol,nrad,r_mc,z_mc,zgrid,xgrid
  use interpolate_module, only: linear_2d_interpolate
  implicit none
  real(p_), intent(in) :: theta, radcor
  real(p_), intent(out) :: r,z
  call linear_2d_interpolate(mpol,nrad, zgrid,xgrid, r_mc, theta,radcor,R)  !uniform 1darray is assumed
  call linear_2d_interpolate(mpol,nrad, zgrid,xgrid, z_mc, theta,radcor,Z)  !uniform 1darray is assumed
end subroutine magnetic_coordinates_to_cylindrical_coordinates


function abs_jacobian_func(theta, radcor) result (z) !!used in generating non-uniformly distributed random numbers that satisfied a probability density function proportional to the |jacobian|
  use constants, only: p_
  use magnetic_coordinates, only: mpol,nrad,zgrid,xgrid, abs_jacobian !as input
  use interpolate_module, only: linear_2d_interpolate
  implicit none
  real(p_) :: radcor,theta,z
  !real(p_) :: abs_jacobian_min,abs_jacobian_max

  call linear_2d_interpolate(mpol,nrad, zgrid, xgrid, abs_jacobian,theta,radcor,z)  !uniform 1d array is assumed

!!$  do i=1,mpol
!!$     do j=1,nrad
!!$        jshift=j
!!$        jac1(i,j)=jacobian0(i,jshift)
!!$     enddo
!!$  enddo
!!$  abs_jacobian_max=maxval(jac1)
  !abs_jacobian_max=maxval(jacobian0(:,1:nrad))
  !  abs_jacobian_min=minval(jac1)

  !  z=(z-abs_jacobian_min)/(abs_jacobian_max-abs_jacobian_min) !turns out to be wrong! a subtle bug which I spent much time in finding, the shift is wrong. This bug is due to my misunderstanding about the rejection method, not due to programming mistakes.
  !z=z/abs_jacobian_max
end function abs_jacobian_func


subroutine calculate_dvol(m, dvol) !one dtheta2 cell can include multiple equilibrium dtheta
  use constants, only: p_
  use domain_decomposition, only: ipol_eq, dtheta2
  use magnetic_coordinates, only: nrad, mpol, nrad, abs_jacobian, dradcor, dtor
  implicit none
  integer, intent(in) :: m
  real(p_), intent(out), allocatable :: dvol(:)
  real(p_) :: jac(-m:mpol+m, nrad), jac0
  integer :: i, i1, i2, s, j, jeq

  allocate(dvol(nrad))

  jac(1:mpol, :) = abs_jacobian(:,:)
  do i = 0, -m, -1
     jac(i,:)= abs_jacobian(mpol-1+i,:)
  enddo
  do i = mpol+1, mpol+m
     jac(i,:)= abs_jacobian(i+1-mpol,:)
  enddo

  if(mod(m,2) == 0) then
     s = m/2
     i1 = ipol_eq - s
     i2 = ipol_eq + s
     do j = 1, nrad
        jeq = j
        jac0 = (sum(jac(i1+1:i2-1, jeq)) + 0.5*(jac(i1, jeq) + jac(i2, jeq)))/m
        dvol(j) = jac0*dradcor*dtheta2*dtor
     enddo
  else
     s = (m-1)/2
     i1 = ipol_eq - s
     i2 = ipol_eq + s
     do j = 1, nrad
        jeq = j
        jac0 = sum(jac(i1:i2, jeq))/m
        dvol(j) = jac0*dradcor*dtheta2*dtor
     enddo
  endif
end subroutine calculate_dvol

end module misc


module table_in_mc !used in computing guiding-center drift
  use constants, only: p_
  implicit none
  save
  real(p_), dimension(:,:), allocatable :: br_mc, bz_mc, bphi_mc, &
       & bp_mc, b_mc, bdgxcgy, & !here bdgxcgy is B0_dot_grad_x_cross_grad_y. Similar names for other cooordinate combination.
       & bdgxcgz, bdgycgz, &
       & w1,w2,w3,w4,w5,w5p,w6,w7,w8,w8p,w9,w10,w12,w13

contains

  subroutine prepare_table_in_mc()
    use constants, only: one
    use domain_decomposition, only: myid
    use magnetic_coordinates, only: mpol, nrad, r_mc,z_mc, jacobian, xgrid, &
         &  GSpsi_prime, grad_psi, grad_alpha, &
         & grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, grad_alpha_dot_grad_theta, &
         & grad_psi_r, grad_psi_z, grad_theta_r, grad_theta_z, &
         & grad_alpha_r, grad_alpha_z
    use magnetic_field, only: br, bz, bphi, b, &
         & b_r, b_z, unitbr_z, unitbz_r, unitbphi_r,unitbphi_z
    implicit none

    real(p_) :: brval,bzval,bphival,bval
    real(p_) :: b_rval,b_zval
    real(p_) :: dradial_dr_func,dradial_dz_func !function names
    real(p_) :: dradial_dr_func2,dradial_dz_func2 !function names
    real(p_) :: dtheta_dr_func,dtheta_dz_func !function names
    real(p_) :: ddelta_dr_func,ddelta_dz_func,ddelta_dr_lsf_midplane_twopi,ddelta_dz_lsf_midplane_twopi !function names
    real(p_) :: dradial_dr_val,dradial_dz_val
    real(p_) :: dtheta_dr_val,dtheta_dz_val
    real(p_) :: ddelta_dr_val,ddelta_dz_val
    real(p_) :: unitbr,unitbz,unitbphi
    real(p_) :: curl_unitb_rcomp,curl_unitb_zcomp,curl_unitb_phicomp
    real(p_) :: unitb_dot_curl_unitb
    real(p_) :: r,z,radcor
    real(p_) :: grad_psi_val,grad_alpha_val, grad_psi_dot_grad_alpha_val
    real(p_) :: grad_psi_dot_grad_theta_val,grad_alpha_dot_grad_theta_val
    real(p_) :: dalpha_dr_val,dalpha_dz_val,dalpha_dphi_val
    integer :: i,j

    allocate(w1(mpol,nrad))
    allocate(w2(mpol,nrad))
    allocate(w3(mpol,nrad))
    allocate(w4(mpol,nrad))
    allocate(w5(mpol,nrad))
    allocate(w5p(mpol,nrad))
    allocate(w6(mpol,nrad))
    allocate(w7(mpol,nrad))
    allocate(w8(mpol,nrad))
    allocate(w8p(mpol,nrad))
    allocate(w9(mpol,nrad))
    allocate(w10(mpol,nrad))
    allocate(w12(mpol,nrad))
    allocate(w13(mpol,nrad))
    allocate(bdgxcgz(mpol,nrad))
    allocate(bdgycgz(mpol,nrad))
    allocate(bdgxcgy(mpol,nrad))
    allocate(b_mc(mpol,nrad))
    allocate(br_mc(mpol,nrad))
    allocate(bz_mc(mpol,nrad))
    allocate(bphi_mc(mpol,nrad))
    allocate(bp_mc(mpol,nrad))


    do j=1,nrad
       radcor=xgrid(j)
       do i=1,mpol
          r=r_mc(i,j)
          z=z_mc(i,j)
          brval=br(r,z)
          bzval=bz(r,z)
          bphival=bphi(r,z)
          bval=sqrt(brval**2+bzval**2+bphival**2)
          b_zval=b_z(r,z)
          b_rval=b_r(r,z)
          unitbr=brval/bval
          unitbz=bzval/bval
          unitbphi=bphival/bval
          curl_unitb_rcomp=-unitbphi_z(r,z)
          curl_unitb_phicomp=unitbr_z(r,z)-unitbz_r(r,z)
          curl_unitb_zcomp=unitbphi_r(r,z)+unitbphi/r

          unitb_dot_curl_unitb=unitbr*curl_unitb_rcomp +unitbphi*curl_unitb_phicomp&
               & +unitbz*curl_unitb_zcomp

          dradial_dr_val = grad_psi_r(i,j)
          dradial_dz_val = grad_psi_z(i,j) 
          dtheta_dr_val = grad_theta_r(i,j)
          dtheta_dz_val = grad_theta_z(i,j)
          ddelta_dr_val = -grad_alpha_r(i,j)
          ddelta_dz_val = -grad_alpha_z(i,j)

          dalpha_dr_val = grad_alpha_r(i,j)
          dalpha_dz_val = grad_alpha_z(i,j)
          dalpha_dphi_val=one

          !        if(i.eq.mpol)    ddelta_dr_val=ddelta_dr_lsf_midplane_twopi(r) !at theta=twopi cut
          !        if(i.eq.mpol)    ddelta_dz_val=ddelta_dz_lsf_midplane_twopi(r) !at theta=twopi cut

!!$!write(*,*) 'i,j=',dradial_dr_val, dradial_dr_func2(r,z),dradial_dz_val, dradial_dz_func2(r,z)
          w1(i,j)=unitb_dot_curl_unitb/bval
          w2(i,j)=unitbr*dtheta_dr_val+unitbz*dtheta_dz_val
          !      w2(i,j)=-(psi_lcfs-psi_axis)/(bval*jacobian(i,j))
          w3(i,j)=curl_unitb_rcomp/bval*dradial_dr_val+curl_unitb_zcomp/bval*dradial_dz_val
          w4(i,j)=curl_unitb_rcomp/bval*dtheta_dr_val+curl_unitb_zcomp/bval*dtheta_dz_val
          w5(i,j)=curl_unitb_phicomp/(bval*r)-curl_unitb_rcomp/bval*ddelta_dr_val-curl_unitb_zcomp/bval*ddelta_dz_val
          w5p(i,j)=curl_unitb_phicomp/(bval*r)
          w6(i,j)=(bphival*b_zval)*dradial_dr_val+(-bphival*b_rval)*dradial_dz_val
          w6(i,j)=w6(i,j)/bval**2
          w7(i,j)=(bphival*b_zval)*dtheta_dr_val+(-bphival*b_rval)*dtheta_dz_val
          w7(i,j)=w7(i,j)/bval**2
          w8(i,j)=(bzval*b_rval-brval*b_zval)/r-(bphival*b_zval)*ddelta_dr_val-(-bphival*b_rval)*ddelta_dz_val
          w8(i,j)=w8(i,j)/bval**2
          w8p(i,j)=(bzval*b_rval-brval*b_zval)/(r*bval**2)
          w9(i,j)=unitbr*b_rval+unitbz*b_zval
          w10(i,j)=(b_rval*curl_unitb_rcomp+b_zval*curl_unitb_zcomp)/bval

!!$        grad_psi(i,j)=sqrt(dradial_dr_val**2+dradial_dz_val**2) !now calculated in another subroutine
!!$        grad_alpha(i,j)=sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)
!!$        grad_psi_dot_grad_alpha(i,j)=dradial_dr_val*ddelta_dr_val+dradial_dz_val*ddelta_dz_val

          b_mc(i,j)=bval
          br_mc(i,j)=brval
          bz_mc(i,j)=bzval
          bphi_mc(i,j)=bphival
          bp_mc(i,j)=sqrt(brval**2+bzval**2)

          grad_psi_val=grad_psi(i,j)
          grad_alpha_val=grad_alpha(i,j)
          grad_psi_dot_grad_alpha_val=grad_psi_dot_grad_alpha(i,j)
          grad_psi_dot_grad_theta_val=grad_psi_dot_grad_theta(i,j)
          grad_alpha_dot_grad_theta_val=grad_alpha_dot_grad_theta(i,j)

          w12(i,j)=GSpsi_prime/bval*(grad_psi_dot_grad_theta_val*grad_psi_dot_grad_alpha_val &
               & -grad_psi_val**2*grad_alpha_dot_grad_theta_val)
          w13(i,j)=GSpsi_prime/bval*(grad_psi_dot_grad_theta_val*grad_alpha_val**2 &
               & -grad_psi_dot_grad_alpha_val*grad_alpha_dot_grad_theta_val)

          bdgxcgz(i,j) = bphival*(dradial_dz_val*dtheta_dr_val-dradial_dr_val*dtheta_dz_val)
          bdgycgz(i,j) = brval*dtheta_dz_val/r &
               & + bphival*(dalpha_dz_val*dtheta_dr_val - dalpha_dr_val*dtheta_dz_val) &
               & + bzval*(-dtheta_dr_val)/r
          bdgxcgy(i,j) = bval**2/GSpsi_prime
       enddo

    enddo

    if(myid.eq.0) call diagnostic()
  contains
    subroutine diagnostic()
      character(8)::filename
      integer:: file_unit,u
      open(newunit=u,file='grad_theta_alpha.txt')
      do j=1,nrad
         radcor=xgrid(j)
         do i=1,mpol
!!$          write(u,*) r_mc(i,j), z_mc(i,j), dtheta_dr_val,dtheta_dz_val,&
!!$               & sqrt(dtheta_dz_val**2+dtheta_dr_val**2), ddelta_dr_val,  ddelta_dz_val,&
!!$               sqrt(ddelta_dr_val**2+ddelta_dz_val**2), sqrt(one/r**2+ddelta_dr_val**2+ddelta_dz_val**2)
!!$          write(u,*) r_mc(i,j), z_mc(i,j),jacobian(i,j),&
!!$               & -w2(i,j)

            !write(u,*) r_mc(i,j), z_mc(i,j),w12(i,j),w13(i,j),w5(i,j)
            write(u,*) r_mc(i,j), z_mc(i,j), w5(i,j), w8(i,j)

            !               & -b_mc(i,j)/GSpsi_prime*w2(i,j)


         enddo
         write(u,*) 
      enddo
      close(u)
!!$  write(filename,'(a4,i4.4)') 'metr',myid
!!$  file_unit=myid+311
!!$  open(file_unit,file=filename)
!!$  do i=1,mpol
!!$     do j=1,nrad
!!$        write(file_unit,'(2i8.4,3(1pe14.5))')  i,j,grad_alpha(i,j),grad_psi_dot_grad_alpha(i,j), grad_psi(i,j)
!!$     enddo
!!$     write(file_unit,*)
!!$  enddo
!!$  close(file_unit)

      ! write(*,*) 'min max w1-4=', maxval(w1),minval(w1),maxval(w2),minval(w2),maxval(w3),minval(w3),maxval(w4),minval(w4)
      !write(*,*) 'min max w5-8=', maxval(w5),minval(w5),maxval(w6),minval(w6),maxval(w7),minval(w7),maxval(w8),minval(w8)
      !write(*,*) 'min max w9-10=', maxval(w9),minval(w9),maxval(w10),minval(w10)

    end subroutine diagnostic

  end subroutine prepare_table_in_mc

end module table_in_mc



module func_in_mc
contains
  function b_mc_func(theta, radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: b_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,b_mc,theta,radcor,f)
  end function b_mc_func

  function br_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: br_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,br_mc,theta,radcor,f)
  end function br_mc_func

  function bz_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: bz_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate

    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bz_mc,theta,radcor,f)
  end function bz_mc_func


  function bphi_mc_func(theta,radcor) result(f)
    use constants,only:p_
    use table_in_mc,only: bphi_mc
    use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_)::theta,radcor,f
    call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bphi_mc,theta,radcor,f)
  end function bphi_mc_func

pure real(p_)  function tor_shift_func(theta, radcor) result(z)
    use constants, only: p_
    use magnetic_coordinates, only: mpol, nrad, zgrid, xgrid, tor_shift_mc
    use interpolate_module, only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: theta, radcor
    call linear_2d_interpolate(mpol, nrad, zgrid, xgrid, tor_shift_mc, theta, radcor, z) !interpolate in magnetic coordinates
  end function tor_shift_func


pure real(p_)   function grad_psi_func(theta,radcor) result(f)
    use constants, only: p_
    use magnetic_coordinates, only: mpol, nrad, zgrid, xgrid, grad_psi
    use interpolate_module, only: linear_2d_interpolate
    implicit none
    real(p_), intent(in) :: theta, radcor
    call linear_2d_interpolate(mpol, nrad, zgrid, xgrid, grad_psi, theta, radcor, f)
  end function grad_psi_func

  pure real(p_) function minor_r_radcor(radcor) result (z)
    use constants,only: p_, two,twopi
    use magnetic_coordinates,only: nrad, xgrid, minor_r_array
    use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_),intent(in):: radcor

    call linear_1d_interpolate(nrad, xgrid, minor_r_array, radcor, z)  
  end function minor_r_radcor

  pure real(p_) function radcor_minor_r(minor_r) result (z)
    use constants,only: p_
    use magnetic_coordinates,only: nrad,xgrid,minor_r_array
    use interpolate_module,only: linear_1d_interpolate_nonuniform
    implicit none
    real(p_), intent(in):: minor_r

    call linear_1d_interpolate_nonuniform(nrad, minor_r_array, xgrid, minor_r, z)

  end function radcor_minor_r


  pure real(p_) function minor_r_prime(radcor) result (z) !derivative of minor_r with respect to the radial coordinate
    use constants, only: p_, two
    use magnetic_coordinates, only: nrad,xgrid,minor_r_prime_array
    use interpolate_module, only: linear_1d_interpolate
    implicit none
    real(p_), intent(in) :: radcor

    call linear_1d_interpolate(nrad, xgrid, minor_r_prime_array, radcor, z)  
  end function minor_r_prime

  
end module func_in_mc

module calculate_toroidal_shift_module
contains
subroutine calculate_toroidal_shift(psival,x_contour,z_contour,np_lcfs,end_i,r,z,tor_shift) !tor_shift=q*delta. zeta=phi-tor_shift
  use constants,only: p_, one_half
  ! use poloidal_flux_2d,only: nx,nz,xarray,zarray
  use magnetic_coordinates,only:sign_of_jacobian,sign_of_GSpsi_prime, itheta0
  use radial_module,only:psi_axis,psi_lcfs
  use domain_decomposition,only:myid
  use magnetic_field, only : psi_z_func,psi_r_func, g_func !toroidal field function
  implicit none
  real(p_),intent(in):: psival
  integer,intent(in):: np_lcfs,end_i
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  real(p_),intent(in):: r,z
  real(p_),intent(out):: tor_shift

  real(p_):: x_mid,z_mid,gx0,dl, g_value
  real(p_):: pfn,mr,costh,sinth,r0
  integer:: i

  g_value=g_func(psival)

  tor_shift=0._p_
  if (end_i.lt.itheta0) then 
     do i=itheta0-1,end_i+1,-1
        x_mid=(x_contour(i)+x_contour(i-1))*one_half
        z_mid=(z_contour(i)+z_contour(i-1))*one_half
        gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
        dl=-sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
        tor_shift=tor_shift+g_value/(x_mid*gx0)*dl
     enddo
     x_mid=(x_contour(end_i+1)+r)*one_half
     z_mid=(z_contour(end_i+1)+z)*one_half
     gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=-sqrt((x_contour(end_i+1)-r)**2+(z_contour(end_i+1)-z)**2)
     tor_shift=tor_shift+g_value/(x_mid*gx0)*dl

  elseif(end_i.ge.itheta0) then
     x_contour(end_i+1)=r !replace No. end_i+1 point by the given point (r,z)
     z_contour(end_i+1)=z
     do i=itheta0, end_i
        x_mid=(x_contour(i)+x_contour(i+1))*one_half
        z_mid=(z_contour(i)+z_contour(i+1))*one_half
        gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
        dl=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
        tor_shift=tor_shift+g_value/(x_mid*gx0)*dl
     enddo
  endif
  !--for testing----analytical formula, for concentric circular configuration
!!$  r0=1.32_p_
!!$  mr=sqrt((x_contour(end_i+1)-r0)**2+(z_contour(end_i+1)-0._p_)**2)
!!$  costh=(x_contour(end_i+1)-r0)/mr
!!$  sinth=z_contour(end_i+1)/mr
!!$  pfn=(psival-psi_axis)/(psi_lcfs-psi_axis)
!!$  tor_shift=2*qfunc(pfn)*atan((r0-mr)/sqrt(r0**2-mr**2)*sinth/(costh+1))
  !----for tesing----------------------

!!$  do i=2,end_i+1
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     tor_shift=tor_shift+g_value/(x_mid*gx0(i-1))*dl(i-1)
!!$  enddo
  tor_shift=tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
end subroutine calculate_toroidal_shift


pure subroutine calculate_toroidal_shift2(j, g_value, x_contour,z_contour, jacob, m, &
     & tor_shift, tor_shift_left_bdry_minus_one)
  !cf. calculate_toroidal_shift, the differences are (1) g_value instead of psival is provided as input;
  !(2) compute a series of tor_shift on a magnetic surface (rather than a single value). The range of theta is [-pi:pi]
  use constants, only: p_, one_half, two
  use magnetic_coordinates, only: dtheta, sign_of_jacobian,sign_of_GSpsi_prime,itheta0, GSpsi_array, GSpsi_prime, qhat
  use domain_decomposition, only: dtheta2, multi_eq_cells
  use magnetic_field, only: psi_z_func, psi_r_func
  implicit none
  integer,intent(in) :: j, m
  real(p_),intent(in) :: g_value
  real(p_),intent(in) :: x_contour(m), z_contour(m), jacob(m)
  real(p_),intent(out) :: tor_shift(m), tor_shift_left_bdry_minus_one
  real(p_) :: gx0, dl, sum, x_mid, z_mid, jac_mid, qhat_mid
  integer :: i

!!$  tor_shift(itheta0) = 0._p_
!!$  do i=itheta0+1, m !for theta in (0:pi]
!!$     x_mid=(x_contour(i-1)+x_contour(i))/two
!!$     z_mid=(z_contour(i-1)+z_contour(i))/two
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
!!$     tor_shift(i) = tor_shift(i-1) + g_value/(x_mid*gx0)*dl
!!$  enddo
!!$
!!$  do i=itheta0-1, 1, -1 !for theta in (0:-pi]
!!$     x_mid=(x_contour(i+1)+x_contour(i))/two
!!$     z_mid=(z_contour(i+1)+z_contour(i))/two
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2) 
!!$     dl=sqrt((x_contour(i)-x_contour(i+1))**2+(z_contour(i)-z_contour(i+1))**2)
!!$     tor_shift(i) = tor_shift(i+1) - g_value/(x_mid*gx0)*dl
!!$  enddo
!!$  tor_shift = tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign

  !-------------another method------------:
  tor_shift(itheta0) = 0._p_
  do i = itheta0 + 1, m !for theta in (0:pi]
     x_mid = (x_contour(i-1) + x_contour(i))/two
     jac_mid = (jacob(i-1) + jacob(i))/two
     qhat_mid = (qhat(i-1,j) + qhat(i,j))/two
     !tor_shift(i) = tor_shift(i-1) - g_value/(x_mid**2)*jac_mid/GSpsi_prime*dtheta
     tor_shift(i) = tor_shift(i-1) + qhat_mid*dtheta
  enddo

  do i = itheta0 - 1, 1, -1 !for theta in (0:-pi]
     x_mid = (x_contour(i+1)+x_contour(i))/two
     jac_mid = (jacob(i+1)+jacob(i))/two
     qhat_mid = (qhat(i+1,j) + qhat(i,j))/two
     !tor_shift(i) = tor_shift(i+1) + g_value/(x_mid**2)*jac_mid/GSpsi_prime*dtheta
     tor_shift(i) = tor_shift(i+1) - qhat_mid*dtheta
  enddo
!---------------------------------------
  !another way of calculating tor_shift
!!$ tor_shift(1)=0._p_
!!$  do i=1,m-1 
!!$     x_mid=(x_contour(i+1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i+1)+z_contour(i))*one_half
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i+1)-x_contour(i))**2+(z_contour(i+1)-z_contour(i))**2)
!!$     tor_shift(i+1)=tor_shift(i)+g_value/(x_mid*gx0)*dl
!!$  enddo
!!$  tor_shift=tor_shift-tor_shift(itheta0)

!!$  tor_shift(1)=0._p_
!!$  do i=2,m
!!$     x_mid=(x_contour(i-1)+x_contour(i))*one_half
!!$     z_mid=(z_contour(i-1)+z_contour(i))*one_half
!!$     gx0=sqrt((psi_r_func(x_mid,z_mid))**2+(psi_z_func(x_mid,z_mid))**2)
!!$     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
!!$     tor_shift(i)=tor_shift(i-1)+g_value/(x_mid*gx0)*dl
!!$  enddo
  !tor_shift = tor_shift*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign

  tor_shift_left_bdry_minus_one = tor_shift(1) - (tor_shift(m) - tor_shift(m-multi_eq_cells))

end subroutine calculate_toroidal_shift2


subroutine calculate_toroidal_shift_at_theta_cut(psival,x_contour,z_contour,np,tor_shift_a,tor_shift_b)
  use constants,only: p_, one_half
  use magnetic_coordinates,only: sign_of_jacobian, sign_of_GSpsi_prime, itheta0
  use domain_decomposition,only:myid
  use magnetic_field, only : psi_z_func,psi_r_func, g_func
  implicit none
  real(p_),intent(in):: psival
  integer,intent(in):: np
  real(p_),intent(in):: x_contour(np),z_contour(np)
  real(p_),intent(out):: tor_shift_a, tor_shift_b
  real(p_):: x_mid,z_mid,gx0,dl, g_value
  integer:: i

  g_value=g_func(psival)

  tor_shift_b=0._p_
  do i=itheta0+1,np ! for location near the cut above the midplane
     x_mid=(x_contour(i)+x_contour(i-1))*one_half
     z_mid=(z_contour(i)+z_contour(i-1))*one_half
     gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
     tor_shift_b=tor_shift_b+g_value/(x_mid*gx0)*dl
  enddo

  tor_shift_a=0._p_
  do i=itheta0,2,-1 ! for location near the cut below the midplane
     x_mid=(x_contour(i)+x_contour(i-1))*one_half
     z_mid=(z_contour(i)+z_contour(i-1))*one_half
     gx0=sqrt((psi_z_func(x_mid,z_mid))**2+(psi_r_func(x_mid,z_mid))**2)
     dl=-sqrt((x_contour(i)-x_contour(i-1))**2+(z_contour(i)-z_contour(i-1))**2)
     tor_shift_a=tor_shift_a+g_value/(x_mid*gx0)*dl
  enddo

  tor_shift_a=tor_shift_a*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
  tor_shift_b=tor_shift_b*(-sign_of_jacobian)/sign_of_GSpsi_prime !include the correct sign
end subroutine calculate_toroidal_shift_at_theta_cut
end module calculate_toroidal_shift_module


subroutine calculate_metric()
  use constants, only: p_, twopi
  use magnetic_coordinates,only: mpol,nrad,r_mc,z_mc,dtheta,dradcor, &
       & tor_shift_mc, tor_shift_left_bdry_minus_one, & !output
       & qhat, jacobian, & !as output
       & pfn_inner,pfn_bdry,GSpsi_array, &
       & zgrid,xgrid, &
       & sign_of_jacobian,sign_of_GSpsi_prime, GSpsi_prime, &
       & grad_psi, grad_alpha,grad_theta,&
       & grad_psi_dot_grad_alpha,grad_psi_dot_grad_theta,grad_alpha_dot_grad_theta !output
  use radial_module, only: psi_lcfs,psi_axis, npsi, sign_bphi, qpsi
  use radial_module, only: q_with_sign, qrad !as output
  use magnetic_field, only : g_func, qfunc
  use calculate_toroidal_shift_module, only : calculate_toroidal_shift2
  use domain_decomposition, only: myid
  use control_parameters, only: diagnosis
  implicit none
  integer :: i, j, ierr
  real(p_) :: g

  call calculate_gradients_of_psi_and_theta() !computing gradients at mc gridpoints

  if(myid==0) write(*,*) 'max_jacobian=', maxval(jacobian), 'min_jacobian=', minval(jacobian)

  sign_of_jacobian=sign(1._p_,jacobian(mpol/3,nrad/3))
  sign_of_GSpsi_prime=sign(1._p_,GSpsi_prime) !radial coordinator is assumed to be pfn

  allocate(q_with_sign(npsi))
  q_with_sign = abs(qpsi)*(-sign_bphi)*sign_of_jacobian*sign_of_GSpsi_prime !include the correct sign
  ! after this call, function "qfunc" is ready to be used
  allocate(qrad(nrad))
  do i = 1, nrad
     qrad(i) = qfunc(xgrid(i))
  enddo
  if(myid==0) write(*,*) 'q_inner=', qrad(1), 'q_bdry=', qrad(nrad)

  allocate(tor_shift_mc(mpol, nrad))
  allocate(tor_shift_left_bdry_minus_one(nrad)) 

  ! local safety factor
  allocate(qhat(mpol,nrad))
  do j = 1, nrad
     g = g_func(GSpsi_array(j))
     do i = 1, mpol
        qhat(i,j) = -g/r_mc(i,j)**2*jacobian(i,j)/GSpsi_prime
     enddo
  enddo

  do j=1,nrad
     g=g_func(GSpsi_array(j))
     call calculate_toroidal_shift2(j, g, r_mc(:,j), z_mc(:,j), jacobian(:,j), mpol, &
          & tor_shift_mc(:,j), tor_shift_left_bdry_minus_one(j))

  enddo

!!$ if(myid.eq.0) call draw_grids_on_theta_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc)
!!$ if(myid.eq.0) call draw_alpha_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc)
!!$ if(myid.eq.0) call draw_alpha_contours_on_a_magnetic_surface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !i.e., magnetic field lines

  call calculate_gradients_of_generalized_toroidal_angle()

  if(myid==0) then
     call diagnostic1()
     call diagnostic2()
  endif

contains
  subroutine diagnostic1()
    integer:: i,u
    write(*,*) 'maximum of tor_shift_mc=', maxval(tor_shift_mc)
    write(*,*) 'minimum of tor_shift_mc=', minval(tor_shift_mc)
    open(newunit=u,file='tor_shift_mc.txt')
    do j=1,nrad
       do i=1,mpol
          write(u,*) r_mc(i,j), z_mc(i,j), tor_shift_mc(i,j), qhat(i,j)
       enddo
       write(u,*)
       write(u,*)
    enddo
    close(u)
  end subroutine diagnostic1

  subroutine diagnostic2()
    use magnetic_coordinates, only : pfn, dl_mc
    use func_in_mc, only: minor_r_radcor
    integer:: i,u

    real(p_) :: s, q2
    open(newunit=u,file="mc_derivatives_mc3.txt")
    do i=1,mpol
       do j=1,nrad
          write(u,'(2i8.4,9(1pe14.5))')  i,j,grad_alpha(i,j),grad_psi_dot_grad_alpha(i,j),&
               & grad_psi(i,j),jacobian(i,j),grad_alpha_dot_grad_theta(i,j),&
               & tor_shift_mc(i,j),grad_psi_dot_grad_theta(i,j),zgrid(i),xgrid(j)
       enddo
       write(u,*)
    enddo
    close(u)
!!$     do j=1,nrad
!!$        write(*,*) j,xgrid(j),minor_r_radcor(xgrid(j)),r_mag_surf0(1,j)-r_axis,&
!!$             & minor_r_prime(xgrid(j))
!!$     enddo
    open(newunit=u,file="q2.txt")
    do j=1,nrad
       s=0
       do i=1,mpol-1
          s = s + g_func(GSpsi_array(j))/(grad_psi(i,j)*(psi_lcfs-psi_axis)*r_mc(i,j))*dl_mc(i,j)
       enddo
       q2 = s/twopi
       write(u,'(20(ES16.4E2))') pfn(j), minor_r_radcor(pfn(j)), q2, qrad(j)
    enddo
    close(u)
  end subroutine diagnostic2

end subroutine calculate_metric


subroutine calculate_gradients_of_psi_and_theta()
  use constants, only: p_, zero, one, two
  use magnetic_coordinates, only: m=>mpol, n=>nrad, r=>r_mc, z=>z_mc, dtheta, dradcor,&
   & jacobian, abs_jacobian, jacobian_av, & !as output
   & rth, zth, rpsi, zpsi, grad_psi, grad_theta, grad_psi_dot_grad_theta, & !as output
   & grad_psi_r, grad_psi_z, grad_theta_r, grad_theta_z !as output
  use domain_decomposition, only: myid
  !use magnetic_field, only: minor_r_prime, minor_r_radcor
  implicit none
  integer :: i,j,u
  real(p_) ::  minor_r_prime_val
  
  allocate(jacobian(m,n))
  allocate(abs_jacobian(m,n))
  allocate(jacobian_av(n))
  allocate(grad_psi(m,n))
  allocate(grad_theta(m,n))
  allocate(grad_psi_dot_grad_theta(m,n))
  allocate(grad_theta_r(m,n))
  allocate(grad_theta_z(m,n))
  allocate(grad_psi_r(m,n))
  allocate(grad_psi_z(m,n))
  allocate(rpsi(m,n))
  allocate(zpsi(m,n))
  allocate(rth(m,n))
  allocate(zth(m,n))

  !call partial_derivative_in_mc(m,n,r,z,dtheta,dradcor,rpsi,rth,zpsi,zth,jacobian)
  call partial_derivative_in_mc2(m,n,r,z,dtheta,dradcor,rpsi,rth,zpsi,zth,jacobian) 

  abs_jacobian = abs(jacobian)
  
  grad_psi=r/abs_jacobian*sqrt(zth**2+rth**2)
  grad_theta=r/abs_jacobian*sqrt(zpsi**2+rpsi**2)
  grad_psi_dot_grad_theta=-(r/jacobian)**2*(zth*zpsi+rth*rpsi)

  grad_psi_r = -r/jacobian*zth
  grad_psi_z = r/jacobian*rth
  grad_theta_r = r/jacobian*zpsi
  grad_theta_z = -r/jacobian*rpsi

  do j = 1,n
     jacobian_av(j) = sum(abs_jacobian(1:m-1,j))/(m-1) !uniform theta grid is assumed
  enddo
!  if(myid.eq.0) call plot_psi_r_z_theta_r_z_mc()

contains
!!$  subroutine plot_psi_r_z_theta_r_z_mc()
!!$    use magnetic_coordinates,only: zgrid,xgrid
!!$    integer:: u
!!$    real(p_):: minor_r_val,minor_r_prime_val
!!$
!!$    open(newunit=u, file='mc_derivatives_mc1.txt')
!!$    do i=1,m
!!$       do j=1,n
!!$          write(u,*) r(i,j),z(i,j), jacobian(i,j)
!!$
!!$          !minor_r_prime_val=minor_r_prime(xgrid(j))
!!$          !write(u,*) r(i,j),z(i,j),jacobian(i,j)/minor_r_prime_val !for equal arc-length poloidal angle, jacobian(i,j)/minor_r_prime_val is equal -R*minor_r.
!!$          !write(u,*) r(i,j),z(i,j),grad_psi_r(i,j)*minor_r_prime_val, grad_psi_z(i,j)*minor_r_prime_val,&
!!$          !write(u,*) r(i,j),z(i,j),grad_theta_r(i,j),grad_theta_z(i,j),grad_theta(i,j),jacobian(i,j)
!!$       enddo
!!$       write(u,*)
!!$    enddo
!!$    close(u)
!!$
!!$    open(newunit=u, file='psi_r_z_theta_r_z_mc_analytic.txt') !for concentric-circular magnetic field (equal-arc-length theta)
!!$    do i=1,m
!!$       do j=1,n
!!$          !minor_r_val=minor_r_radcor(xgrid(j))
!!$          write(u,*) r(i,j),z(i,j),cos(zgrid(i)),sin(zgrid(i)),&
!!$               &-sin(zgrid(i))/minor_r_val,cos(zgrid(i))/minor_r_val,-R(i,j)*minor_r_val
!!$       enddo
!!$       write(u,*)
!!$    enddo
!!$    close(u)
!!$  end subroutine plot_psi_r_z_theta_r_z_mc
end subroutine calculate_gradients_of_psi_and_theta


subroutine calculate_gradients_of_generalized_toroidal_angle() 
  use constants,only: p_, zero,one,two
  use control_parameters, only : diagnosis
  use magnetic_coordinates,only: m=>mpol, n=>nrad, r=>r_mc, z=>z_mc, &
  & tor_shift_mc,jacobian,dtheta,dradcor, &
  & rpsi,zpsi,rth,zth, &
  & grad_alpha_r, grad_alpha_z, grad_alpha_phi, & !as output 
  & grad_alpha, grad_psi_dot_grad_alpha, grad_alpha_dot_grad_theta, & !as output
  & grad_alpha_r_left_bdry_minus_one, grad_alpha_z_left_bdry_minus_one !as output 
  use domain_decomposition,only: myid, dtheta2
  implicit none
  real(p_) :: tor_shift_psi(m,n), tor_shift_th(m,n)
  real(p_) :: tor_shift_psi_left_bdry_minus_one(n)
  integer :: i,j

  allocate(grad_alpha_r(m,n))
  allocate(grad_alpha_z(m,n))
  allocate(grad_alpha_phi(m,n))
  allocate(grad_alpha(m,n))
  allocate(grad_psi_dot_grad_alpha(m,n))
  allocate(grad_alpha_dot_grad_theta(m,n))
  allocate(grad_alpha_r_left_bdry_minus_one(n))
  allocate(grad_alpha_z_left_bdry_minus_one(n))

  call partial_derivative_of_tor_shift_in_mc2(m,n,tor_shift_mc,dtheta,dradcor,tor_shift_psi,tor_shift_th, &
       & tor_shift_psi_left_bdry_minus_one) 
  grad_alpha_r = tor_shift_psi*r/jacobian*zth - tor_shift_th*r/jacobian*zpsi
  grad_alpha_z = tor_shift_th*r/jacobian*rpsi - tor_shift_psi*r/jacobian*rth
  grad_alpha_phi = one/r
  grad_alpha=sqrt(grad_alpha_phi**2+grad_alpha_r**2+grad_alpha_z**2)
  grad_psi_dot_grad_alpha=-r/jacobian*zth*grad_alpha_r+r/jacobian*rth*grad_alpha_z
  grad_alpha_dot_grad_theta=grad_alpha_r*r/jacobian*zpsi-grad_alpha_z*r/jacobian*rpsi

  !i=m-1 !wrong
  i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
  do j=1,n
     grad_alpha_r_left_bdry_minus_one(j)=tor_shift_psi_left_bdry_minus_one(j)*r(i,j)/jacobian(i,j)*zth(i,j)&
          & -tor_shift_th(i,j)*r(i,j)/jacobian(i,j)*zpsi(i,j)
     grad_alpha_z_left_bdry_minus_one(j)=tor_shift_th(i,j)*r(i,j)/jacobian(i,j)*rpsi(i,j)&
          & -tor_shift_psi_left_bdry_minus_one(j)*r(i,j)/jacobian(i,j)*rth(i,j)
  enddo
  !if(myid==0 .and. (diagnosis .eqv. .true.)) call plot_alpha_r_z_mc()
  if(myid==0 .and. (diagnosis .eqv. .true.)) call verification2()
  
contains
  ! subroutine plot_alpha_r_z_mc()
  !   use magnetic_coordinates,only: xgrid,zgrid
  !   use magnetic_field,only: minor_r_radcor, qfunc

  !   integer :: u
  !   real(p_) :: theta,minor_r,major_r0,major_r,q,local_q,dq_dr,factor
  !   real(p_) :: ddelta_minor_r,alpha_r,alpha_z
  !   real(p_) :: tmp,tmp2,dpsi_dr
  !   real(p_), parameter :: g0=1.32*1.91,psi_axis=0.,  psi_lcfs=  0.146018378_p_
  !   open(newunit=u,file='mc_derivatives_mc2.txt')
  !   do i=1,m
  !      do j=1,n
  !         !write(u,*) r(i,j),z(i,j),grad_alpha_r(i,j),grad_alpha_z(i,j), grad_alpha(i,j)
  !         write(u,*) r(i,j),z(i,j),tor_shift_th(i,j),tor_shift_psi(i,j)
  !         !write(u,*) r(i,j),z(i,j),tor_shift_mc(i,j)
  !      enddo
  !      write(u,*)
  !      write(u,*)
  !   enddo

  !   i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
  !   do j=1,n
  !      write(u,*) r(i,j),z(i,j), tor_shift_th(i,j),tor_shift_psi_left_bdry_minus_one(j)
  !   enddo
  !   close(u)
  !   !write(*,*) 'maxval(grad_alpha)=',maxval(grad_alpha),'minval(grad_alpha)=',minval(grad_alpha)
  !   open(newunit=u,file='alpha_r_z_mc_analytic.txt') !for concentric-circular magnetic field (equal-arc-length theta and geometric minor radius as the flux surface label)
  !   do i=1,m
  !      theta=zgrid(i)
  !      do j=1,n
  !         major_r0=1.32_p_!for DIII-D cyclone base case parameter
  !         minor_r=minor_r_radcor(xgrid(j))
  !         major_r=major_r0+minor_r*cos(theta)
  !         q=qfunc(xgrid(j))
  !         local_q=q*sqrt(major_r0**2-minor_r**2)/major_r
  !         dq_dr=0.78*1.4/0.24_p_ !for DIII-D cyclone base case parameter, a linear q profile is assumed, the same as what Ben did
  !         factor=(major_r0-minor_r)/sqrt(major_r0**2-minor_r**2)*tan(theta/two)
  !         ddelta_minor_r=two*dq_dr*atan(factor)+two*q/(1+factor**2)*tan(theta/2._p_)*&
  !              & (-major_r0)/((major_r0+minor_r)*sqrt(major_r0**2-minor_r**2))
  !         alpha_r=-ddelta_minor_r*cos(theta)+local_q*sin(theta)/minor_r
  !         alpha_z=-local_q*cos(theta)/minor_r-ddelta_minor_r*sin(theta)
  !         !write(u,*) r(i,j),z(i,j), alpha_r,alpha_z
  !         tmp=(major_r0-minor_r)/sqrt(major_r0**2-minor_r**2)
  !         tmp2=tan(theta/2)
  !         !dpsi_dr=one/minor_r_prime(xgrid(j)) 
  !         dpsi_dr=g0*minor_r/(q*sqrt(major_r0**2-minor_r**2)*(psi_lcfs-psi_axis)) !another way to calculate dpsi_dr
  !         !write(u,*) r(i,j),z(i,j), local_q, 2*q*tmp*(tan(theta/2)**2/2+0.5)/(tmp**2*tan(theta/2)**2+1)
  !         write(u,*) r(i,j),z(i,j), local_q, ddelta_minor_r/dpsi_dr !local_q is equal to ddelta_dtheta
  !         !write(u,*) r(i,j),z(i,j), (2*atan(tmp2*(major_r0 - minor_r)/sqrt(major_r0**2 - minor_r**2))*dq_dr&
  !         !     &  + 2*(tmp2*minor_r*(major_r0 - minor_r)/(major_r0**2 - minor_r**2)**(3./2) &
  !         !     & - tmp2/sqrt(major_r0**2 - minor_r**2))*q/(tmp2**2*(major_r0 - minor_r)**2&
  !         !     & /(major_r0**2 - minor_r**2) + 1))/dpsi_dr !ddelta/dradcor=ddelta/dr*(dr/dradcor), where the formula for computing ddelta/dr is obtained by using sympy

  !         !write(u,*) r(i,j),z(i,j), 2*q*atan(tmp*tan(theta/2))
  !      enddo
  !      write(u,*)
  !      write(u,*)
  !   enddo
  !   close(u)
  ! end subroutine plot_alpha_r_z_mc

  subroutine verification2() !verify grad_psi_cross_grad_alpha=B0/GSpsi_prime
    use magnetic_coordinates,only: xgrid, grad_psi_r,grad_psi_z, GSpsi_prime
    use math,only:  cross_product_in_cartesian
    use magnetic_field,only : br,bz,bphi
    real(p_):: ax,ay,az, cx,cy,cz,dx,dy,dz
    real(p_):: brval,bzval,bphival
    integer:: u,i,j
    open(newunit=u,file='cross_product_comparision.txt')

!    i=1
!    i=m

 i=m-NINT(dtheta2/dtheta) !dtheta2 is the grid spacing for the perturbations and dtheta is the equilibrium grid spacing.
    do j=1,n
       ax=grad_psi_r(i,j)
       ay=0._p_
       az=grad_psi_z(i,j)

!!$       cx=grad_alpha_r(i,j)
!!$       cy=grad_alpha_phi(i,j)
!!$       cz=grad_alpha_z(i,j)

       cx=grad_alpha_r_left_bdry_minus_one(j)
       cy=grad_alpha_phi(i,j)
       cz=grad_alpha_z_left_bdry_minus_one(j)

       call cross_product_in_cartesian(ax,ay,az,cx,cy,cz,dx,dy,dz) !grad_psi_cross_grad_alpha

        brval=br(r(i,j),z(i,j))
        bzval=bz(r(i,j),z(i,j))
        bphival=bphi(r(i,j),z(i,j))
       write(u,*) dx,dy,dz, brval/GSpsi_prime,bphival/GSpsi_prime,bzval/GSpsi_prime
    enddo

close(u)
  end subroutine verification2
end subroutine calculate_gradients_of_generalized_toroidal_angle


subroutine partial_derivative_in_mc(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)
  !calculate the partial derivative of R and Z with respect to the magnetic cooordinates (psi,theta)
  !jacob is also calculated in this subroutine
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: r(m,n),z(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: rpsi(m,n),rth(m,n),zpsi(m,n),zth(m,n),jacob(m,n)
  real(p_):: tmp0
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        rpsi(i,j)=(r(i,j+1)-r(i,j-1))/(two*dpsi)
        zpsi(i,j)=(z(i,j+1)-z(i,j-1))/(two*dpsi)
     enddo

     !use linear interpolation to get the value  j=n
     tmp0=(r(i,n)-r(i,n-1))/dpsi
     rpsi(i,n)=two*tmp0-rpsi(i,n-1)
     tmp0=(z(i,n)-z(i,n-1))/dpsi
     zpsi(i,n)=two*tmp0-zpsi(i,n-1)

     !use linear interpolation to get the value j=1
     tmp0=(r(i,2)-r(i,1))/dpsi
     rpsi(i,1)=two*tmp0-rpsi(i,2)

     tmp0=(z(i,2)-z(i,1))/dpsi
     zpsi(i,1)=two*tmp0-zpsi(i,2)

  enddo

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        rth(i,j)= (r(i+1,j)-r(i-1,j))/(two*dtheta)
        zth(i,j)=(z(i+1,j)-z(i-1,j))/(two*dtheta)
     enddo

     !use peroidic property of r and z to calculate the partial derivative for boundary points at theta=0 and 2pi
     rth(1,j)=(r(2,j)-r(m-1,j))/(two*dtheta)
     zth(1,j)=(z(2,j)-z(m-1,j))/(two*dtheta)
     rth(m,j)=rth(1,j)
     zth(m,j)=zth(1,j)
  enddo

  !calculate the Jacobian:
  do i=1,m
!          do j=2,n !the jacobian at the magnetic axis is zero
     do j=1,n
        jacob(i,j)=r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,fai)
     enddo
     !jacob(i,n)=two*jacob(i,n-1)-jacob(i,n-2)
     !use linear interpolation to get the value of jacobian at the magnetic surface near the magnetic axis
  enddo

  !write(*,*) jacob(1:m,5)
end subroutine partial_derivative_in_mc

subroutine partial_derivative_in_mc2(m,n,r,z,dtheta,dpsi,rpsi,rth,zpsi,zth,jacob)
  !calculate the partial derivative of R and Z with respect to the magnetic cooordinates (psi,theta)
  !jacob is also calculated in this subroutine
  use constants,only: p_
  use magnetic_coordinates, only : radcor=>xgrid, theta=>zgrid
  use splines, only: spline3ders
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: r(m,n), z(m,n)
  real(p_), intent(in) :: dtheta, dpsi
  real(p_), intent(out) :: rpsi(m,n), rth(m,n), zpsi(m,n), zth(m,n), jacob(m,n)
  integer :: i, j

  do i = 1, m
     call spline3ders(radcor, r(i,:), radcor, dynew=rpsi(i,:))
     call spline3ders(radcor, z(i,:), radcor, dynew=zpsi(i,:))
  enddo

  do j=1,n
     call spline3ders(theta, r(:,j), theta, dynew=rth(:,j))
     call spline3ders(theta, z(:,j), theta, dynew=zth(:,j))
  enddo

  do i=1,m
     do j=1,n
        jacob(i,j) = r(i,j)*(rth(i,j)*zpsi(i,j)-rpsi(i,j)*zth(i,j))   !Jacobain of coordinate system (psi,theta,phi)
     enddo
     !jacob(i,1) = 2*jacob(i,2) - jacob(i,3)
  enddo
  
end subroutine partial_derivative_in_mc2


subroutine partial_derivative_of_tor_shift_in_mc(m,n,tor_shift,dtheta,dpsi,tor_shift_psi,tor_shift_th,&
     & tor_shift_psi_left_bdry_minus_one)
  !calculate the partial derivative with respect to the magnetic cooordinates (psi,theta)
  use constants,only:p_
  use constants,only:zero,one,two,twopi,one_half
  use magnetic_coordinates,only: tor_shift_left_bdry_minus_one
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: tor_shift(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: tor_shift_psi(m,n),tor_shift_th(m,n)
  real(p_),intent(out):: tor_shift_psi_left_bdry_minus_one(n)
  real(p_):: tmp0,twopi_q
  integer:: i,j

  do i=1,m  
     do j=2,n-1 !use center difference scheme for inner points
        tor_shift_psi(i,j)=(tor_shift(i,j+1)-tor_shift(i,j-1))/(two*dpsi)
     enddo
     !use linear interpolation to get the value at j=n
     tmp0=(tor_shift(i,n)-tor_shift(i,n-1))/dpsi
     tor_shift_psi(i,n)=two*tmp0-tor_shift_psi(i,n-1)
     !use linear interpolation to get the value j=1
     tmp0=(tor_shift(i,2)-tor_shift(i,1))/dpsi
     tor_shift_psi(i,1)=two*tmp0-tor_shift_psi(i,2)
  enddo

  do j=2,n-1
     tor_shift_psi_left_bdry_minus_one(j)=&
          & (tor_shift_left_bdry_minus_one(j+1)-tor_shift_left_bdry_minus_one(j-1))/(two*dpsi)
  enddo
  !use linear interpolation to get the value at j=n
  tmp0=(tor_shift_left_bdry_minus_one(n)-tor_shift_left_bdry_minus_one(n-1))/dpsi
  tor_shift_psi_left_bdry_minus_one(n)=two*tmp0-tor_shift_psi_left_bdry_minus_one(n-1)
  !use linear interpolation to get the value j=1
  tmp0=(tor_shift_left_bdry_minus_one(2)-tor_shift_left_bdry_minus_one(1))/dpsi
  tor_shift_psi_left_bdry_minus_one(1)=two*tmp0-tor_shift_psi_left_bdry_minus_one(2)

  do j=1,n
     do i=2,m-1 !use center difference scheme for inner points
        tor_shift_th(i,j)= (tor_shift(i+1,j)-tor_shift(i-1,j))/(two*dtheta)
     enddo
     !for boundary points at theta cut
!!$     twopi_q=tor_shift(m,j)
!!$     tor_shift_th(1,j)=(tor_shift(2,j)+twopi_q-tor_shift(m-1,j))/(two*dtheta)
!!$     tor_shift_th(m,j)=tor_shift_th(1,j)
     !tor_shift_th(1,j)=(tor_shift(2,j)-tor_shift_left_bdry_minus_one(j))/(two*dtheta), wrong! since left_bdry_minus_one is defined on the grids for perturbations, which are different from equilibrium grids.
     tor_shift_th(1,j)=two*tor_shift_th(2,j)-tor_shift_th(3,j)   !use linear interpolation to get the value at i=1 and i=m
     tor_shift_th(m,j)=tor_shift_th(1,j) !tor_shift_th is periodic

  enddo
end subroutine partial_derivative_of_tor_shift_in_mc


subroutine partial_derivative_of_tor_shift_in_mc2(m,n,tor_shift,dtheta,dpsi,tor_shift_psi,tor_shift_th,&
     & tor_shift_psi_left_bdry_minus_one)
  !calculate the partial derivative with respect to the magnetic cooordinates (psi,theta)
  use constants,only:p_,zero,one,two,twopi,one_half
  use magnetic_coordinates,only: tor_shift_left_bdry_minus_one, xgrid, zgrid
  use splines, only: spline3ders
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: tor_shift(m,n)
  real(p_),intent(in):: dtheta,dpsi
  real(p_),intent(out):: tor_shift_psi(m,n),tor_shift_th(m,n)
  real(p_),intent(out):: tor_shift_psi_left_bdry_minus_one(n)
  integer:: i,j
  real(p_) :: tmp(m,n), tmp0(n)

  do i = 1, m
     call spline3ders(xgrid, tor_shift(i,:), xgrid, tmp(i,:), tor_shift_psi(i,:), tmp(i,:))
  enddo

  do j = 1, n
     call spline3ders(zgrid, tor_shift(:,j), zgrid, tmp(:,j), tor_shift_th(:,j), tmp(:,j))
  enddo

  call spline3ders(xgrid, tor_shift_left_bdry_minus_one(:), xgrid, &
       & tmp0(:), tor_shift_psi_left_bdry_minus_one(:), tmp0(:))


end subroutine partial_derivative_of_tor_shift_in_mc2



function jacobian_func(theta,radcor) result (z) 
  use constants,only: p_
  use magnetic_coordinates,only: mpol,nrad,zgrid,xgrid,jacobian !as input
  use interpolate_module,only: linear_2d_interpolate

  implicit none
  real(p_)::radcor,theta,z
  call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,theta,radcor,z)  !uniform 1darray is assumed
end function jacobian_func
subroutine mapping_table_for_cylindrical_to_magnetic_coordinates()
  !Prepare numerical tables radcor(R,Z), theta(R,Z), and tor_shift(R,Z).
  !Boris pusher works in cylindrical coordinates, so we need these pre-mapping data to interpolate the Boris orbit into magnetic coordinates.
  !Later I decided that tor_shift(R,Z) is not used in the rest of the program. The value of tor_shift and the various gradients of psi, theta, and alpha will be computed by interpolating the numerical table in magnetic coordinates (numerical talbes are created in "calculate_metric" routine)
  use constants, only: p_, pi
  use boundary, only: np_lcfs !,x_lcfs,z_lcfs
  use mapping_module, only: nx_mapping, nz_mapping
  use mapping_module, only: r_cyl, z_cyl, dr, dz, i0, j0 !as output
  use mapping_module, only: radcor, theta_a, theta_b, tor_shift_a, tor_shift_b !as output
  use domain_decomposition, only: myid
  use math, only: pnpoly
  use magnetic_field, only: psi_func
  implicit none
  logical :: within_region(nx_mapping,nz_mapping)
  integer :: i, j, inout1(nx_mapping,nz_mapping), inout2(nx_mapping,nz_mapping)
  real(p_) :: theta(nx_mapping,nz_mapping), tor_shift(nx_mapping,nz_mapping)
  real(p_) :: r_inner_surf(np_lcfs), z_inner_surf(np_lcfs), r_outer_surf(np_lcfs), z_outer_surf(np_lcfs)

  call choose_boundary_magnetic_surfaces_for_the_mapping(r_inner_surf,z_inner_surf,r_outer_surf,z_outer_surf)
  call create_cylindrical_grids(r_outer_surf,z_outer_surf,np_lcfs,nx_mapping,nz_mapping,r_cyl,z_cyl,dr,dz,i0,j0)

  do i=1,nx_mapping
     do j=1,nz_mapping
!!$        call PNPOLY(r_cyl(i),z_cyl(j),r_mc(:,nrad),z_mc(:,nrad),mpol,INOUT1(i,j))
!!$        call PNPOLY(r_cyl(i),z_cyl(j),r_mc(:,2),z_mc(:,2),mpol,INOUT2(i,j))
        call PNPOLY(r_cyl(i),z_cyl(j),r_outer_surf,z_outer_surf,np_lcfs,INOUT1(i,j))
        call PNPOLY(r_cyl(i),z_cyl(j),r_inner_surf,z_inner_surf,np_lcfs,INOUT2(i,j))
     enddo
  enddo

  within_region = .true. !whether a point is within the specifed region.
  do i = 1, nx_mapping 
     do j = 1, nz_mapping
        if((inout1(i,j).eq. -1) .or. (inout2(i,j).eq.1)) within_region(i,j) = .false.
     enddo
  enddo

 if(myid.eq.0) call diagnostic1()
  ! call arrange_lcfs(x_lcfs,z_lcfs,np_lcfs,r_axis,z_axis)

  radcor = 0._p_ !initialize
  theta = 0._p_
  tor_shift = 0._p_

  do i = 1, nx_mapping
     do j = 1, nz_mapping
        if(i<i0 .and. j.eq.j0) then !theta cut, special treatment is needed here
           !done at the end of this subroutine
        else
           if(within_region(i,j).eqv..true.)  call mapping(r_cyl(i),z_cyl(j),radcor(i,j),theta(i,j),tor_shift(i,j)) 
        endif
     enddo
  enddo

  theta_a = theta
  theta_b = theta
  tor_shift_a = tor_shift
  tor_shift_b = tor_shift

  j = j0 !special treatment is needed at theta cut
  do i = 1, i0-1
     theta_a(i,j) = -pi
     theta_b(i,j) = pi
     call mapping2(r_cyl(i), z_cyl(j), radcor(i,j), tor_shift_a(i,j), tor_shift_b(i,j))
  enddo

!!$  tor_shift_b=tor_shift
!!$  do i=i0,nx_mapping
!!$     j=j0 !re-calculate tor_shift_angle at low-field-side midplane (2pi*q, instead of zero) and store them in tor_shift_b array
!!$     if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$  enddo
!!$  tor_shift_b=tor_shift
!!$  do i=1,nx_mapping
!!$     do j=1,nz_mapping
!!$        if(i<i0 .and. j.eq.j0) then !the value of theta and tor_shift at the high-field-side midplane
!!$           if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$           if(within_region(i,j).eqv..true.)   call mapping2(r_cyl(i),z_cyl(j),tor_shift_b(i,j))
!!$        endif
!!$     enddo
!!$  enddo
  if(myid.eq.0) call diagnostic2()
contains
  subroutine diagnostic1()
    integer:: u1, u2
    open(newunit=u1,file='mapping_in.txt')
    open(newunit=u2,file='mapping_out.txt')
    do i = 1, nx_mapping
       do j = 1, nz_mapping
          if(within_region(i,j) .eqv. .false.) then 
             write(u2,*) r_cyl(i),z_cyl(j) !out
          else
             write(u1,*) r_cyl(i),z_cyl(j) !in
          endif
       enddo
    enddo
    close(u1)
    close(u2)
  end subroutine diagnostic1
  subroutine diagnostic2()
    integer:: u
    open(newunit=u,file='mapping_table.txt')
    do i=1,nx_mapping
       do j=1,nz_mapping 
          !do j=j0,j0
          if(within_region(i,j).eqv..true.) then
             write(u,*) r_cyl(i), z_cyl(j), theta_a(i,j), theta_b(i,j), tor_shift_a(i,j), tor_shift_b(i,j)
             !write(u,*) r_cyl(i),z_cyl(j),psi_func(r_cyl(i),z_cyl(j))
          else
             write(u,*) r_cyl(i),z_cyl(j), 'NaN',' NaN',' NaN',' NaN'
          endif
       enddo
       write(u,*) 
    enddo
    close(u)
    !  write(*,*) 'maximum of tor_shift=',maxval(tor_shift_a),'minimum of tor_shift=',minval(tor_shift_a)
    !  write(*,*) 'maximum of tor_shift=',maxval(tor_shift_b),'minimum of tor_shift=',minval(tor_shift_b)
  end subroutine diagnostic2
end subroutine mapping_table_for_cylindrical_to_magnetic_coordinates


subroutine mapping(r,z,radcor,theta,tor_shift)
  !given (R,Z), this subroutine finds the magnetic surface that passes through the point and calculates its poloidal angle and toroidal shift
  use constants,only : p_, zero,one,two,twopi,pi
  use math, only : arc_length
  use boundary, only: x_lcfs,z_lcfs,np_lcfs
  use radial_module,only:r_axis,z_axis,psi_axis,psi_lcfs
  use magnetic_field,only:radcor_as_func_of_pfn
  use control_parameters,only: poloidal_angle_type
  use domain_decomposition,only:myid
  use calculate_toroidal_shift_module
  use contour_mod,only : find_contour
  use magnetic_field, only : psi_func,psi_gradient_func, b
  implicit none
  real(p_),intent(in):: r,z
  real(p_),intent(out):: radcor,theta,tor_shift

  real(p_)::psival
  real(p_):: x_contour(np_lcfs+1),z_contour(np_lcfs+1)
  real(p_)::dl(np_lcfs), sum

  real(p_),parameter:: xacc=1.0d-6 !tolerance used in bi-section root-finder
  real(p_):: x1,x2,z1,z2

  real(p_):: slope(np_lcfs),slope2(np_lcfs)
  real(p_):: zfunc,xfunc !equation of the straight line (in poloidal plane) that passing throught the magnetic axis point and one point on LCFS
  integer:: i,end_i
  real(p_):: value1,value2,value3, rmid,zmid,normalization,tmpx,tmpz
  real(p_),parameter:: large_number=1d30

  psival=psi_func(r,z)
  radcor=radcor_as_func_of_pfn((psival-psi_axis)/(psi_lcfs-psi_axis))

  call find_contour(psival,x_contour,z_contour)

  call arc_length(x_contour,z_contour,np_lcfs,dl)
!!$     sum=0.
!!$     do i=1,np_lcfs-1
!!$        sum=sum+dl(i)
!!$     enddo
!!$     circumference=sum

  normalization=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
        normalization=normalization+dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,np_lcfs !finish a full poloidal circle integration to get the normalization factor
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        normalization=normalization+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'straight-field-line') then
     do i=2,np_lcfs
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        normalization=normalization+dl(i-1)/(rmid*psi_gradient_func(rmid,zmid))
     enddo

  elseif(poloidal_angle_type .eq. 'Boozer') then
     do i=2,np_lcfs
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        normalization=normalization+dl(i-1)*b(rmid,zmid)**2*rmid/psi_gradient_func(rmid,zmid)
     enddo

  else
     stop 'please choose poloidal angle type between equal-arc/equal-volume/straight-field-line, in mapping'
  endif

  call locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
  ! if(myid.eq.0) call diagnostic1()
  tmpx=x_contour(end_i+1) !backup the original value
  tmpz=z_contour(end_i+1)
  x_contour(end_i+1)=r !replace No. end_i+1 point by the given point (r,z)
  z_contour(end_i+1)=z
  dl(end_i)=sqrt((x_contour(end_i+1)-x_contour(end_i))**2+(z_contour(end_i+1)-z_contour(end_i))**2)

  !calculate poloidal angle 
  theta=0._p_
  if(poloidal_angle_type .eq. 'equal-arc') then
     do i=2,end_i+1
        theta=theta+dl(i-1) !equal-arc-length poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'equal-volume') then
     do i=2,end_i+1
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        theta=theta+dl(i-1)*rmid/psi_gradient_func(rmid,zmid) !equal-volume poloidal angle
     enddo
  elseif(poloidal_angle_type .eq. 'straight-field-line') then
     do i=2,end_i+1
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        theta=theta+dl(i-1)/(rmid*psi_gradient_func(rmid,zmid))
     enddo
  elseif(poloidal_angle_type .eq. 'Boozer') then
     do i=2,end_i+1
        rmid=0.5_p_*(x_contour(i-1)+x_contour(i))
        zmid=0.5_p_*(z_contour(i-1)+z_contour(i))
        theta=theta+dl(i-1)*b(rmid,zmid)**2*rmid/psi_gradient_func(rmid,zmid)
     enddo

  else
     stop 'please choose poloidal angle type between equal-arc and equal-volume'
  endif
  theta=theta*twopi/normalization-pi !normalized to the range [-pi:pi]

  x_contour(end_i+1)=tmpx !restore to the original value
  z_contour(end_i+1)=tmpz
  call calculate_toroidal_shift(psival,x_contour,z_contour,np_lcfs,end_i,r,z,tor_shift) !calculate toroidal shift (\int_0_theta{q_hat dtheta}) which is needed in the definition of the generalized toroidal angle
contains
  subroutine diagnostic1()
    integer:: u
    open(newunit=u,file='tmp_contour')
    do i=1,end_i
       write(u,*) x_contour(i),z_contour(i)
    enddo
    write(u,*) 
    write(u,*) 
  end subroutine diagnostic1
end subroutine mapping

subroutine mapping2(r,z,radcor,tor_shift_a,tor_shift_b)
  !given (R,Z), this subroutine find the magnetic surface that passes throught the point and calculates its toroidal angle shift
  !the same as mapping(), the difference is that this only takes care of the total_tor_shift at the theta cut, see the comments where this subroutine is called
  use constants,only:p_
  use constants,only:zero,one,two,twopi,pi
  use boundary, only: x_lcfs,z_lcfs,np_lcfs
  use radial_module,only:r_axis,z_axis,psi_axis,psi_lcfs
  use calculate_toroidal_shift_module
  use contour_mod,only : find_contour
  use magnetic_field, only : psi_func, radcor_as_func_of_pfn
  implicit none
  real(p_),intent(in):: r,z
  real(p_),intent(out):: radcor,tor_shift_a,tor_shift_b
  real(p_)::psival
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  real(p_)::dl(np_lcfs-1)
  integer:: i,end_i

  psival=psi_func(r,z)
  radcor=radcor_as_func_of_pfn((psival-psi_axis)/(psi_lcfs-psi_axis))
  call find_contour(psival,x_contour,z_contour)

  !  call calculate_toroidal_shift_total(psival,x_contour,z_contour,np_lcfs,dl,tor_shift_a) 
  call calculate_toroidal_shift_at_theta_cut(psival,x_contour,z_contour,np_lcfs,tor_shift_a,tor_shift_b) !calculate toroidal shift which is needed in the definition of the generalized toroidal angle
end subroutine mapping2


subroutine locate_poloidal_index(r,z,x_lcfs,z_lcfs,np_lcfs,end_i) !poloidal index of point (r,z) is between end_i and end_i+1
  use constants,only:p_
  use constants,only:twopi,pi
  use radial_module,only:r_axis,z_axis
  implicit none
  real(p_),intent(in):: r,z
  integer,intent(in):: np_lcfs
  real(p_),intent(in):: x_lcfs(np_lcfs),z_lcfs(np_lcfs)
  integer,intent(out):: end_i !!poloidal index of point (r,z) is between end_i and end_i+1
  real(p_):: xn,zn,xn0,zn0,angle0,angle1,angle2
  integer:: i

  xn0=r-r_axis
  zn0=z-z_axis
  angle0=acos(xn0/sqrt(xn0*xn0+zn0*zn0))
  if(zn0.le.0.0) angle0=-angle0 !theta in [-pi:pi]

  do i=1,np_lcfs-1
     xn=x_lcfs(i)-r_axis
     zn=z_lcfs(i)-z_axis
     angle1=acos(xn/sqrt(xn**2+zn**2)) 
     if(zn<0.0) angle1=-angle1
     if(i.eq.1) angle1=-pi

     xn=x_lcfs(i+1)-r_axis
     zn=z_lcfs(i+1)-z_axis
     angle2=acos(xn/sqrt(xn**2+zn**2))
     if(zn<0.0) angle2=-angle2 !theta in [-pi:pi]
     if(i+1.eq.np_lcfs) angle2=pi

     if((angle0-angle1)*(angle0-angle2).le.0._p_) exit

  enddo


  end_i=i

end subroutine locate_poloidal_index

subroutine choose_boundary_magnetic_surfaces_for_the_mapping(r_inner_surf,z_inner_surf,r_outer_surf,z_outer_surf)
  use constants, only: p_, two
  use magnetic_coordinates,only: pfn_inner, pfn_bdry !boundary for the region in which magnetic_coordinates are constructed
  use radial_module, only: r_axis, z_axis, psi_axis, psi_lcfs
  use boundary, only: x_lcfs, z_lcfs, np_lcfs
  use contour_mod, only: find_contour
  implicit none
  real(p_), intent(out):: r_inner_surf(np_lcfs), z_inner_surf(np_lcfs), r_outer_surf(np_lcfs), z_outer_surf(np_lcfs)
  real(p_) :: psi_val, psi_bdry, psi_inner

  psi_inner = psi_axis + pfn_inner*(psi_lcfs-psi_axis) 
  psi_val = (psi_axis+psi_inner)/two !the inner boundary for doing the mapping
  call find_contour(psi_val, r_inner_surf, z_inner_surf)

  psi_bdry = psi_axis+pfn_bdry*(psi_lcfs-psi_axis)
  psi_val = (psi_lcfs+psi_bdry)/two !the outer boundary for doing the mapping, which is chosen to be between lcfs and the outer boundary of the region in which magnetic_coordinates are availabe
  call find_contour(psi_val, r_outer_surf, z_outer_surf)

end subroutine choose_boundary_magnetic_surfaces_for_the_mapping

subroutine create_cylindrical_grids(r_outer_surf,z_outer_surf,np_lcfs,nx,nz,r,z,dr,dz,i0,j0)
  !create rectangular box (with cylindrical grids in it) in the poloidal plane with the boundary flux surface within the box and (r_axis,z_axis) is exactly on a grid point
  use constants, only: p_
  !  use magnetic_coordinates,only:r_mc, z_mc,mpol,nrad
  use radial_module, only: r_axis,z_axis
  implicit none
  integer, intent(in) :: np_lcfs, nx, nz
  real(p_), intent(out) :: r(nx), z(nz), dr, dz
  integer, intent(out) :: i0, j0 ! RZ index of the magnetic axis
  real(p_) :: r_min, r_max, z_min, z_max, r_width, z_width
  !  real(p_):: r_mag_surf1(mpol), z_mag_surf1(mpol)
  real(p_) :: r_outer_surf(np_lcfs), z_outer_surf(np_lcfs)
  integer :: i,j,nxp,nzp
  real(p_) :: rp(nx-1),zp(nz-1)

  nxp = nx-1 !using a reduced number, so that I can append the array with an additional element
  nzp = nz-1 

!!$  do i=1,mpol !select the boundary magnetic surface
!!$     r_mag_surf1(i)=r_mc(i,nrad)
!!$     z_mag_surf1(i)=z_mc(i,nrad)
!!$  enddo

  r_min = minval(r_outer_surf)
  r_max = maxval(r_outer_surf)
  r_width = r_max - r_min
  z_min = minval(z_outer_surf)
  z_max = maxval(z_outer_surf)
  z_width = z_max - z_min
  dr=r_width/(nxp-1)
  dz=z_width/(nzp-1)

  i0 = floor((r_axis-r_min)/dr)+1 !i index near the magnetic axis
  j0 = floor((z_axis-z_min)/dz)+1 !j index near the magnetic axis
  rp(i0) = r_axis !set (i0,j0) to be on the magnetic axis, since there is a floor operation when getting (i0,j0), this shifts the box to the low-field-side (lfs) and upward
  zp(j0) = z_axis 

  !I do the above steps because I want the magnetic axis to lie exactly on a grid so that I can get grids that exactly represent the midplane.
  do i=i0-1,1,-1 !the grids at the high-field-side (hfs) of the magnetic axis
     rp(i)=rp(i+1)-dr
  enddo
  do i=i0+1,nxp,1 !the girds at the low-field-side (lfs) of the magnetic axis
     rp(i)=rp(i-1)+dr
  enddo

  do j=j0-1,1,-1 !the grids below the midplane
     zp(j)=zp(j+1)-dz
  enddo
  do j=j0+1,nzp,1 !the grids above the midplane
     zp(j)=zp(j-1)+dz
  enddo

  !add an additional element. This is needed because setting the magnetic axis to be on a grid usually shifts the box to the lfs, which will make some points on the hfs not included in the box
  r(1) = rp(1) - dr
  do i = 1, nxp
     r(i+1) = rp(i)
  enddo
  z(1) = zp(1) - dz !similar reason as mentioned above
  do j = 1, nzp
     z(j+1) = zp(j)
  enddo

  i0 = i0+1 !the i index of the magnetic axis at the new array r
  j0 = j0+1 !the j index of the magnetic axis at the new array z
  block
    use domain_decomposition, only : myid
    integer :: u
    if(myid==0) then
       open(newunit=u,file='rzgrid.txt')
       do i =1,nx
          do j =1,nz
             write(u,*) r(i), z(j)
          enddo
       enddo
       close(u)
    end if
  endblock
end subroutine create_cylindrical_grids


module map_to_mc
contains
  subroutine interpolate_from_cylindrical_to_magnetic_coordinates(r0, z0, theta0, tor_shift0)
    !calculate the magnetic coordinats (theta0,tor_shift0) of (R0,Z0) by using interpolation,
    !radcor0 is not calculated in this subroutine because I can use another reliable way to calculate radcor0 from (r0,z0), i.e., radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))
    !to use this iterpolation, first to make sure that (r0,z0) is within the specfied region
    use constants,only:p_,two,twopi
    use mapping_module,only: r_cyl,z_cyl,dr,dz,radcor,theta_a,theta_b,tor_shift_a,i0,j0,tor_shift_b
    use interpolate_module,only: linear_2d_interpolate_kernel
    implicit none
    real(p_),intent(in):: r0,z0
    real(p_),intent(out):: theta0,tor_shift0
    real(p_):: radcor_tmp(2,2),theta_tmp(2,2),tor_shift_tmp(2,2),qval
    integer::i,j,ii,jj
    real(p_):: radcor_as_func_of_pfn,pfn_func

    i=floor((r0-r_cyl(1))/dr)+1
    j=floor((z0-z_cyl(1))/dz)+1

    !  2D interpolations to get radcor0, theta0, and tor_shift0
!!$  do ii=1,2
!!$     do jj=1,2
!!$        radcor_tmp(ii,jj)=radcor(i+ii-1,j+jj-1)
!!$     enddo
!!$  enddo
!!$
!!$  call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),radcor_tmp,r0,z0,radcor0)

    !  radcor0=radcor_as_func_of_pfn(pfn_func(r0,z0))


    do ii=1,2
       do jj=1,2
          theta_tmp(ii,jj)=theta_a(i+ii-1,j+jj-1)
       enddo
    enddo

!!$ if(i>i0 .and. j+1.eq.j0) then !handle the boundary at the low-field-side midplane
!!$theta_tmp(1,2)=twopi 
!!$theta_tmp(2,2)=twopi 
!!$endif
    if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side (theta cut)
       theta_tmp(1,1)=theta_b(i,j0)
       theta_tmp(2,1)=theta_b(i+1,j0)
    endif
    call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),theta_tmp,r0,z0,theta0)

    do ii=1,2
       do jj=1,2
          tor_shift_tmp(ii,jj)=tor_shift_a(i+ii-1,j+jj-1)
       enddo
    enddo

!!$if(i>i0 .and. j+1.eq.j0) then !handle the boundary at the low-field-side midplne
!!$!qval=2.3074808101594884*1.0004
!!$!tor_shift_tmp(1,2)=two*tor_shift(i,j)-tor_shift(i,j-1) 
!!$!tor_shift_tmp(1,2)=qval*twopi
!!$!tor_shift_tmp(2,2)=qval*twopi
!!$tor_shift_tmp(1,2)=tor_shift_b(i,j0)
!!$tor_shift_tmp(2,2)=tor_shift_b(i+1,j0)
!!$endif

    if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side midplane (theta cut)
       tor_shift_tmp(1,1)=tor_shift_b(i,j0)
       tor_shift_tmp(2,1)=tor_shift_b(i+1,j0)
    endif

    call linear_2d_interpolate_kernel(r_cyl(i),z_cyl(j),tor_shift_tmp,r0,z0,tor_shift0)

    !alpha0=phi0-tor_shift0

  end subroutine interpolate_from_cylindrical_to_magnetic_coordinates


  pure  subroutine interpolate_from_cylindrical_to_magnetic_coordinates1(r0, z0, theta0)
    use constants, only: p_
    use mapping_module, only: r_cyl, z_cyl, dr, dz, radcor, theta_a, theta_b, i0, j0
    use interpolate_module
    implicit none
    real(p_), intent(in) :: r0, z0
    real(p_), intent(out) :: theta0
    real(p_) :: tmp(2,2)
    integer :: i, j, ii, jj

    i = floor((r0-r_cyl(1))/dr)+1
    j = floor((z0-z_cyl(1))/dz)+1

    do ii=1,2
       do jj=1,2
          tmp(ii,jj) = theta_a(i+ii-1,j+jj-1)
       enddo
    enddo

    if(i<i0 .and. j.eq.j0) then !handle the boundary at the high-field-side (theta cut)
       tmp(1,1) = theta_b(i,j0)
       tmp(2,1) = theta_b(i+1,j0)
    endif
    
    call linear_2d_interpolate_kernel(r_cyl(i), z_cyl(j), tmp, r0, z0, theta0)

  end subroutine interpolate_from_cylindrical_to_magnetic_coordinates1
end module map_to_mc
module my_FFTW3
  use constants, only: p_
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  save
  type(C_PTR) :: plan_toroidal, plan_radial, plan_toroidal_backward, plan_radial_backward
  type(C_PTR) :: plan_sin, plan_sin_backward

  complex(p_), allocatable :: in1(:), out1(:)
  complex(p_), allocatable :: in2(:), out2(:)
  complex(p_), allocatable :: in4(:), out4(:)           
  real(p_), allocatable :: in3(:), out3(:)

contains
  subroutine  initialize_fft
    use magnetic_coordinates,only: m=>mtor, n=>nrad, mpol2
    allocate(in1(0:m-1))
    allocate(out1(0:m-1))           
    allocate(in2(0:n-1))
    allocate(out2(0:n-1))
    allocate(in4(0:mpol2-1))
    allocate(out4(0:mpol2-1))

    call dfftw_plan_dft_1d(plan_toroidal, m, in1, out1, FFTW_FORWARD, FFTW_measure)
    call dfftw_plan_dft_1d(plan_radial, n, in2, out2, FFTW_FORWARD, FFTW_MEASURE)
    call dfftw_plan_dft_1d(plan_toroidal_backward, m, in1, out1, FFTW_backward, FFTW_measure)
    plan_radial_backward = fftw_plan_dft_1d(n, in2, out2, FFTW_backward, FFTW_measure)

    allocate(in3(n-2)) 
    allocate(out3(n-2))
    plan_sin = fftw_plan_r2r_1d(n-2, in3, out3, FFTW_RODFT00, FFTW_MEASURE)

  end subroutine initialize_fft

end module my_FFTW3

module transform_module
implicit none

contains
  subroutine oned_fourier_transform1(s, s_fft, m, n) !calculating 1d DFT of s(:,:) along the first dimension
    use, intrinsic :: iso_c_binding
    use my_FFTW3, only: plan_toroidal, in1,out1
    use constants,only:p_
    implicit none
    include 'fftw3.f03'
    integer, intent(in) :: m, n
    real(p_), intent(in) :: s(0:m-1, n)
    complex(p_), intent(out) :: s_fft(0:m-1, n)
    integer :: j

    do j = 1, n  ! using openmp to parallize this loop gave wrong results (on Edison machine), strange!, possibly indicating fftw_execute_dft is not thread-safe.
       in1(:) = s(:,j) !copy in, meanwhile convert real array to complex array
       call fftw_execute_dft(plan_toroidal, in1(:), out1(:))     !Fourier transformation along the first dimension
       s_fft(:,j) = out1(:) !copy out
    enddo

!tests indicate the following is slower than the above, do not know why.
!!$    s_complex=s !convert real array to complex array
!!$    do j=1,n
!!$       call fftw_execute_dft(plan_toroidal, s_complex(:,j), s_fft(:,j))     !Fourier transformation along the first dimension
!!$    enddo

  end subroutine oned_fourier_transform1

  subroutine oned_DFT_parallel_version(s, s_fft, m, n) !DFT of s(:,:) along the 1st dimension
   !when call this subroutine using different input, make sure that n is identical because the "save" used here
    use, intrinsic :: iso_c_binding
    use my_FFTW3,only: plan_toroidal, in1, out1
    use constants,only:p_
    use domain_decomposition,only: TCLR, ntube, grid_comm
    use mpi
    implicit none
    include 'fftw3.f03'
    integer, intent(in) :: m,n
    real(p_), intent(in) :: s(0:m-1,n)
    complex(p_), intent(out) :: s_fft(0:m-1,n)
    integer:: i,j,ierr
    logical, save :: is_first=.true.
    integer, save :: my_part_start, my_part_end, spacing, my_range
    integer, allocatable, save:: recvcounts(:), displacement(:)
    complex(p_), allocatable,save:: my_s_fft(:,:)

    if (is_first.eqv..true.) then !needs to be genreated only for the first time and used multiple times
       spacing=n/ntube
       my_part_start=TCLR*spacing+1
       my_part_end=my_part_start+(spacing-1)
       if(TCLR.eq.(ntube-1)) my_part_end=n !the last process handles all the remainder part, in the case that n is not a perfect multiple of ntube
       my_range=my_part_end-my_part_start+1
       allocate(my_s_fft(0:m-1,my_range))
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:)=spacing*m
       recvcounts(ntube-1)= recvcounts(ntube-1)+(n-spacing*ntube)*m !last process contains additional elements
       do i=0,ntube-1
          displacement(i)=i*spacing*m
       enddo
       is_first=.false.
    endif

    do j = my_part_start,my_part_end  ! using openmp to parallize this loop gave wrong results (on Edison machine), strange!, possibly indicating fftw_execute_dft is not thread-safe.
       in1(:)=s(:,j) !copy in, meanwhile convert real array to complex array
       call fftw_execute_dft(plan_toroidal, in1(:), out1(:))     !FFT along the 1st dimension
       my_s_fft(:,j-my_part_start+1) = out1(:) !copy out
    enddo

    if(ntube.gt.1) then
       call mpi_allgatherv(my_s_fft, m*my_range, MPI_DOUBLE_COMPLEX, &
            &     s_fft, recvcounts, displacement, MPI_DOUBLE_COMPLEX, grid_comm, ierr)
    else
       s_fft = my_s_fft
    endif
  end subroutine oned_DFT_parallel_version

  subroutine oned_backward_DFT_parallel_version(field_fft, field, m, n) !radial task decomposion at a z gridpoint
    !when call this subroutine using different input (field_fft), make sure that n is identical because the "save" property used here
    use, intrinsic :: iso_c_binding
    use my_FFTW3, only: plan_toroidal_backward
    use constants,only:p_
    use domain_decomposition,only: TCLR, ntube, grid_comm
    use mpi
    implicit none
    include 'fftw3.f03' 
    integer, intent(in) :: m, n
    complex(p_), intent(in) :: field_fft(0:m-1,n)
    real(p_), intent(out) :: field(0:m-1,n)
    complex(p_) :: in1(0:m-1), out1(0:m-1)  
    integer :: i, j, ierr
    logical, save :: is_first = .true.
    integer, save :: my_part_start, my_part_end, spacing, my_range
    integer, allocatable, save:: recvcounts(:), displacement(:)
    real(p_), allocatable, save :: my_field(:,:)

    if (is_first.eqv..true.) then !needs to be genreated only for the first time and used multiple times
       is_first = .false.
       spacing = n/ntube
       my_part_start = TCLR*spacing+1
       my_part_end = my_part_start+(spacing-1)
       !the last process handles all the remainder part, in the case that n is not a perfect multiple of ntube
       if(TCLR.eq.(ntube-1)) my_part_end=n 
       my_range=my_part_end-my_part_start+1
       allocate(my_field(0:m-1,my_range))
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:)=spacing*m
       recvcounts(ntube-1)= recvcounts(ntube-1)+(n-spacing*ntube)*m !last process contains additional elements
       do i=0,ntube-1
          displacement(i)=i*spacing*m
       enddo
    endif

    do j = my_part_start, my_part_end  !each processor handles its radial range
       in1(:) = field_fft(:,j) !copy in
       call fftw_execute_dft(plan_toroidal_backward, in1(:), out1(:)) !fourier transform along toroidal direction
       my_field(:,j-my_part_start+1) = real(out1(:))/m !copy out and inclued the 1/m factor
    enddo

    if(ntube.gt.1) then
       call mpi_allgatherv(my_field, m*my_range, MPI_DOUBLE, &
            &  field, recvcounts, displacement, MPI_DOUBLE, grid_comm, ierr)
    else
       field = my_field
    endif
  end subroutine oned_backward_DFT_parallel_version


  subroutine oned_backward_fourier_transform1(field_fft, field, m, n) 
    use, intrinsic :: iso_c_binding
    use my_FFTW3, only: plan_toroidal_backward
    use constants, only: p_
    implicit none
    include 'fftw3.f03' 
    integer, intent(in) :: m, n
    complex(p_), intent(in) :: field_fft(0:m-1, n)
    real(p_), intent(out) :: field(0:m-1, n)
    complex(p_) :: in(0:m-1, n), out(0:m-1, n)  
    integer :: j

    in = field_fft
    do j = 1, n ! DFT along the first direction
       call fftw_execute_dft(plan_toroidal_backward, in(:,j), out(:,j))
    enddo

    field = real(out)/m
  end subroutine oned_backward_fourier_transform1



  subroutine oned_fourier_transform2(s,s_fft,m,n) !calculating 1d DFT of s(:,:) along the second dimension
    use my_FFTW3
    use constants,only:p_
    implicit none

    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    type(C_PTR) :: plan
    integer:: i

    !Fourier transformation along the second dimension
    in=s
    plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan, in(i,:), out(i,:))
    enddo
    call fftw_destroy_plan(plan)  
    s_fft=out
  end subroutine oned_fourier_transform2


  subroutine twod_fourier_transform(s,s_fft,m,n) !calculating 2d DFT of source, tested
    use my_FFTW3
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    integer:: i,j

    !Fourier transformation along the first dimension
    in=s
    !plan_toroidal = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_FORWARD,FFTW_ESTIMATE)
    do j=0,n-1
       call fftw_execute_dft(plan_toroidal, in(:,j), out(:,j))
    enddo
    !call fftw_destroy_plan(plan_toroidal)  

    !  in=out

    !  call toroidal_filter(in,out,m,n)

!!$  write(*,*) 'given by fftw'
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j)),i=0,m-1)
!!$  enddo
!!$do j=0,n-1
!!$call my_fft(in(:,j),out(:,j),m)
!!$enddo
!!$ write(*,*) 'given by my_fft'
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j)),i=0,m-1)
!!$  enddo

    !Fourier transformation along the second dimension
    in=out
    !plan_radial = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan_radial, in(i,:), out(i,:))
    enddo
    !call fftw_destroy_plan(plan_radial)  
    s_fft=out

  end subroutine twod_fourier_transform


subroutine twod_fourier_transform_xz(s,s_fft,m,n) !calculating 2d DFT of source, tested
    use my_FFTW3
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(in):: s(0:m-1,0:n-1)
    complex(p_),intent(out):: s_fft(0:m-1,0:n-1)
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    integer:: i,j
type(C_PTR) :: plan_z
    !Fourier transformation along the first dimension
    in=s
    do j=0,n-1
       call fftw_execute_dft(plan_radial, in(:,j), out(:,j))
    enddo

    !Fourier transformation along the second dimension
    in=out
    plan_z = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_FORWARD,FFTW_ESTIMATE)
    do i=0,m-1
       call fftw_execute_dft(plan_z, in(i,:), out(i,:))
    enddo
    call fftw_destroy_plan(plan_z)
    s_fft=out

  end subroutine twod_fourier_transform_xz


  subroutine twod_inverse_fourier_transform(field_fft,field,m,n) !tested
    use my_FFTW3
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    !  complex(C_DOUBLE_COMPLEX),intent(in):: field_fft(0:m-1,0:n-1)
    complex(p_),intent(in):: field_fft(0:m-1,0:n-1)
    real(p_),intent(out):: field(0:m-1,0:n-1)
    !  complex(C_DOUBLE_COMPLEX) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
    integer:: i,j

    in=field_fft
    !plan_toroidal_backward = fftw_plan_dft_1d(m, in(:,1),out(:,1), FFTW_backward,FFTW_ESTIMATE)
    do j=0,n-1 !fourier transform along the first direction
       call fftw_execute_dft(plan_toroidal_backward, in(:,j), out(:,j))
    enddo
    !call fftw_destroy_plan(plan_toroidal_backward)  

    in=out
    !plan_radial_backward = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_backward,FFTW_ESTIMATE)
    do i=0,m-1 !fourier transform along the second direction
       call fftw_execute_dft(plan_radial_backward, in(i,:), out(i,:))
    enddo
    !call fftw_destroy_plan(plan_radial_backward)  

    field=real(out)/(m*n)
!!$  !check whether the results are the same
!!$  write(*,*) 'original data='
!!$  do j=0,n-1
!!$     write(*,*) (s(i,j),i=0,m-1)
!!$  enddo
!!$  write(*,*) 'DFT+backward-DFT data='
!!$  do j=0,n-1
!!$     write(*,*) (real(out(i,j))/(m*n),i=0,m-1)
!!$  enddo
  end subroutine twod_inverse_fourier_transform


subroutine simple_fft(in,out,m) !for test
  use constants,only:p_
  use constants,only:twopi
  implicit none
  integer,intent(in):: m
  complex(p_),intent(in)::in(0:m-1)
  complex(p_),intent(out)::out(0:m-1)
  complex(p_),parameter::ii=(0.0_p_,1._p_)
  complex(p_):: sum
  integer:: i,ip

  do i=0,m-1
     sum=0._p_
     do ip=0,m-1
        sum=sum+in(ip)*exp(-twopi*ii/m*ip*i)
     enddo
     out(i)=sum
  enddo

end subroutine simple_fft

subroutine my_inverse_fft(in,out,m) !for test
  use constants,only:p_
  use constants,only:twopi
  implicit none
  integer,intent(in):: m
  complex(p_),intent(in)::in(0:m-1)
  complex(p_),intent(out)::out(0:m-1)
  complex(p_),parameter::ii=(0.0_p_,1._p_)
  complex(p_):: sum
  integer:: i,ip

  do i=0,m-1
     sum=0._p_
     do ip=0,m-1
        sum=sum+in(ip)*exp(twopi*ii/m*ip*i)
     enddo
     out(i)=sum
  enddo
out=out/m
end subroutine my_inverse_fft





subroutine oned_backward_fourier_transform2(field_fft,field,m,n) 
    use my_FFTW3
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  complex(p_),intent(in):: field_fft(0:m-1,0:n-1)
  real(p_),intent(out):: field(0:m-1,0:n-1)
  type(C_PTR) :: plan
  complex(p_) :: in(0:m-1,0:n-1), out(0:m-1,0:n-1)  
  integer:: i

  in=field_fft
  plan = fftw_plan_dft_1d(n, in(1,:),out(1,:), FFTW_backward,FFTW_ESTIMATE)
  do i=0,m-1 !fourier transform along the second direction
     call fftw_execute_dft(plan, in(i,:), out(i,:))
  enddo
  call fftw_destroy_plan(plan)  
  field=real(out)/n
end subroutine oned_backward_fourier_transform2

subroutine oned_sine_transform2(s, s_dst, m, n) !calculating 1d DST of s(:,:) along the second dimension
  use my_FFTW3, only: plan_sin, in3, out3
  use constants,only:p_
  implicit none
  integer, intent(in) :: m, n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  real(p_), intent(out) :: s_dst(0:m-1,0:n-1)
 ! type(C_PTR) :: plan !a pointer, which needs to be distoried, otherwise may cause memory leak
  integer:: i
  !plan =fftw_plan_r2r_1d(n, in(1,:), s_dst(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i = 0, m-1
     in3(:) = s(i, :)
     !call dfftw_execute_r2r(plan_sin,s(i,:),s_dst(i,:)) !DST along the second dimension
     call dfftw_execute_r2r(plan_sin, in3, out3) !DST along the second dimension
     s_dst(i, :) = out3(:)
  enddo
!  call fftw_destroy_plan(plan)  

end subroutine oned_sine_transform2

subroutine oned_inverse_sine_transform2(s_dst,s,m,n) !computing 1d inverse DST of s(:,:) along the second dimension
  use my_FFTW3, only: plan_sin, in3, out3
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: s_dst(0:m-1,0:n-1)
  real(p_),intent(out):: s(0:m-1,0:n-1)
!  real(p_) :: in(0:m-1,0:n-1)
!  type(C_PTR) :: plan
  integer:: i

  !in=s_dst
  !plan =fftw_plan_r2r_1d(n, in(1,:), s(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     in3(:)=s_dst(i,:)
    ! call dfftw_execute_r2r(plan_sin,s_dst(i,:),s(i,:))   !DST along the second dimension
     call dfftw_execute_r2r(plan_sin,in3,out3)   !DST along the second dimension
     s(i,:)=out3(:)
  enddo
  !call fftw_destroy_plan(plan)  
  s=s/(2*(n+1))
end subroutine oned_inverse_sine_transform2


subroutine dst_dft(s,s_spectrum,m,n)
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  implicit none

  include 'fftw3.f03'
  integer,intent(in):: m,n
  real(p_),intent(in):: s(0:m-1,0:n-1)
  complex(p_),intent(out):: s_spectrum(0:m-1,0:n-1)
  real(p_):: in1(0:m-1,0:n-1),s_dst(0:m-1,0:n-1)
  complex(p_):: in2(0:m-1,0:n-1)
  type(C_PTR) :: plan
  integer:: i,j

  in1=s
  plan =fftw_plan_r2r_1d(n,   in1(1,:), s_dst(1,:), FFTW_RODFT00, FFTW_ESTIMATE)
  do i=0,m-1
     call dfftw_execute_r2r(plan,in1(i,:), s_dst(i,:)) !DST along the second dimension
  enddo
  call fftw_destroy_plan(plan)  

  in2=s_dst
  plan = fftw_plan_dft_1d(m,  in2(:,1), s_spectrum(:,1), FFTW_FORWARD,FFTW_ESTIMATE)
  do j=0,n-1
     call fftw_execute_dft(plan, in2(:,j), s_spectrum(:,j))   !DFT along the first dimension
  enddo
  call fftw_destroy_plan(plan)
!s_spectrum=s_spectrum/((n+1)*m) !the corresponding expansion coeficient is the dst_dft devided by (n+1)*m, see my notes on Fourier analysis
end subroutine dst_dft


end module transform_module


module fourn_module
contains
  subroutine twod_fourier_transform_nr(s,s_fft,m,n) !tested on Feb 27, 2018, note that the convension of FFT and inverse FFT definition is differnt between FFTW and Numerical-recipe book. I forgot about this when comparing the result given by this subroutine with that by FFTW subroutine, and it took about one hour for me to finally figure out why they are different.
    use constants,only:p_
!    use transform_module,only: my_fft
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(in):: s(m,n)
    complex(p_),intent(out):: s_fft(m,n)
    real(p_):: data(2*m*n)
    integer,parameter:: ndim=2
    complex(p_),parameter:: ii=(0._p_,1._p_)
    complex(p_):: in(m,n),out(m,n)
    integer:: i,j,k,nn(ndim),isign

    do j=1,n
       do i=1,m
          k=i+m*(j-1)
          data(2*k-1)=s(i,j)  !real part
          data(2*k)=0._p_ !imaginary part
       enddo
    enddo

!---for testing--------passed successfully
!!$    in=s
!!$    do j=1,n
!!$       call my_fft(in(:,j),out(:,j),m)
!!$    enddo
!!$    !Fourier transformation along the second dimension
!!$    in=out
!!$    do i=1,m
!!$       call my_fft(in(i,:),out(i,:),n)
!!$    enddo
!!$    s_fft=out
!----for testing----
    nn(1)=m
    nn(2)=n
    isign=-1 !isign=-1 corresponds to inverse FFT in numerical-recipe-book, but corresponds to the forward FFT in FFTW
    call fourn(data,nn,ndim,isign)

    do j=1,n
       do i=1,m
          k=i+m*(j-1)
          s_fft(i,j)=data(2*k-1)+ii*data(2*k)
       enddo
    enddo
 end subroutine twod_fourier_transform_nr


subroutine twod_inverse_fourier_transform_nr(s_fft,s,m,n) 
  use constants,only:p_
  implicit none
  integer,intent(in):: m,n
  complex(p_),intent(in):: s_fft(m,n)
  real(p_),intent(out):: s(m,n)
  real(p_):: data(2*m*n)
  integer,parameter:: ndim=2
  integer:: i,j,k,nn(ndim),isign

  do j=1,n
     do i=1,m
        k=i+m*(j-1)
        data(2*k-1)=real(s_fft(i,j))  !real part
        data(2*k)=imag(s_fft(i,j)) !imaginary part
     enddo
  enddo

  nn(1)=m
  nn(2)=n
  isign=1 !isign=1 corresponds to FFT in numerical-recipe-book, but corresponds the inverse FFT in FFTW
  call fourn(data,nn,ndim,isign)
  do i=1,m
     do j=1,n
        k=i+m*(j-1)
        s(i,j)=data(2*k-1)
     enddo
  enddo
s=s/(m*n)
end subroutine twod_inverse_fourier_transform_nr

SUBROUTINE fourn(data,nn,ndim,isign)
  use constants,only:p_
  implicit none
  INTEGER isign,ndim,nn(ndim)
  REAL(p_):: data(*)
  !Replaces data by its ndim -dimensional discrete Fourier transform, if isign is input as 1.
  !nn(1:ndim) is an integer array containing the lengths of each dimension (number of complex values), which MUST all be powers of 2. 
  !data is a real array of length twice the product of these lengths, in which the data are stored as in a multidimensional complex FORTRAN array
  !If isign is input as -1, data is replaced by its inverse transform times the product of the lengths of all dimensions. 
  !code from the numerical recipe book
  INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,&
       & ip3,k1,k2,n,nprev,nrem,ntot
  REAL(p_):: tempi,tempr
  real(p_):: theta,wi,wpi,wpr,wr,wtemp !Double precision for trigonometric re-currences.

  ntot=1
  do idim=1,ndim !Compute total number of complex values.
     ntot=ntot*nn(idim)
  enddo
  nprev=1
  do idim=1,ndim !Main loop over the dimensions.
     n=nn(idim)
     nrem=ntot/(n*nprev)
     ip1=2*nprev
     ip2=ip1*n
     ip3=ip2*nrem
     i2rev=1
     do  i2=1,ip2,ip1 !This is the bit-reversal section of the routine.
        if(i2.lt.i2rev)then
           do  i1=i2,i2+ip1-2,2
              do  i3=i1,ip3,ip2
                 i3rev=i2rev+i3-i2
                 tempr=data(i3)
                 tempi=data(i3+1)
                 data(i3)=data(i3rev)
                 data(i3+1)=data(i3rev+1)
                 data(i3rev)=tempr
                 data(i3rev+1)=tempi
              enddo
           enddo
        endif
        ibit=ip2/2
1       if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
           i2rev=i2rev-ibit
           ibit=ibit/2
           goto 1
        endif
        i2rev=i2rev+ibit
     enddo
     ifp1=ip1 !Here begins the Danielson-Lanczos section of the routine.
2    if(ifp1.lt.ip2)then
        ifp2=2*ifp1
        theta=isign*6.28318530717959d0/(ifp2/ip1) !Initialize for the trig. recur-rence.
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do  i3=1,ifp1,ip1
           do  i1=i3,i3+ip1-2,2
              do  i2=i1,ip3,ifp2
                 k1=i2 !Danielson-Lanczos formula:
                 k2=k1+ifp1
                 tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                 tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                 data(k2)=data(k1)-tempr
                 data(k2+1)=data(k1+1)-tempi
                 data(k1)=data(k1)+tempr
                 data(k1+1)=data(k1+1)+tempi
              enddo
           enddo
           wtemp=wr !Trigonometric recurrence.
           wr=wr*wpr-wi*wpi+wr
           wi=wi*wpr+wtemp*wpi+wi
        enddo
        ifp1=ifp2
        goto 2
     endif
     nprev=n*nprev
  enddo
  return
END SUBROUTINE fourn
end module fourn_module
module fk_module !fully-kinetic for ions
  use constants,only:p_, zero
  implicit none
  save
  real(p_) :: mass_i,charge_i,dtao_fk
  real(p_) :: ni0,ti0,kappa_ni,kappa_ti !ti0 is a typical value of ion temperature in the compuational region, in flux tube model, ti0 is the temperature at the reference magnetic surface
  real(p_) :: vt_i,vmin_i,vmax_i
  real(p_) ::normalizing_factor
  real(p_) :: omegan_fk, tn_fk, vn_fk
  
  integer :: ntouch_bdry_i=0, total_ntouch_bdry_i !numbe of markers that touch the boundary in each process and all processes, respectiveyl
  integer :: total_nmarker_i  ! total number of ion markers (including the particles in all the processors).
  integer :: nmarker_i_per_cell  ! total number of ion markers (including the particles in all the processors).
  integer :: nmarker_i !particle number in a single processor, its value will be differnt for differnt processors and at differnt time

  real(p_),allocatable :: ps_vol_i(:) !determined by the initial loading, is constant for each marker over the time-evoultion
  real(p_),allocatable :: w_i(:)  !weight of ion markers
  real(p_),allocatable :: w_i_mid(:)  !weight of ion markers at t_{n+1/2}
  real(p_),allocatable :: w_i_star(:)  !weight of ion markers

  real(p_),allocatable :: r_i(:),z_i(:),phi_i(:) !Cylindrical coordinates at time t_{n},  unit Ln, rad
  real(p_),allocatable :: r_i_old(:),z_i_old(:),phi_i_old(:) !Cylindrical coordinates at time t_{n},  unit Ln, rad, temporary working arrays for computing averaing
  real(p_),allocatable :: r_i_mid(:),z_i_mid(:),phi_i_mid(:) !Cylindrical coordinates at time t_{n+1/2}, unit Ln, rad, 
  real(p_),allocatable :: radcor_i(:),theta_i(:), alpha_i(:),tor_shift_i(:) !magnetic coordinates at integer-time-step, alpha is the generalized toroidal angle
  real(p_),allocatable :: radcor_i_mid(:),theta_i_mid(:),alpha_i_mid(:) !magnetic coordinates at half time-step

  real(p_),allocatable :: vr_i(:),vz_i(:),vphi_i(:) ! projection of velocity at t_{n-1/2} to the basis cylindrical vector at integer-time-step t_{n}, unit vn_fk=ln/tn_fk, 
  real(p_),allocatable :: vr_i_old(:),vz_i_old(:),vphi_i_old(:) !temporary working arrays for computing averaing
  real(p_),allocatable :: vr_i_integer_mid(:),vz_i_integer_mid(:),vphi_i_integer_mid(:) !projection of velocity at t_{n} to the basis cylindrical vector at t_{n+1/2}, initial condition for the second boris pusher, unit vn_fk=ln/tn_fk, 

  real(p_),allocatable :: vr_i_mid(:),vz_i_mid(:),vphi_i_mid(:) !projection of velocity at t_{n+1/2} to the local basis vectors at t_{n+1/2}
  real(p_),allocatable :: vr_i_integer(:),vz_i_integer(:),vphi_i_integer(:) !projection of velocity at t_{n} to the local basis vectors at t_{n}

  real(p_),allocatable :: v_i(:),vpar_i(:),vx_i(:),vy_i(:) !vx is defined by vx=v_dot_grad_x, vy is defined by vy=v_dot_grad_y,note that grad_x and grad_y are not perpendicular to each other
  real(p_),allocatable :: grad_psi_i(:),grad_alpha_i(:), grad_psi_dot_grad_alpha_i(:),bval_i(:)
  real(p_),allocatable :: v_i_mid(:),vpar_i_mid(:),vx_i_mid(:),vy_i_mid(:)
  real(p_),allocatable :: grad_psi_i_mid(:),grad_alpha_i_mid(:), grad_psi_dot_grad_alpha_i_mid(:),bval_i_mid(:)

  logical,allocatable :: touch_bdry_i(:),active_i(:) !indicates whether the orbit of a marker touches the boundary
  logical,allocatable :: touch_bdry_i_mid(:),active_i_mid(:)
  real(p_),dimension(:,:),allocatable :: my_den_i_left, my_den_i_right !fk ion density

  integer  :: ion_spatial_loading_scheme  !1=>uniform in (psi,theta,alpha) coordinates; 2=>uniform in real space
  integer  :: ion_velocity_loading_scheme !1=>uniform in (v,theta_v,phi_v) coordinates; 2=>Isotropic Gaussian in (v,theta_v,phi_v)
  integer :: fk_nonlinear
contains
  subroutine initialize_fk()
    use domain_decomposition,only: numprocs, myid
    use magnetic_coordinates, only : nrad, mpol2,mtor
    namelist/fk_nmlt/mass_i,charge_i,ni0,ti0,kappa_ni,kappa_ti,nmarker_i_per_cell, &
         & ion_spatial_loading_scheme, ion_velocity_loading_scheme, fk_nonlinear
    integer:: fixed_large_size, u

     open(newunit=u,file='input.nmlt')
     read(u,fk_nmlt)
     close(u)
     if(myid==0)  write(*,fk_nmlt)

    
    total_nmarker_i=nmarker_i_per_cell*nrad*mpol2*mtor
    if(myid.eq.0) write(*,*) 'total number of ions=',     total_nmarker_i
    nmarker_i=total_nmarker_i/numprocs !nmarker_i initially store the number of markers initially loaded per processor (i.e.total_nmarker_i/numprocs), latter actual number of markers per proc will be assigned to nmarker_i, the value of which will be differnt for differnt processors and at differnt time
    fixed_large_size=(total_nmarker_i/numprocs)*3/2 !the number of particle per proc after re-arranging the particles between the processors may exceed the number of original loaded particles per proc (i.e., total_nmarker_i/numprocs), increasing the array length by a factor of 3/2 is needed to make sure that the array is big enough to contain all the particles that belong to the domain for which the processor is responsible.
    !  write(*,*), 'nmarker_i, fixed_large_size=',nmarker_i, fixed_large_size

    allocate(radcor_i(fixed_large_size)) 
    allocate(theta_i(fixed_large_size))
    allocate(alpha_i(fixed_large_size))  
    allocate(tor_shift_i(fixed_large_size))

    allocate(radcor_i_mid(fixed_large_size)) 
    allocate(theta_i_mid(fixed_large_size))
    allocate(alpha_i_mid(fixed_large_size))

    allocate(v_i(fixed_large_size))
    allocate(vr_i(fixed_large_size))
    allocate(vz_i(fixed_large_size))
    allocate(vphi_i(fixed_large_size))

    allocate(r_i(fixed_large_size))
    allocate(z_i(fixed_large_size))
    allocate(phi_i(fixed_large_size))
    allocate(r_i_mid(fixed_large_size)) 
    allocate(z_i_mid(fixed_large_size)) 
    allocate(phi_i_mid(fixed_large_size))

    allocate(w_i(fixed_large_size)) 
    allocate(w_i_mid(fixed_large_size)) 
    allocate(w_i_star(fixed_large_size)) 

    allocate(ps_vol_i(fixed_large_size))
    allocate(active_i(fixed_large_size)) !whether particles are within computational boundary
    allocate(active_i_mid(fixed_large_size)) !whether particles are within computational boundary
    allocate(touch_bdry_i(fixed_large_size)) !whether particles are within computational boundary
    allocate(touch_bdry_i_mid(fixed_large_size)) !whether particles are within computational boundary

    allocate(vr_i_integer(fixed_large_size)) 
    allocate(vz_i_integer(fixed_large_size)) 
    allocate(vphi_i_integer(fixed_large_size)) 

    allocate(vr_i_integer_mid(fixed_large_size)) 
    allocate(vz_i_integer_mid(fixed_large_size)) 
    allocate(vphi_i_integer_mid(fixed_large_size)) 

    allocate(vr_i_mid(fixed_large_size)) 
    allocate(vz_i_mid(fixed_large_size)) 
    allocate(vphi_i_mid(fixed_large_size)) 

    allocate(vpar_i(fixed_large_size)) !velocity components in magnetic coordinates
    allocate(vx_i(fixed_large_size)) 
    allocate(vy_i(fixed_large_size)) 
    allocate(grad_psi_i(fixed_large_size)) 
    allocate(grad_alpha_i(fixed_large_size)) 
    allocate(grad_psi_dot_grad_alpha_i(fixed_large_size)) 
    allocate(bval_i(fixed_large_size)) 

    allocate(r_i_old(fixed_large_size)) 
    allocate(z_i_old(fixed_large_size)) 
    allocate(phi_i_old(fixed_large_size)) 

    allocate(vr_i_old(fixed_large_size)) 
    allocate(vz_i_old(fixed_large_size)) 
    allocate(vphi_i_old(fixed_large_size)) 

    allocate(v_i_mid(fixed_large_size))
    allocate(vpar_i_mid(fixed_large_size))
    allocate(vx_i_mid(fixed_large_size))
    allocate(vy_i_mid(fixed_large_size))

    allocate(grad_psi_i_mid(fixed_large_size))
    allocate(grad_alpha_i_mid(fixed_large_size))
    allocate( grad_psi_dot_grad_alpha_i_mid(fixed_large_size))
    allocate(bval_i_mid(fixed_large_size))

    allocate(my_den_i_left(mtor,nrad), source=zero)
    allocate(my_den_i_right(mtor,nrad), source=zero)
  
  end subroutine initialize_fk

end module fk_module



module sort_ions

contains
subroutine sort_ions_according_to_poloidal_location(theta)
  use constants,only:p_
  use constants,only: twopi
  use pputil
  use fk_module
  implicit none
  real(p_),intent(in):: theta(:)
  integer:: ierr,np_old,np_new
  !assign particles to the different processors according to their theta coordinates, using the subroutines provided in pputil_yj.f90

  np_old=nmarker_i
  call init_pmove(theta(:),np_old,twopi,ierr)

  call pmove(ps_vol_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(w_i_star(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(r_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

 call pmove(r_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phi_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(radcor_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(theta_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(radcor_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(theta_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(alpha_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(vr_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(vr_i_integer_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i_integer_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i_integer_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vr_i_integer(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i_integer(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i_integer(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(vr_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
!!$
  call pmove(v_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vx_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vy_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_alpha_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_dot_grad_alpha_i(:),    np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(bval_i(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove(v_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vpar_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vx_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vy_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_alpha_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(grad_psi_dot_grad_alpha_i_mid(:),    np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(bval_i_mid(:),     np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call pmove2(active_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(active_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(touch_bdry_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove2(touch_bdry_i_mid(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  nmarker_i=np_new

  !     call check_domain_particles(theta,nmarker_i)
end subroutine sort_ions_according_to_poloidal_location

end module sort_ions

subroutine check_domain_particles(theta,nmarker_i) !pass the test, comfirming domain decomposition is consistent with particles grouping
  use constants,only:p_
  use domain_decomposition,only:theta_start,dtheta2
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: theta(nmarker_i)
integer:: k

do k=1,nmarker_i
if(theta(k)<theta_start .or. theta(k)>theta_start+dtheta2) write(*,*) 'warningg*** particle not in domain'
enddo

end subroutine check_domain_particles


subroutine load_fk()
  !spatial location in magnetic coordinates (psi,theta,phi) and then transform to cylindrical coordinates
  use constants,only:p_, one,twopi,pi,two,kev,fourpi
 
  use magnetic_coordinates,only: radcor_low=>xlow,radcor_upp=>xupp,vol, &
       &  mpol,nrad,xgrid,zgrid,jacobian,toroidal_range,jacobian
  use misc, only: magnetic_coordinates_to_cylindrical_coordinates
  use fk_module,only: total_nmarker_i,mass_i,ti0, vn_fk
  use fk_module,only: nmarker_i,radcor_i,theta_i,active_i ! as output
  use fk_module,only:active_i_mid,touch_bdry_i_mid,touch_bdry_i
  use fk_module,only: alpha_i,tor_shift_i !only allocate these array, their values are not assigned in the present subroutine
  use fk_module,only: r_i,z_i,phi_i,v_i,vr_i,vz_i,vphi_i !as output, loading using magnetic coordinates, the Boris pusher works in cylindrical coordinates, therefore, we need to transform from mag. cor. to cylin. cor.
  use fk_module,only: r_i_old,z_i_old,phi_i_old
  use fk_module,only: ps_vol_i,normalizing_factor !as output
  use fk_module,only: w_i,w_i_mid,w_i_star !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vpar_i, vx_i,vy_i  !only allocate the array, the value is not set in this subroutine
  use fk_module,only: grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i !only allocate the array, the value is not set in this subroutine
  use fk_module,only: r_i_mid,z_i_mid,phi_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vr_i_old,vz_i_old,vphi_i_old
  use fk_module,only: vr_i_mid,vz_i_mid,vphi_i_mid !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vr_i_integer,vz_i_integer,vphi_i_integer !only allocate the array, the value is not set in this subroutine
  use fk_module,only: vr_i_integer_mid,vz_i_integer_mid,vphi_i_integer_mid
  use fk_module,only: vt_i,vmin_i,vmax_i !as output, vt_i in SI unit
  use fk_module, only: ion_spatial_loading_scheme, ion_velocity_loading_scheme
  use pputil !containing subroutines that sort particles into differnt processors
  use domain_decomposition,only: numprocs,myid
  use math, only :  random_yj
  use misc, only: abs_jacobian_func
  use interpolate_module

  implicit none

  integer:: iseed
  integer,parameter:: max_try=10000
  real(p_):: radcor_val,theta_val,rannum1,rannum2,rannum3,tmp

  integer:: i,ierr,j,file_unit
  !  integer:: status(MPI_STATUS_SIZE)
  character(5):: filename
  real(p_):: pos1,pos2,jacobian_val
  real(p_) :: abs_jacobian_max
  integer:: np_old,np_new
  real(p_):: vt,vmin,vmax,v_val,maxwellian_func_ion
  !real(p_),allocatable:: theta_v(:),phi_v(:)
  !real(p_),allocatable::vx(:),vy(:),vz(:),tmp_array(:)
  real(p_)::vx(nmarker_i),vy(nmarker_i),vz(nmarker_i),tmp_array(3*nmarker_i)
  real(p_)::maxwellian_min,maxwellian_max

  !  allocate(pitch_angle_i(fixed_large_size))
  !  allocate(gyro_angle_i(fixed_large_size))

  !  allocate(theta_v(nmarker_i)) !local array
  !  allocate(phi_v(nmarker_i)) !local array
!!$  allocate(vx(nmarker_i)) !local array
!!$  allocate(vy(nmarker_i)) !local array
!!$  allocate(vz(nmarker_i)) !local array
!!$  allocate(tmp_array(3*nmarker_i))

abs_jacobian_max=maxval(abs(jacobian(:,1:nrad)))
  !  radcor_min=minval(xgrid)
  !  radcor_max=maxval(xgrid)


  ! ---random generator, when use MPI_send to generate iseed for other processes, it is actual a sequence generator,instead of parallel generator
!!$  if ( myid .eq. 0 ) then ! master generates random numbers first, others wait in line
!!$     iseed = 0
!!$  else 
!!$     call MPI_Recv(iseed, 1, MPI_INT, myid-1, 1, MPI_COMM_WORLD, status,ierr) !other processes wait to receive the iseed
!!$  endif

  iseed=-(1777+myid*3) !set the iseed in different procs, when using this, it is a parallel generator, but the random numbers in different procs may be related if the iseed chosen for differnt procs is not good enough
  !  write(*,*) 'myid=',myid, 'iseed=',iseed

  ! now generate the random numbers
  tmp = random_yj(iseed) !just to trigger the use of the iseed, the generated random number is not used, 

  if(ion_spatial_loading_scheme.eq.1) then
     do i=1,nmarker_i
        rannum1 = random_yj(0) !0 means using last random number as iseed 
        rannum2 = random_yj(0) 
        radcor_val=radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
        !theta_val=rannum2*twopi  !theta in [0:2*pi]
        theta_val=(rannum2-0.5_p_)*twopi !theta in [-pi:+pi]
        radcor_i(i)=radcor_val
        theta_i(i)=theta_val
     enddo
  elseif (ion_spatial_loading_scheme.eq.2) then

     do i=1,nmarker_i     
        do j=1,max_try !rejection method to generate nonuniform random numbers
           rannum1 = random_yj(0) !0 means using last random number as iseed 
           rannum2 = random_yj(0) !use last random number as iseed
           !rannum3 = random_yj(0) !use last random number as iseed
           radcor_val=radcor_low+(radcor_upp-radcor_low)*rannum1 !scale the random number to the range [radcor_low: radcor_upp]
           theta_val=-pi+rannum2*twopi
           pos1=abs_jacobian_func(theta_val,radcor_val)
           !       write(*,*) 'abs_jacobian_func(theta_val,radcor_val)=',pos1
           pos2 = random_yj(0) !use last random number as iseed   
           pos2=pos2*abs_jacobian_max !scaled to the range [0: abs_jacobian_max]
           if(pos1<pos2) then
              cycle
           else
              radcor_i(i)=radcor_val
              theta_i(i)=theta_val
              !phi_i(i)=toroidal_range*rannum3
              exit
           endif
        enddo
        !     if(j.eq.max_try+1) stop "***stop**, rejection method is not successful in generating distribution"
        !     write(*,*) 'j=',j
     enddo
  else
     stop 'please specify a loading scheme for the spatial distribution ion markers'
  endif

  do i=1,nmarker_i
     call magnetic_coordinates_to_cylindrical_coordinates(theta_i(i),radcor_i(i),r_i(i),z_i(i)) !to get the corresponding (R,Z) coordinates
  enddo

  do i=1,nmarker_i !setting toroidal coordinate of particles
     rannum3 = random_yj(0) !use last random number as iseed
     phi_i(i)=toroidal_range*rannum3
  enddo

  !setting velocity
  vt=sqrt(two*ti0*kev/mass_i)
  vmin=-3._p_*vt/vn_fk !normalized by vn_fk
  vmax=+3._p_*vt/vn_fk !normalized by vn_fk
  maxwellian_max=maxwellian_func_ion(0._p_)

!if(myid.eq.0) write(*,*)  'vt=',vt, 'vmin=',vmin,'vmax=',vmax, 'vmin*vn_fk=',vmin*vn_fk, 'vt/vn_fk=',vt/vn_fk

  if(ion_velocity_loading_scheme.eq.1) then !using uniform loading in v, instead of Gaussian
     do i=1,3*nmarker_i        
        rannum1 = random_yj(0) !0 means using last random number as iseed 
        v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
        tmp_array(i)=v_val
     enddo
  elseif (ion_velocity_loading_scheme.eq.2) then
     do i=1,3*nmarker_i        
        do j=1,max_try !rejection method to generate nonuniform random numbers
           rannum1 = random_yj(0) !0 means using last random number as iseed 
           v_val=vmin+rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
           pos1=maxwellian_func_ion(v_val*vn_fk)
           pos2 = random_yj(0) !0 means using last random number as iseed   
           pos2=pos2*maxwellian_max !scaled to [0,maxwellian_max]
           if(pos1<pos2) then
              cycle
           else
              tmp_array(i)=v_val
              exit
           endif
        enddo
!        if(myid.eq.0) write(*,*) 'j=',j
     enddo
  else
     stop 'please specify a loading scheme for the velocity distribution ion markers'
  endif

  do i=1,nmarker_i
     vx(i)=tmp_array(i)
     vy(i)=tmp_array(i+nmarker_i)
     vz(i)=tmp_array(i+2*nmarker_i)
  v_i(i)=sqrt(vx(i)*vx(i)+vy(i)*vy(i)+vz(i)*vz(i))
  enddo

!if(myid.eq.0) call calculate_possibility_density(vz,nmarker_i,100,vmin,vmax)

  do i=1,nmarker_i
     vz_i(i)=vz(i)
     vr_i(i)=vx(i)*cos(phi_i(i))+vy(i)*sin(phi_i(i))
     vphi_i(i)=-vx(i)*sin(phi_i(i))+vy(i)*cos(phi_i(i))
  enddo

!  v_i=sqrt(vr_i*vr_i+vz_i*vz_i+vphi_i*vphi_i)
  !if(myid.eq.3) call calculate_possibility_density(v_i,nmarker_i,100,vmin,vmax)

  !  v_i=v_i/vn_fk !normalized by vn_fk
!!$  do i=1,nmarker_i !setting direction of velocity
!!$     rannum1 = random_yj(0) !0 means using last random number as iseed
!!$     theta_v(i)=pi*rannum1
!!$     rannum1 = random_yj(0) !0 means using last random number as iseed
!!$     phi_v(i)=twopi*rannum1
!!$  enddo
  !transform to components in cylindrical coordinates
!!$  do i=1,nmarker_i !velocity components in a constant Cartesian coordinate system
!!$     vz_i(i)=v_i(i)*cos(theta_v(i))
!!$     vx(i)=v_i(i)*sin(theta_v(i))*cos(phi_v(i))
!!$     vy(i)=v_i(i)*sin(theta_v(i))*sin(phi_v(i))
!!$  enddo

!!$  do i=1,nmarker_i !projected onto the basis vectors of cylindrical coordinates
!!$     vr_i(i)=vx(i)*cos(phi_i(i))+vy(i)*sin(phi_i(i))
!!$     vphi_i(i)=vy(i)*cos(phi_i(i))-vx(i)*sin(phi_i(i))
!!$  enddo

  if(ion_spatial_loading_scheme.eq.1 .and. ion_velocity_loading_scheme.eq.1) then
     !normalizing_factor=total_nmarker_i/(twopi*toroidal_range*(radcor_upp-radcor_low)*twopi*pi*(vmax-vmin))
     normalizing_factor=total_nmarker_i/(twopi*toroidal_range*(radcor_upp-radcor_low)*(vmax-vmin)**3)
     do i=1,nmarker_i
        call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor)
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor)
     enddo
  elseif(ion_spatial_loading_scheme.eq.1 .and. ion_velocity_loading_scheme.eq.2) then
     !normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*twopi*pi*vt/vn_fk*sqrt(pi)/two)
     normalizing_factor=total_nmarker_i/((radcor_upp-radcor_low)*twopi*toroidal_range*(sqrt(pi)*vt/vn_fk)**3)
     do i=1,nmarker_i
        call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,theta_i(i),radcor_i(i),jacobian_val)
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))*abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk))
        ps_vol_i(i)=abs(jacobian_val)/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk))
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.1) then
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*(vmax-vmin))
     normalizing_factor=total_nmarker_i/(vol*(vmax-vmin)**3)
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor)
        ps_vol_i(i)=one/(normalizing_factor)
     enddo

  elseif(ion_spatial_loading_scheme.eq.2 .and. ion_velocity_loading_scheme.eq.2) then
     !  normalizing_factor=total_nmarker_i/(vol*fourpi*vt/vn_fk*sqrt(pi)/two) !wrong
     !normalizing_factor=total_nmarker_i/(vol*twopi*pi*vt/vn_fk*sqrt(pi)/two) !wrong again
     normalizing_factor=total_nmarker_i/(vol*(sqrt(pi)*vt/vn_fk)**3) !corrected
     do i=1,nmarker_i
        !ps_vol_i(i)=v_i(i)**2*sin(theta_v(i))/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk)) !wrong
        ps_vol_i(i)=one/(normalizing_factor*maxwellian_func_ion(v_i(i)*vn_fk))
     enddo
  endif


!!$   iseed=next_seed
!!$  if (myid .ne. numprocs-1) then
!!$     call MPI_Send(iseed, 1, MPI_INT, myid+1, 1, MPI_COMM_WORLD,ierr)  !send the iseed to next process
!!$  endif

  !assign the loaded particles to the corresponding processors, using the subroutines provided in pputil_yj.f90
  np_old=nmarker_i
  call init_pmove(theta_i(:),np_old,twopi,ierr)
  call pmove(theta_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(radcor_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(r_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(z_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(phi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(v_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vr_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vz_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(vphi_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit
  call pmove(ps_vol_i(:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  nmarker_i=np_new

  !  call some_test3(myid,numprocs)

!!$  write(filename,'(a1,i4.4)') 'i',myid
!!$  open(newunit=file_unit, file=filename)
!!$  do i=1,nmarker_i
!!$     write(file_unit,'(4(1pe14.5),i6,i3)')  radcor_i(i) ,theta_i(i),r_i(i),z_i(i),i,myid
!!$     !     write(file_unit,*)  vr_i(i),vphi_i(i),vz_i(i)
!!$  enddo
!!$  close(file_unit)

  active_i=.true. ! initially, all markers are active, i.e., within the computational region
  touch_bdry_i=.false.
  active_i_mid=.true. ! initially, all markers are active, i.e., within the computational region
  touch_bdry_i_mid=.false.

  vmin_i=vmin
  vmax_i=vmax
  vt_i=vt
end subroutine load_fk

function maxwellian_func_ion(v) result(z)
use constants,only:p_
  use constants,only:two,kev
  use fk_module,only: mass_i,ti0 !as input
implicit none
real(p_):: v,z

z=exp(-mass_i*v*v/(two*ti0*kev)) !the normalizing factor (mi/(twopi*Ti*kev))^(3/2) is not included


end function maxwellian_func_ion

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
    character(100) :: density_file, density_radcor, temperature_file, temperature_radcor
    real(p_) :: density_unit, temperature_unit
    namelist/adiabatic_electron_nmlt/density_file, density_unit, density_radcor, &
         & temperature_file, temperature_unit, temperature_radcor
    integer :: u
    
    open(newunit=u,file='input.nmlt')
    read(u, adiabatic_electron_nmlt)
    close(u)

    call ne_object%set_profile(density_file, density_unit, 1.0d0, density_radcor) 
    call te_object%set_profile(temperature_file, temperature_unit, kev, temperature_radcor) 

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
    !call alpha_ecrit%set_profile("cfetr/alpha_ecrit.txt", Mev, kev, "poloidal-flux-sqrt")
    !call alpha_normc%set_profile("cfetr/alpha_normc.txt", one, one, "poloidal-flux-sqrt")

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
module gk_module
  use constants, only: p_, twopi, elementary_charge, kev
  use normalizing, only : ln, bn
  implicit none
  save
  integer :: nsm !total number of species, get its value from the input namelist
  integer :: nmmax !maximal number of markers in a processor
  integer, allocatable ::  nm_gk(:) !actual number of markers in a processor, fluctuates with time
  integer, allocatable :: total_nm_gk(:)  !total number of markers (in all the processors) for a species
  integer, allocatable :: gk_spatial_loading_scheme(:), gk_velocity_loading_scheme(:)
  integer, allocatable :: gk_nonlinear(:)
  real(p_), allocatable :: mass_gk(:), charge_gk(:), charge_sign_gk(:)
  real(p_), allocatable :: tgk0(:) !in kev, a typical value of temperature in the compuational region
  real(p_), allocatable :: ngk0(:) !in m^(-3) a typical value of density in the compuational region
  real(p_), allocatable :: vn_gk(:), dtao_gk(:)
  logical, allocatable :: gk_flr(:)
  integer, parameter ::  gyro_npt = 4 !number of points on a gyro-ring
  real(p_) :: w_unit
  real(p_), allocatable :: xgc(:,:), zgc(:,:), ygc(:,:) !magnetic coordinates of guiding centers
  real(p_), allocatable :: vpar_gk(:,:), w_gk(:,:) !weight of gk markers

  real(p_), allocatable :: xgc_mid(:,:), zgc_mid(:,:), ygc_mid(:,:)
  real(p_), allocatable :: vpar_gk_mid(:,:), w_gk_mid(:,:) !vlaues at half-time-step
  
  real(p_), allocatable :: mu_gk(:,:), v_gk(:,:)
  real(p_), allocatable :: ps_vol_gk(:,:)
  logical, allocatable :: lost_gc(:,:) !indicates whether a guiding-center location exceed the radial boundary
  !computed on fly, i.e., temporary working arrays. Hence, no time subfix (_mid), and no need to sort
  real(p_), allocatable :: x_ring(:,:,:), y_ring(:,:,:), z_ring(:,:,:)   
  !---for testing
  !integer :: ntouch_bdry_gc=0, total_ntouch_bdry_gc !numbe of markers that touch the boundary in each process and all processes

contains
  subroutine initialize_gk(nsm, baxis, minor_a, dt_omega_i_axis)
    use domain_decomposition,only : numprocs, myid
    use magnetic_coordinates, only : nrad, mpol2, mtor, pfn, pfn_inner, pfn_bdry
    use gk_radial_profiles, only : initialize_gk_radial_profiles
    use gk_profile_funcs, only: gkt_func, gkn_func, gkdtdx_func, gkdndx_func
    integer, intent(in) :: nsm
    real(p_), intent(in) :: baxis, minor_a, dt_omega_i_axis
    real(p_), allocatable ::  omega_gk_axis(:), rho_gk(:)
    character(100) :: density_file(nsm), density_radcor(nsm), temperature_file(nsm), temperature_radcor(nsm)
    real(p_) :: density_unit(nsm), temperature_unit(nsm), pfn0
    character(len=100) :: format1
    integer, allocatable  :: nm_gk_per_cell(:)
    namelist/gk_nmlt/ mass_gk, charge_gk, gk_flr, nm_gk_per_cell, density_file, density_unit, &
         & density_radcor, temperature_file, temperature_unit, temperature_radcor, &
         & gk_spatial_loading_scheme, gk_velocity_loading_scheme, gk_nonlinear
    integer :: u, i, k

    allocate(mass_gk(nsm), charge_gk(nsm), charge_sign_gk(nsm))
    allocate(ngk0(nsm), tgk0(nsm))
    allocate(omega_gk_axis(nsm), dtao_gk(nsm), vn_gk(nsm))
    allocate(rho_gk(nsm))
    allocate(nm_gk_per_cell(nsm))
    allocate(total_nm_gk(nsm))
    allocate(nm_gk(nsm))
    allocate(gk_flr(nsm))
    allocate(gk_nonlinear(nsm))
    allocate(gk_spatial_loading_scheme(nsm))
    allocate(gk_velocity_loading_scheme(nsm))

    open(newunit=u,file='input.nmlt')
    read(u,gk_nmlt)
    close(u)
    if(myid==0)  write(*,gk_nmlt)
    charge_gk = charge_gk*elementary_charge
    charge_sign_gk(:) = sign(1._p_, charge_gk(:))
    vn_gk(:) = ln/(twopi/(bn*abs(charge_gk(:))/mass_gk(:)))
    omega_gk_axis=abs(baxis*charge_gk)/mass_gk

    total_nm_gk(:) = nm_gk_per_cell(:) * nrad * mpol2 * mtor
    nm_gk(:) = total_nm_gk(:)/numprocs !the number of markers initially loaded per processor
    !Later, nm_gk will be set to the actual number of markers per proc 
    !the value of which will be differnt for differnt processors and at differnt time

    !marker number in one process fluctuates with time, so choose a large number (fixed) to be safe
    nmmax = int((maxval(total_nm_gk)/numprocs)*2.5) 

    allocate(xgc(nmmax, nsm)) 
    allocate(zgc(nmmax,nsm))
    allocate(ygc(nmmax,nsm))

    allocate(v_gk(nmmax, nsm))
    allocate(vpar_gk(nmmax, nsm))
    allocate(mu_gk(nmmax, nsm))
    allocate(w_gk(nmmax, nsm)) 
    allocate(w_gk_mid(nmmax, nsm)) 
    allocate(ps_vol_gk(nmmax, nsm))

    allocate(xgc_mid(nmmax,nsm)) 
    allocate(zgc_mid(nmmax,nsm))
    allocate(ygc_mid(nmmax,nsm))
    allocate(vpar_gk_mid(nmmax,nsm))
    allocate(x_ring(gyro_npt, nmmax, nsm))
    allocate(y_ring(gyro_npt, nmmax, nsm))
    allocate(z_ring(gyro_npt, nmmax, nsm))
    allocate(lost_gc(nmmax,nsm))
    if(myid.eq.0) write(*,*) 'marker number of each gk species=',total_nm_gk

    call initialize_gk_radial_profiles(nsm, density_file, density_unit, density_radcor, &
         &   temperature_file, temperature_unit, temperature_radcor)



    do i = 1, nsm
       pfn0 =  0.0
       tgk0(i) = gkt_func(pfn0, i) !typical temperature for a species
       ngk0(i) = gkn_func(pfn0, i) !typical number density for a species
    enddo
    format1='(A20,20ES18.4)'
    rho_gk = sqrt(tgk0*kev/mass_gk)/omega_gk_axis
    if(myid==0) write(*, format1) 'rho_gk(:) in meter=', rho_gk
    if(myid==0) write(*, format1) 'a/rho_gk(:)=', minor_a/rho_gk
    !time step in unit of the gyro-period twopi/abs(omegan_gk(:)):
    dtao_gk(:) = dt_omega_i_axis*(bn*abs(charge_gk(:))/mass_gk(:))/omega_gk_axis(1)/twopi
    if(myid==0) write(*, format1) 'dt/tn_gk=',dtao_gk
    if(myid==0) write(*, format1) 'tgk0=', tgk0

    block !diagnostic information
      use magnetic_field, only: tfn_func_pfn
      use func_in_mc, only: minor_r_radcor, minor_r_prime
      character(len=100) :: fn
      integer :: u, ns
      real(p_) :: tfn, c

      if(myid==0) then
         do ns = 1, nsm
            fn = 'profiles_nsx.txt'
            write(fn(12:12),'(i1.1)') ns
            open(newunit=u, file=fn)
            do i= 1, size(pfn)
               tfn = tfn_func_pfn(pfn(i))
               c = minor_r_prime(pfn(i))
               write(u,'(80ES18.6E4)') pfn(i), tfn, minor_r_radcor(pfn(i)), &
                    & gkn_func(pfn(i), ns), gkt_func(pfn(i), ns), &
                    & gkdndx_func(pfn(i), ns)/c,   gkdtdx_func(pfn(i), ns)/c
            enddo
            close(u)
         enddo
      endif
    endblock

  end subroutine initialize_gk

  subroutine sort_gk_markers(ns, theta, step)
    !assign particles to different processors according to their poloidal coordinates theta
    use constants, only: twopi
    use pputil, only: init_pmove, pmove, ppexit, pmove2
    implicit none
    integer, intent(in) :: ns
    real(p_), intent(in) :: theta(:)
    integer, intent(in) :: step
    integer :: np_old, np_new, ierr

    np_old = nm_gk(ns)

    call init_pmove(theta(:), np_old, twopi, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(xgc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(zgc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ygc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(vpar_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(mu_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(v_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove2(lost_gc(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(ps_vol_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit
    call pmove(w_gk(:,ns), np_old, np_new, ierr)
    if (ierr.ne.0) call ppexit

    if (step==1) then !need sorting only at the half time-step
       call pmove(xgc_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(zgc_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(ygc_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(vpar_gk_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
       call pmove(w_gk_mid(:,ns), np_old, np_new, ierr)
       if (ierr.ne.0) call ppexit
    endif

    nm_gk(ns) = np_new
  end subroutine sort_gk_markers

end module gk_module
module load_gk_mod
contains
  subroutine load_gk(ns, nm, xgc, zgc, ygc, vpar_gk,mu_gk,v_gk, ps_vol_gk, w_gk)
    use constants,only: p_, one,twopi,pi,two,kev,fourpi, zero
    use magnetic_coordinates,only: xlow, xupp, vol, &
         & mpol,nrad, zgrid, xgrid, tor_shift_mc, &
         & jacobian, toroidal_range, toroidal_range
    use gk_module,only: total_nm_gk, mass_gk, tgk0, w_unit, &
         & vn_gk, gk_spatial_loading_scheme, gk_velocity_loading_scheme 
    use magnetic_field, only : pfn_func
    use gk_profile_funcs,only : gkn_func
    use domain_decomposition,only: numprocs, myid
    use misc, only: magnetic_coordinates_to_cylindrical_coordinates
    use interpolate_module
    use math, only: random_yj
    use misc, only: abs_jacobian_func
    use pputil
    use math,only: shift_toroidal
    use func_in_mc, only:  br_mc_func, bz_mc_func, bphi_mc_func
    use map_to_mc, only: interpolate_from_cylindrical_to_magnetic_coordinates
    implicit none
    integer, intent(in) :: ns
    integer, intent(inout) :: nm
    real(p_), intent(out) :: xgc(:), zgc(:), ygc(:)
    real(p_), intent(out) :: vpar_gk(:), mu_gk(:), v_gk(:), ps_vol_gk(:), w_gk(:)

    integer :: iseed
    integer, parameter :: max_try=10000 !used in rejection method
    real(p_) :: radcor_val,theta_val,rannum1,rannum2, rannum3, z1, z2, tmp, scalar
    integer :: i, j, ierr, nmax
    real(p_) :: pos1,pos2, pos_max, w0
    real(p_) :: abs_jacobian_max, jacobian_val
    integer :: np_old, np_new
    real(p_) :: vt, vmin, vmax, v_val, probability_max, maxwellian_max
    real(p_) :: theta_v(nm), phi_v(nm), theta_v0
    real(p_) :: vx(nm), vy(nm), vz(nm), tmp_array(3*nm), vx0,vy0,vz0, v0
    real(p_) :: rp_e(nm), zp_e(nm), phip_e(nm) !particle locations
    real(p_) :: brval, bphival, vr, vphi, angle, bval,bxval,byval,bzval
    real(p_) :: tor_shift_e,   normalizing_factor, t(mpol, 1:nrad), tmax
    real(p_) :: rg_e(nm), zg_e(nm), phig_e(nm), gaussian_envelope

    abs_jacobian_max = maxval(abs(jacobian(:,1:nrad)))

    do i =1, mpol
       do j = 1, nrad
          t(i,j) = abs(jacobian(i,j))*gkn_func(xgrid(j), ns)
       enddo
    enddo
    tmax = maxval(t)

    !set different initial seeds for different myid or ns, so that the random number series are not identical
    iseed = 27777 + myid*10000 + ns*1111
    tmp = random_yj(iseed) !to trigger the use of the iseed, the generated random number is not used, 

    select case(gk_spatial_loading_scheme(ns))
    case(1) ! Uniform loading in (radcor,theta). My testing indicates this loading give more accurate Monte-Carlo integration than uniform loading in real space
       do i = 1, nm     
          rannum1 = random_yj(0) !0 means using last random number as iseed 
          rannum2 = random_yj(0) 
          xgc(i)  = xlow + (xupp-xlow)*rannum1 !scale the random number to the range [xlow: xupp]
          zgc(i) = (rannum2 - 0.5_p_)*twopi
       enddo

    case(2) !uniform in real space
       do i=1,nm
          do j=1,max_try !rejection method to generate nonuniform random numbers
             rannum1 = random_yj(0) 
             rannum2 = random_yj(0) 
             radcor_val = xlow +(xupp-xlow)*rannum1 !scaled to the range [xlow: xupp]
             theta_val = (rannum2-0.5_p_)*twopi !scaled to the range [-pi:pi]
             pos1 = abs_jacobian_func(theta_val, radcor_val)
             pos2 = random_yj(0)
             pos2=pos2*abs_jacobian_max !scaled to the range [0: abs_jacobian_max]
             if(pos1<pos2) then
                cycle
             else
                xgc(i) = radcor_val
                zgc(i) = theta_val
                exit
             endif
          enddo
          if(j.eq.max_try+1) stop "***stop**, rejection method falied"
       enddo
    case(3) !nonuniform in real space,  proportional to the equilibrium density n0(x)
       do i=1,nm
          do j=1,max_try 
             rannum1 = random_yj(0)
             rannum2 = random_yj(0) 
             radcor_val = xlow +(xupp-xlow)*rannum1 
             theta_val = (rannum2-0.5_p_)*twopi 
             pos1 = abs_jacobian_func(theta_val, radcor_val)*gkn_func(radcor_val, ns)
             pos2 = random_yj(0)
             pos2=pos2*tmax
             if(pos1<pos2) then
                cycle
             else
                xgc(i) = radcor_val
                zgc(i) = theta_val
                exit
             endif
          enddo
          if(j.eq.max_try+1) stop "***stop**, rejection method falied"
       enddo
    case default
       stop 'please specify a loading scheme for the spatial distribution of gk markers'
    end select

    do i=1,nm !set toroidal coordinates of markers
       rannum1=random_yj(0)
       ygc(i) = toroidal_range*rannum1
    enddo

    do i=1,nm
       call linear_2d_interpolate(mpol, nrad, zgrid, xgrid,tor_shift_mc,zgc(i),xgc(i),tor_shift_e) 
       phip_e(i) = ygc(i) + tor_shift_e
    enddo

    vt = sqrt(two*tgk0(ns)*kev/mass_gk(ns)) 

    select case(gk_velocity_loading_scheme(ns))
    case(1) !uniform in velocity using Cartesian coordinates (vx, vy, vz)
       vmin = -3_p_*vt/vn_gk(ns) 
       vmax = +3_p_*vt/vn_gk(ns)
       do i = 1, 3*nm
          rannum1=random_yj(0) 
          tmp_array(i) = vmin + rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
       enddo

       do i = 1, nm
          vx(i) = tmp_array(i)
          vy(i) = tmp_array(i+nm)
          vz(i) = tmp_array(i+2*nm)
       enddo
       v_gk(1:nm) = sqrt(vx(:)**2 + vy(:)**2 + vz(:)**2)

    case(2) !maxwellian in veloicty using spherical coordinates
       probability_max = probability_sphere(vt, ns)
       vmax= 5.2_p_*vt/vn_gk(ns) !I found the range is important in enusuring the code does not blowup
       do i =1, nm !set velocity magnitude
          do j = 1, max_try !rejection method
             rannum1 = random_yj(0) 
             v0 = rannum1*vmax
             pos1 = probability_sphere(v0*vn_gk(ns), ns)
             rannum2 = random_yj(0)
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
             rannum1 = random_yj(0) 
             theta_v0 = rannum1*pi
             pos1 = sin(theta_v0)
             rannum2 = random_yj(0)
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
          rannum1 = random_yj(0)
          phi_v(i) = twopi*rannum1
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
!!$             rannum1 = random_yj(0)
!!$             v_val = vmin + rannum1*(vmax-vmin) !scale the random number to [vmin:vmax]
!!$             pos1 = maxwellian0(v_val*vn_gk(ns), ns)
!!$             pos2=random_yj(0) 
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
          rannum1 = random_yj(0)
          rannum2 = random_yj(0)
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

    do i = 1, nm
       call magnetic_coordinates_to_cylindrical_coordinates(zgc(i), xgc(i), rp_e(i), zp_e(i))
    enddo

    do i = 1, nm
       brval = br_mc_func(zgc(i), xgc(i))
       bzval = bz_mc_func(zgc(i), xgc(i))
       bphival = bphi_mc_func(zgc(i), xgc(i))
       angle = atan2(vy(i),vx(i))
       vr = vx(i)*cos(angle)+ vy(i)*sin(angle)
       vphi = -vx(i)*sin(angle)+vy(i)*cos(angle)
       call particle_to_guiding_center_location(ns, rp_e(i), phip_e(i), zp_e(i), vr,vphi,vz(i),brval,bphival,bzval, &
            & rg_e(i),phig_e(i), zg_e(i))

       xgc(i) = pfn_func(rg_e(i), zg_e(i)) !calculate the radial coordinate
       call interpolate_from_cylindrical_to_magnetic_coordinates(rg_e(i), zg_e(i), zgc(i), tor_shift_e)
       !ygc(i) = phig_e(i) + tor_shift_e ! a bug?
       ygc(i) = phig_e(i) - tor_shift_e
    enddo


    !if(myid.eq.0) call calculate_possibility_density(vpar_gk,nm,100,minval(vpar_gk),maxval(vpar_gk))
    !if(myid.eq.0) call calculate_possibility_density(mu_gk,nm,100,minval(mu_gk),maxval(mu_gk))

    if(gk_spatial_loading_scheme(ns)==1 .and. gk_velocity_loading_scheme(ns)==1) then
       do i = 1, nm
          call linear_2d_interpolate(mpol, nrad, zgrid, xgrid, jacobian, zgc(i), xgc(i), jacobian_val)
          ps_vol_gk(i) = abs(jacobian_val)*twopi*toroidal_range*(xupp-xlow)*(vmax-vmin)**3/total_nm_gk(ns)
       enddo

    elseif(gk_spatial_loading_scheme(ns)==1 .and. (gk_velocity_loading_scheme(ns)==2 .or. gk_velocity_loading_scheme(ns)==3)) then
       normalizing_factor=total_nm_gk(ns)/(twopi*toroidal_range*(xupp-xlow)*(sqrt(pi)*vt/vn_gk(ns))**3)
       do i = 1, nm
          call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,jacobian,zgc(i),xgc(i),jacobian_val)
          ps_vol_gk(i) = abs(jacobian_val)/(normalizing_factor*maxwellian0(v_gk(i)*vn_gk(ns),ns))
       enddo
    elseif(gk_spatial_loading_scheme(ns)==2 .and. gk_velocity_loading_scheme(ns)==1) then
       normalizing_factor = total_nm_gk(ns)/(vol*(vmax-vmin)**3)
       do i=1,nm
          ps_vol_gk(i)=one/(normalizing_factor)
       enddo
    elseif((gk_spatial_loading_scheme(ns)==2) .and. (gk_velocity_loading_scheme(ns)==2 .or. gk_velocity_loading_scheme(ns)==3)) then
       normalizing_factor = total_nm_gk(ns)/(vol*(sqrt(pi)*vt/vn_gk(ns))**3) 
       do i = 1, nm
          ps_vol_gk(i) = one/(normalizing_factor*maxwellian0(v_gk(i)*vn_gk(ns),ns))
       enddo
    endif

    do i = 1, nm !set initial perturbation
       w0 = ps_vol_gk(i) * f0(xgc(i), v_gk(i), ns) / w_unit
       gaussian_envelope = exp(-((xgc(i)-xgrid(nrad/2))/((xgrid(1)-xgrid(nrad))/4))**2)
       w_gk(i) = w0*(random_yj(0)-0.5)*1d-4*gaussian_envelope
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

  pure real(p_) function f0(x,v,ns) result (z) ! v in unit of vn, f0 in unit 1/(Ln**3*vn**3)
    use constants, only: p_, two, twopi, kev
    use normalizing, only: Ln
    use gk_module, only: mass_gk, vn_gk
    use gk_profile_funcs, only: gkt_func, gkn_func
    implicit none
    integer, intent(in) :: ns
    real(p_), intent(in) :: x, v
    real(p_) :: v_si, te

    v_si = v*vn_gk(ns)
    te = gkt_func(x,ns)*kev
    z = gkn_func(x,ns)*sqrt((mass_gk(ns)/(twopi*te))**3)*exp(-mass_gk(ns)*v_si**2/(two*te))
    z = z*(vn_gk(ns)**3*Ln**3)
  end function f0

  pure real(p_)  function maxwellian0(v, ns) result(z) !v in SI unit
    use constants, only: p_, two, kev
    use gk_module, only: mass_gk, tgk0 !as input
    implicit none
    real(p_), intent(in) :: v
    integer, intent(in) :: ns
    z = exp(-mass_gk(ns)*v*v/(two*tgk0(ns)*kev))
  end function maxwellian0

  pure real(p_)  function probability_sphere(v, ns) result(z) !v in SI unit
    use constants, only: p_
    implicit none
    real(p_), intent(in) :: v
    integer, intent(in) :: ns
    z = v**2*maxwellian0(v,ns)
  end function probability_sphere


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
subroutine field_lines_analyse()
  use constants,only:p_
  use constants,only: two
  use radial_module,only: z_axis
  implicit none
  
  integer:: n_tor_loop,max_npt_along_field_line,krad !used in field line tracing module
!  namelist /field_line_tracing_nl/  n_tor_loop,max_npt_along_field_line,krad
  real(p_),allocatable:: r_start(:),z_start(:),phi_start(:) !starting point of field lines
  real(p_),allocatable:: r_poincare(:,:),z_poincare(:,:),phi_poincare(:,:) !pointcare points
  integer:: j,k
  integer,allocatable::nloop_actual(:)

n_tor_loop=10
max_npt_along_field_line=8000
krad=1

  
!!$  open(31,file='input.nmlt')
!!$  read(31,field_line_tracing_nl)
!!$  close(31)
!!$  write(*,field_line_tracing_nl)

  allocate(r_start(krad))
  allocate(z_start(krad))
  allocate(phi_start(krad))
  allocate(r_poincare(n_tor_loop+1,krad))
  allocate(z_poincare(n_tor_loop+1,krad))
  allocate(phi_poincare(n_tor_loop+1,krad))
  allocate(nloop_actual(krad))

!  !$omp parallel do
  do k=1,krad
     r_start(k)=2.15_p_+(k-1)*0.15d0/(20-1)
     z_start(k)=z_axis
     phi_start(k)=0._p_
     call field_line_tracing(r_start(k),z_start(k),phi_start(k), max_npt_along_field_line,n_tor_loop,&
          & r_Poincare(:,k),z_Poincare(:,k),phi_Poincare(:,k),nloop_actual(k))
  enddo
!  !$omp end parallel do

  open(26,file='poincare.txt')
  do k=1,krad
     do j=1,nloop_actual(k)
        write(26,*) r_Poincare(j,k),z_Poincare(j,k),phi_Poincare(j,k)
     enddo
     write(26,*)
     write(26,*)
  enddo
  close(26)

!!$  do k=1,krad
!!$     call draw_magnetic_surface(r_start(k),z_start(k),'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
!!$  enddo


end subroutine field_lines_analyse



subroutine field_line_tracing(r0,z0,phi0,npt,n_tor_loop,r_poincare,z_poincare,phi_poincare,nloop_actual)
  !given coordinates (R,Z,phi), this subroutine follows the field lines passing through this point until it has finish n_tor_loop toroidal loop or exceeds the specifed maximum number of points along the field-line, npt. This subroutine also calculates the safety factor of the field-line found.
  use constants,only:p_,two,twopi,one_half
  use boundary,only: np_lcfs,x_lcfs,z_lcfs,nlim,rlim,zlim !use to check whether field line touch the boundary
 use magnetic_field, only : br,bz,bphi
 implicit none

  real(p_),intent(in):: r0,z0,phi0
  integer,intent(in)::npt,n_tor_loop
  real(p_):: r(npt),z(npt),phi(npt)
  real(p_),intent(out):: r_poincare(n_tor_loop+1),z_poincare(n_tor_loop+1),phi_poincare(n_tor_loop+1)
  integer,intent(out):: nloop_actual
  real(p_),parameter:: step=1d-3 !meter, trial of dr or dz step
  real(p_):: brval,bzval,bphival,bpolval,dr,dz,dphi
  
  real(p_):: r_mid,z_mid,dl_pol,qval
  logical:: loss
  integer:: j,k,jj

  k=1 !Poincare points
  r_poincare(k)=r0 !Poincare points
  z_poincare(k)=z0
  phi_poincare(k)=phi0

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  loss=.false.
  do j=1,npt-1
     !2nd Runge-Kutta
     brval=    br(r(j),z(j))
     bzval=    bz(r(j),z(j))
     bphival=bphi(r(j),z(j))
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr=step*one_half
        if(brval.lt.0._p_) dr=-step*one_half
        dz=bzval/brval*dr
     else
        dz=step*one_half
        if(bzval.lt.0._p_) dz=-step*one_half
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/bpolval*dl_pol/r(j)

     !first step:
     r_mid=r(j)+dr
     z_mid=z(j)+dz
     brval=    br(r_mid,z_mid)
     bzval=    bz(r_mid,z_mid)
     bphival=bphi(r_mid,z_mid)
     bpolval=sqrt(brval**2+bzval**2)


     if(abs(bzval).lt.abs(brval)) then
        dr=step
        if(brval.lt.0._p_) dr=-step
        dz=bzval/brval*dr
     else
        dz=step
        if(bzval.lt.0._p_) dz=-step
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/(bpolval*r_mid)*dl_pol

     r(j+1)=r(j)+dr
     z(j+1)=z(j)+dz
     phi(j+1)=phi(j)+dphi

     call  check_whether_field_line_touch_boundary(r(j+1),z(j+1),phi(j+1),x_lcfs,z_lcfs,np_lcfs,loss)
     if (loss.eqv. .true.) exit

     if(abs(floor(abs(phi(j)-phi0)/twopi)-floor(abs(phi(j+1)-phi0)/twopi)).eq.1) then ! finish one toroidal turn
        !   write(*,*) 'j=',j,'k=',k, phi(j),phi(j+1)
        k=k+1
        r_poincare(k)=(r(j)+r(j+1))/two
        z_poincare(k)=(z(j)+z(j+1))/two
        phi_poincare(k)=((phi(j)+phi(j+1))/two)/twopi
     endif

     if(abs(phi0-phi(j+1))/twopi.ge.n_tor_loop) exit

  enddo

 if(j.eq.npt) then
     open(76,file='bad_line.txt')
     do jj=1,j-1
     write(76,*) r(jj),z(jj),phi(jj)
     enddo
     close(76)
     call safety_factor_a_field_line(r,z,phi,j,qval)
     call draw_magnetic_surface(r0,z0,'ref_field_line.txt') !draw the magnetic surface which passes through (r0,z0)
     stop 'max number of tracing steps of field line is exceeded before achiving the specified number of toroidal loop'
endif

  nloop_actual=k
  write(*,*) 'nloop_actual=',nloop_actual, 'actual step along field line=',j

  call safety_factor_a_field_line(r,z,phi,j+1,qval)

!!$  open(76,file='field_line.txt')
!!$  do jj=1,j+1
!!$     write(76,*) r(jj),z(jj),phi(jj)
!!$  enddo
!!$  close(76)


!call check_field_line_in_field_aligned_coordinates(r,z,phi,j,qval)

end subroutine field_line_tracing


subroutine field_line_tracing_simplified(r0,z0,phi0) !see the comments in subroutine "field_line_tracing"
  use constants,only:p_,two,twopi,one_half
use magnetic_field, only :  br,bz,bphi
  implicit none
  real(p_),intent(in):: r0,z0,phi0
  integer,parameter::npt=3000
  real(p_):: r(npt),z(npt),phi(npt)
  real(p_),parameter:: step=1d-3 !meter, trial of dr or dz step
  real(p_):: brval,bzval,bphival,bpolval,dr,dz,dphi
  real(p_):: r_mid,z_mid,dl_pol,qval
  integer:: j,u

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  do j=1,npt-1 !2nd Runge-Kutta
     brval=    br(r(j),z(j))
     bzval=    bz(r(j),z(j))
     bphival=bphi(r(j),z(j))
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr=step*one_half
        if(brval.lt.0._p_) dr=-step*one_half
        dz=bzval/brval*dr
     else
        dz=step*one_half
        if(bzval.lt.0._p_) dz=-step*one_half
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/bpolval*dl_pol/r(j)

     !first step of 2nd Runge-Kutta:
     r_mid=r(j)+dr
     z_mid=z(j)+dz
     brval=    br(r_mid,z_mid)
     bzval=    bz(r_mid,z_mid)
     bphival=bphi(r_mid,z_mid)
     bpolval=sqrt(brval**2+bzval**2)

     if(abs(bzval).lt.abs(brval)) then
        dr=step
        if(brval.lt.0._p_) dr=-step
        dz=bzval/brval*dr
     else
        dz=step
        if(bzval.lt.0._p_) dz=-step
        dr=brval/bzval*dz
     endif

     dl_pol=sqrt(dr**2+dz**2)
     dphi=bphival/(bpolval*r_mid)*dl_pol

     r(j+1)=r(j)+dr
     z(j+1)=z(j)+dz
     phi(j+1)=phi(j)+dphi

  enddo

  open(newunit=u,file='field_line')
  do j=1,npt
     write(u,*) phi(j),z(j),r(j)
  enddo
  close(u)
end subroutine field_line_tracing_simplified


subroutine field_line_tracing0(r0,z0,phi0, npt, dphi, r,z,phi)
  use constants,only : p_, two, twopi, one_half
  use magnetic_field, only :  br, bz, bphi
  implicit none
  real(p_),intent(in) :: r0,z0,phi0
  real(p_),intent(in) :: dphi
  integer, intent(in) :: npt
  real(p_), intent(out) :: r(npt),z(npt),phi(npt)
  real(p_) :: brval,bzval,bphival
  real(p_) :: dr, dz, r_mid, z_mid
  integer :: j

  r(1)=r0
  z(1)=z0
  phi(1)=phi0

  do j=1,npt-1 !2nd Runge-Kutta
     brval=    br(r(j),z(j))
     bzval=    bz(r(j),z(j))
     bphival=bphi(r(j),z(j))

     dr= brval/bphival*(r(j)*dphi)
     dz= bzval/bphival*(r(j)*dphi)

     r_mid=r(j)+dr*one_half !first step of 2nd Runge-Kutta
     z_mid=z(j)+dz*one_half
     brval=    br(r_mid,z_mid)
     bzval=    bz(r_mid,z_mid)
     bphival=bphi(r_mid,z_mid)

     dr= brval/bphival*(r_mid*dphi)
     dz= bzval/bphival*(r_mid*dphi)
    
     r(j+1)=r(j)+dr !second step of 2nd Runge-Kutta
     z(j+1)=z(j)+dz
     phi(j+1)=phi(j)+dphi
  enddo
end subroutine field_line_tracing0

subroutine check_field_line_in_field_aligned_coordinates(r,z,phi,npt,qval)
  use constants,only:p_
  use constants,only: two,twopi
  use radial_module,only:z_axis
  use mapping_module,only: nx_mapping ,j0,r_cyl,tor_shift_b
  use map_to_mc, only : interpolate_from_cylindrical_to_magnetic_coordinates
  use magnetic_field, only : psi_func, pfn_func, radcor_as_func_of_pfn
  implicit none
  integer,intent(in):: npt
  real(p_),intent(in):: r(npt),z(npt),phi(npt),qval
  real(p_):: radcor(npt),theta(npt),alpha(npt),tor_shift(npt)
  real(p_)::x1,x2,y1,y2,z1,z2,dx,dy,dz,dl(npt)
  real(p_):: sum=0._p_,real_shift,total_shift,tmp_array(nx_mapping)
  integer:: j,kk



  do j=1,npt
     radcor(j)=radcor_as_func_of_pfn(pfn_func(r(j),z(j))) !get radial coordinate
     call interpolate_from_cylindrical_to_magnetic_coordinates(r(j),z(j),theta(j),tor_shift(j))
  enddo

  dl(1)=0._p_
  do j=2,npt
     x1=r(j-1)*cos(phi(j-1))
     x2=r(j)*cos(phi(j))
     y1=r(j-1)*sin(phi(j-1))
     y2=r(j)*sin(phi(j))
     z1=z(j-1)
     z2=z(j)
     dx=x2-x1
     dy=y2-y1
     dz=z2-z1
     dl(j)=dl(j-1)+sqrt(dx*dx+dy*dy+dz*dz)
  enddo

  open(11,file='field_line_in_field_aligned_co.txt')
  do j=2,npt

!!$     if(abs(theta(j)-theta(j-1)) .ge. twopi*0.9) then !indecate finishing one poloidal loop
!!$        !total_shift=tor_shift(j)*(z_axis-z(j-1))/(z(j)-z_axis)+tor_shift(j-1)
!!$        !total_shift=qval*twopi*1.0004
!!$        !total_shift=2.35227*twopi
!!$        do kk=1,nx_mapping
!!$           tmp_array(kk)= tor_shift_b(kk,j0)
!!$        enddo
!!$        call linear_1d_interpolate(nx_mapping,r_cyl,tmp_array,r(j),total_shift)
!!$        sum=sum+total_shift
!!$     endif
!!$     alpha(j)=phi(j)-(tor_shift(j)+sum)
   call accumulate_tor_shift(theta(j-1),theta(j),r(j),tor_shift(j),real_shift)
    alpha(j)=phi(j)-real_shift 
     write(11,*) dl(j),radcor(j),theta(j),alpha(j), tor_shift(j), phi(j),r(j),z(j)
  enddo
  close(11)
end subroutine check_field_line_in_field_aligned_coordinates


subroutine safety_factor_a_field_line(r,z,phi,npt,qval)
  !given a field line, calculate its safety factor
  use constants, only: p_, two, twopi
  use magnetic_field, only : psi_func
  use magnetic_field, only : qfunc0
  implicit none
  integer, intent(in) :: npt
  real(p_), intent(in) :: r(npt), z(npt), phi(npt)
  real(p_), intent(out) :: qval
  real(p_) :: phi_old
  integer :: j, npass

  npass = 0
  do j=1,npt-1
     if(z(j)*z(j+1).lt.0) then !indicates one midplane-crossing
        npass = npass+1
        if(npass == 1) phi_old=(phi(j)+phi(j+1))/two
     endif
     if(npass == 3) then !indicates that the line has finished one poloidal period
        qval=abs(phi_old-phi(j))/twopi
        write(*,*) 'safety factor of field line passing (r,z)',r(1),z(1),'is', qval,&
             & 'q value specified in gfile =', qfunc0(psi_func(r(1),z(1)))
        exit
     endif
  enddo

  write(*,*) 'toroidal loops the field line travels=',(phi(npt)-phi(1))/twopi

end subroutine safety_factor_a_field_line


subroutine check_whether_field_line_touch_boundary(r,z,phi,rlim,zlim,nlim,loss)
  use constants,only:p_
  use math, only : pnpoly
!  use boundary,only: nlim,rlim,zlim 
  implicit none
  real(p_),intent(in):: r,z,phi
  integer,intent(in):: nlim
  real(p_),intent(out):: rlim(nlim),zlim(nlim)
  logical,intent(out):: loss
  integer:: inout

  call PNPOLY(r,z,rlim,zlim,nlim,INOUT) !find out wheter the point (r,z) is within the limiter
  !        if (inout.eq.1) then !within the LCFS
  if (inout.eq.-1 .or.inout.eq.0) then !the particle is out of the limiter
     write(*,*) '==>This field line touches the limiter at (R,Z,phi)=', r,z,phi
     loss=.true.
     !stop
  else
     loss=.false.
  endif

end subroutine


  subroutine accumulate_tor_shift(theta_old,theta_new,r,tor_shift,real_shift) !,kt)
    use constants,only:p_
    use constants,only: twopi
    use mapping_module,only: nx_mapping,tor_shift_b,j0,r_cyl
  use interpolate_module,only: linear_1d_interpolate
    implicit none
    real(p_),intent(in):: theta_old,theta_new,r,tor_shift
    real(p_),intent(out)::real_shift
    real(p_),save::sum=0._p_
    real(p_):: tmp_array(nx_mapping),twopi_q
!integer,intent(in):: kt
    integer::kk
    if(abs(theta_old-theta_new) .ge. twopi*0.9) then !indicate finishing one poloidal loop
       do kk=1,nx_mapping
          tmp_array(kk)= tor_shift_b(kk,j0)
       enddo
       call linear_1d_interpolate(nx_mapping,r_cyl,tmp_array,r,twopi_q) !the result is twopi*q, I use this instead of directly using twopi*q because the latter may cause some cancellation problem
       sum=sum+twopi_q*sign(1._p_,theta_old-theta_new)
!write(*,*) 'sum=',sum, 'twopi_q=',twopi_q,'tor_shift=',tor_shift,'sum+tor_shift/ real_shift',sum+tor_shift, 'kt=',kt
    endif

    real_shift=sum+tor_shift


    !real_shift=sum-(twopi_q-tor_shift)
  end subroutine accumulate_tor_shift


  subroutine draw_magnetic_surface(r0,z0,filename) !draw the magnetic surface which passes through (r0,z0)
  use constants,only:p_
  use constants,only:zero,one,two,twopi
  use boundary,only:np_lcfs,x_lcfs,z_lcfs
  use radial_module,only: r_axis,z_axis
  use contour_mod,only : find_contour
  use magnetic_field, only : psi_func
  implicit none

  real(p_),intent(in):: r0,z0
  character(*),intent(in)::  filename
  real(p_):: psival
  real(p_):: x_contour(np_lcfs),z_contour(np_lcfs)
  integer:: i,u
  
  psival=psi_func(r0,z0)

  call find_contour(psival,x_contour,z_contour)

  open(newunit=u,file=filename)
  do i=1,np_lcfs
     write(u,*) x_contour(i),z_contour(i)
  enddo
  close(u)
end subroutine draw_magnetic_surface
module deposit_gk_module
contains
  
  subroutine deposit_gk(ns, lost, vpar_gk, w_gk, x_ring, y_ring, z_ring, &
       & density_left, density_right, jpar_left, jpar_right)
    !with markers' (v,x,w) given, do the deposition to get density and jpar on grids
    use constants, only: p_, one
    use magnetic_coordinates, only: n=>nrad,m=>mtor, xgrid, ygrid, dtor, dradcor, xlow, xupp
    use domain_decomposition, only: dtheta2, theta_start
    use gk_module, only: nm_gk, gk_flr, gyro_npt, charge_gk, vn_gk
    use normalizing, only: vu,  qu
    implicit none
    integer,intent(in)  :: ns
    logical, intent(in) :: lost(:)
    real(p_),intent(in) :: vpar_gk(:), w_gk(:), x_ring(:,:), y_ring(:,:), z_ring(:,:)
    real(p_), intent(inout) :: density_left(:,:), density_right(:,:)
    real(p_), intent(inout) :: jpar_left(:,:), jpar_right(:,:)
    real(p_) :: cz1, cz2, cy1, cy2, cx1, cx2, kernel, kernel2
    integer  :: i,j,k, iz, i_plus1, j_plus1, kr, npt0, nm

    nm= nm_gk(ns)
    
    if(gk_flr(ns) .eqv. .false.) then
       npt0 = 1
    else
       npt0 = gyro_npt
    endif

    do k = 1, nm
       if (lost(k) .eqv. .true.) cycle
       kernel = (charge_gk(ns)/qu)*w_gk(k)/npt0 !for computing charge density
       kernel2 = kernel*vpar_gk(k)*(vn_gk(ns)/vu) !for computing parallel current

       do kr = 1, npt0 !loop over the points on a gyro-ring
          i=floor((y_ring(kr,k)-ygrid(1))/dtor+1) 
          cy1=(y_ring(kr,k)-ygrid(i))/dtor
          cy2=one-cy1

          j=floor((x_ring(kr,k)-xgrid(1))/dradcor+1)
          cx1= (x_ring(kr,k)-xgrid(j))/dradcor
          cx2=one-cx1

          cz1=(z_ring(kr,k)-theta_start)/dtheta2
          cz2=one-cz1

          i_plus1=i+1
          if(i.eq.m) i_plus1=1 !periodic condition
          j_plus1=j+1
          if(j.eq.n) cycle ! particle is out of radial computational region

          density_left(i,j) = density_left(i,j) +kernel*cz2*cy2*cx2
          density_left(i_plus1,j) = density_left(i_plus1,j)  +kernel*cz2*cy1*cx2
          density_left(i,j_plus1) = density_left(i,j_plus1) +kernel*cz2*cy2*cx1
          density_left(i_plus1,j_plus1) = density_left(i_plus1,j_plus1)+kernel*cz2*cy1*cx1

          density_right(i,j) = density_right(i,j)+kernel*cz1*cy2*cx2
          density_right(i_plus1,j) = density_right(i_plus1,j)+kernel*cz1*cy1*cx2
          density_right(i,j_plus1) = density_right(i,j_plus1)+kernel*cz1*cy2*cx1
          density_right(i_plus1,j_plus1) = density_right(i_plus1,j_plus1)+kernel*cz1*cy1*cx1

          jpar_left(i,j) = jpar_left(i,j) + kernel2*cz2*cy2*cx2
          jpar_left(i_plus1,j) = jpar_left(i_plus1,j) + kernel2*cz2*cy1*cx2
          jpar_left(i,j_plus1) = jpar_left(i,j_plus1) + kernel2*cz2*cy2*cx1
          jpar_left(i_plus1,j_plus1) = jpar_left(i_plus1,j_plus1)+kernel2*cz2*cy1*cx1

          jpar_right(i,j) = jpar_right(i,j) + kernel2*cz1*cy2*cx2
          jpar_right(i_plus1,j) = jpar_right(i_plus1,j)+ kernel2*cz1*cy1*cx2
          jpar_right(i,j_plus1) = jpar_right(i,j_plus1) + kernel2*cz1*cy2*cx1
          jpar_right(i_plus1,j_plus1) = jpar_right(i_plus1,j_plus1)+kernel2*cz1*cy1*cx1          
       enddo

!!$       do kr = 1, npt0 !nearest neighbour interpolation, turns out to give more noisy results
!!$          i = nint((y_ring(kr,k)-ygrid(1))/dtor)+1
!!$          j = nint((x_ring(kr,k)-xgrid(1))/dradcor)+1
!!$          iz = nint((z_ring(kr,k)-theta_start)/dtheta2)
!!$
!!$          if(iz==0) then
!!$             density_left(i,j) = density_left(i,j) + kernel
!!$             jpar_left(i,j) = jpar_left(i,j) + kernel2
!!$          else
!!$             density_right(i,j) = density_right(i,j) + kernel
!!$             jpar_right(i,j) = jpar_right(i,j) + kernel2
!!$          endif
!!$       enddo

    enddo
  end subroutine deposit_gk

end module deposit_gk_module
module deposit_fk_module
  implicit none
contains
  subroutine deposit_fk(nmarker_i,active_i,radcor_i,theta_i,alpha_i,phi_i,w_i) !given (r,v,w) of markers, do the deposition to get perpendicular currents and number density on spatial grids
    !periodic toroidal boundary condition and the poloidal boundary condtion (connection condition along the field line) are taken into account
    !    use mpi
    use constants,only:p_
    use normalizing,only: nu
    use constants,only: one,two
    use fk_module,only: ni0, my_den_i_left, my_den_i_right !as output
    use magnetic_coordinates,only:m=>mtor,n=>nrad
    use magnetic_coordinates,only:xgrid,zgrid,ygrid,dtor,dradcor,dtheta,jacobian
    !use magnetic_coordinates,only:phi_grid_left,phi_grid_right
    use domain_decomposition,only: dtheta2,theta_start,ipol_eq, multi_eq_cells


    integer,intent(in)::nmarker_i
    logical,intent(in):: active_i(nmarker_i)
    real(p_),intent(in)::radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i),phi_i(nmarker_i)
    real(p_),intent(in):: w_i(nmarker_i)
    real(p_):: coeff_theta_1,coeff_theta_2,coeff_alpha_1,coeff_alpha_2,coeff_radcor_1,coeff_radcor_2
    real(p_):: dv1, dv2
    integer:: k,i,j, jeq,i_plus_one,j_plus_one

    !before the deposition, set the array to zero
    my_den_i_left=0._p_
    my_den_i_right=0._p_

    !$omp parallel do private(coeff_theta_1,coeff_theta_2, coeff_alpha_1, coeff_alpha_2, &
    !$omp & coeff_radcor_1, coeff_radcor_2,i,j,i_plus_one,j_plus_one,dphi_left00,dphi_left01,dphi_left10,dphi_left11,&
    !$omp & dphi_right00, dphi_right01, dphi_right10,dphi_right11,kernel) !!&
!!!!$omp & reduction(+:my_jr_left,my_jr_right,my_jphi_left,my_jphi_right,my_jz_left,my_jz_right,my_den_i_left,my_den_i_right)
    do k=1,nmarker_i !for each marker, deposit it to the corresponding grids
       if(active_i(k).eqv..false.) cycle ! markers outside the computational region do not contribute density or current to any grids
       !determine the interpolating coeefficeint according to the location of a marker
       coeff_theta_1=(theta_i(k)-theta_start)/dtheta2
       coeff_theta_2=one-coeff_theta_1

       !call location(m,ygrid,alpha_i(k),i)
       i=floor((alpha_i(k)-ygrid(1))/dtor+1) !uniform xarray is assumed, otherwise we need to call location() subroutine to locate xval
       coeff_alpha_1=(alpha_i(k)-ygrid(i))/dtor
       coeff_alpha_2=one-coeff_alpha_1
       !if(myid.eq.2) write(*,*) 'alpha, i=',i, 'k=',k
       !call location(n, xgrid, radcor_i(k),j)
       j=floor((radcor_i(k)-xgrid(1))/dradcor+1)
       coeff_radcor_1= (radcor_i(k)-xgrid(j))/dradcor
       coeff_radcor_2=one-coeff_radcor_1

       i_plus_one=i+1
       if(i.eq.m) i_plus_one=1 !periodic condition in the toroidal direction
       j_plus_one=j+1
       if(j.eq.n) j_plus_one=j !all the active markers are actually within the boundary flux surface labeled by n. If j.eq.n, then the marker is exactly on the boundary flux surface, this code line is needed to avoid exeeding array bounds in case that the marker is exactly on the boundary flux surface or slightly outside of it due to numerical truncation errors
       !if(j.eq.n) j_plus_one=1 !periodic condition. is disabled later in this subroutine, now using fixed zero boundary condition along the radial direction. In the past, I used the periodic bounary condition because I want to use Fourier transform along the radial direction. However this periodic radial boundary condition is not reasonable. The reason is as follows: Here the radial direction is d/dpsi in (psi,theta,alpha) coordinates. This direction is a combinition of the toroidal and the usual radial direction and thus is an artificial direction introduced when using (psi,theta,alpha) coordinates, and does not has any physical reason to justify that a peroidic condition should be satisfied along this artificial direction.

       !density, only used in isotropic fluid electrons model or adiabatic electrons model to provide delta_ne, which is assumed to be equal to delta_ni
       my_den_i_left(i,j)=my_den_i_left(i,j) +w_i(k)*coeff_theta_2*coeff_alpha_2*coeff_radcor_2
       my_den_i_right(i,j)=my_den_i_right(i,j)+w_i(k)*coeff_theta_1*coeff_alpha_2*coeff_radcor_2
       my_den_i_left(i_plus_one,j)= my_den_i_left(i_plus_one,j)+w_i(k)*coeff_theta_2*coeff_alpha_1*coeff_radcor_2
       my_den_i_right(i_plus_one,j)=my_den_i_right(i_plus_one,j)+w_i(k)*coeff_theta_1*coeff_alpha_1*coeff_radcor_2
       my_den_i_left(i,j_plus_one)=my_den_i_left(i,j_plus_one)+w_i(k)*coeff_theta_2*coeff_alpha_2*coeff_radcor_1
       my_den_i_right(i,j_plus_one)=my_den_i_right(i,j_plus_one)+w_i(k)*coeff_theta_1*coeff_alpha_2*coeff_radcor_1
       my_den_i_left(i_plus_one,j_plus_one)=my_den_i_left(i_plus_one,j_plus_one)+w_i(k)*coeff_theta_2*coeff_alpha_1*coeff_radcor_1
       my_den_i_right(i_plus_one,j_plus_one)=my_den_i_right(i_plus_one,j_plus_one)+w_i(k)*coeff_theta_1*coeff_alpha_1*coeff_radcor_1
    enddo
    !$omp end parallel do

    do j=1,n  !divided by the space volume of a cell, to give the current denisty, note that this 'cell' is the cell defined by PIC (i.e., grid is the center of the cell). (while grids are the boundaries of the "mpi cell" for grouping and pushing particles)
       jeq=j
       dv1=abs(jacobian(ipol_eq,jeq))*dradcor*dtheta2*dtor !volume of the cell (the center of the cell is the grid)
       dv2=abs(jacobian(ipol_eq+multi_eq_cells,jeq))*dradcor*dtheta2*dtor !volume of the cell
       my_den_i_left(:,j)  = my_den_i_left (:,j)/dv1
       my_den_i_right(:,j) = my_den_i_right(:,j)/dv2
    enddo

    my_den_i_left  = my_den_i_left /nu
    my_den_i_right = my_den_i_right/nu
  end subroutine deposit_fk
end module deposit_fk_module
module communication_connection
  !communicate field value between neighbour cells, and handle the connecting condition across the theta cut
  use constants, only: p_
  implicit none
contains

  subroutine connect_condition_across_theta_cut(u, direction)
    use interpolate_module, only: linear_1d_interpolate
    use magnetic_coordinates, only: ygrid, tor_shift_mc, mpol, toroidal_range
    use math, only: shift_toroidal
    real(p_), intent(inout) :: u(:,:)
    integer, intent(in) :: direction ! alowable values: +1 or -1
    !If direction == +1, then follow magnetic field line along the positive direction of theta
    !If direction==-1, then follow magnetic field line along the negative direction of theta
    real(p_), allocatable :: u_old(:)
    real(p_) :: y0, twopiq
    integer :: i, j, m, n

    m = size(u,1) !toroidal
    n = size(u,2) !radial
    allocate(u_old(m+1))

    do j = 1, n !radial
       u_old(1:m) = u(:,j)
       u_old(m+1) = u(1,j) ! toroidal periodic condition

       twopiq = tor_shift_mc(mpol,j) - tor_shift_mc(1,j)
       do i = 1, m
          y0 = ygrid(i) + twopiq*direction
          call shift_toroidal(y0, toroidal_range) 
          call linear_1d_interpolate(m+1, ygrid, u_old, y0, u(i,j))  
       enddo
    enddo
  end subroutine connect_condition_across_theta_cut


  subroutine merge_source(my_jpar_left, my_jpar_right, jpar)
    use magnetic_coordinates, only: m=>mtor, n=>nrad
    use domain_decomposition, only: GRID_COMM, TUBE_COMM, GCLR, my_left, my_right
    use mpi
    real(p_), intent(in) :: my_jpar_left(m,n), my_jpar_right(m,n)
    real(p_), intent(out) :: jpar(m,n)
    real(p_) :: jpar_left(m,n), jpar_right(m,n), jpar_left0(m,n)
    integer :: status(MPI_STATUS_SIZE), ierr

    !summing over all those procs in a z-cell
    call MPI_ALLREDUCE(my_jpar_left,  jpar_left,  m*n, MPI_REAL8, MPI_SUM, GRID_COMM, ierr)
    call MPI_ALLREDUCE(my_jpar_right, jpar_right, m*n, MPI_REAL8, MPI_SUM, GRID_COMM, ierr)

    call MPI_Sendrecv(jpar_right, m*n, MPI_real8, my_right, 3, &
         &            jpar_left0, m*n, MPI_real8, my_left,  3, Tube_comm, status, ierr)

    if(GCLR==0) call connect_condition_across_theta_cut(jpar_left0, direction=-1) 
    jpar(:,:) = jpar_left(:,:) + jpar_left0(:,:) !add the contribution from the neighbour cell

  end subroutine merge_source

  subroutine update_scalar_at_right_boundary(s) 
    !Every proc is response for one cell which has two boundary grids.
    !Only the field on left-boundary-grid is computed by the present proc.
    !The field on the right-boundary is received from the neighbour proc. This subroutine handle this communication and the edge case for gclr=mpol2-1
    !Note that the definition of cell here is different from the the definition of cell in PIC.
    !The grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.
    use mpi
    use domain_decomposition, only: myid, GCLR, Tube_comm, my_left, my_right
    use magnetic_coordinates, only : mpol2
    real(p_), intent(inout) :: s(:,:,:)
    integer :: status(MPI_STATUS_SIZE), ierr, m, n

    m = size(s,1) !toroidal
    n = size(s,2) !radial
    call MPI_Sendrecv(s(:,:,1),  m*n,  MPI_real8, my_left,  41, &
         &            s(:,:,2),  m*n,  MPI_real8, my_right, 41, Tube_COMM, status, ierr)
    !special treatment at theta cut, to handle the phi-grid mismatch
    if(GCLR==mpol2-1) call connect_condition_across_theta_cut(s(:,:,2), direction=1)

  end subroutine update_scalar_at_right_boundary


  subroutine update_derivatives_at_right_boundary(ax, ay, az)
    !communication between neighbour cells: Every proc is response for one cell which has two boundary grids.
    !Only the field on left-boundary-grid is computed by the present proc.
    !So the value on the right boundary is obtained by communicating with neighbour process.
    !field value at right-boundary of the present cell is needed when pushing particle weights
    use mpi
    use magnetic_coordinates, only: m=>mtor, n=>nrad, mpol2
    use domain_decomposition, only: myid, GCLR, Tube_comm, ntube, my_left, my_right
    real(p_), intent(inout) ::  ax(:,:,:), ay(:,:,:), az(:,:,:)
    integer :: status(MPI_STATUS_SIZE), ierr

    call MPI_Sendrecv(ax(:,:,1), m*n,  MPI_real8, my_left,  4, &
         &            ax(:,:,2), m*n,  MPI_real8, my_right, 4, Tube_COMM,status,ierr)

    call MPI_Sendrecv(ay(:,:,1), m*n,  MPI_real8, my_left,  5, &
         &            ay(:,:,2), m*n,  MPI_real8, my_right, 5, Tube_COMM,status,ierr)

    call MPI_Sendrecv(az(:,:,1), m*n,  MPI_real8, my_left,  6, &
         &            az(:,:,2), m*n,  MPI_real8, my_right, 6, Tube_COMM,status,ierr)

    if(GCLR == mpol2-1) then !special treatment at theta cut
       call connect_condition_across_theta_cut(ax(:,:,2), direction=1)
       call connect_condition_across_theta_cut(ay(:,:,2), direction=1) 
       call connect_condition_across_theta_cut(az(:,:,2), direction=1) 
       call derivative_transform_at_theta_cut(ax(:,:,2), ay(:,:,2)) !transform the components at theta=-pi to theta=+pi
    endif
  end subroutine update_derivatives_at_right_boundary


  subroutine derivative_transform_at_theta_cut(ax, ay)
    !Transform the components of gradient of a scalar field from theta=-pi to theta=pi
    !The basis vectors grad_alpha, in terms of which the gradient of apara/potential are decomposed, is discontinuous across the theta-cut. Therefore the components at one side need to be transformed to the components on the other side
    use magnetic_coordinates, only: mpol, nrad, grad_psi, grad_psi_dot_grad_alpha, &
         & grad_psi_r, grad_psi_z, grad_alpha_r, grad_alpha_z
    real(p_), intent(in) :: ay(:,:)
    real(p_), intent(inout) :: ax(:,:)
    real(p_) :: grad_psi0, gxdgy0, gxdgy1 !gxdgy1 is the value of gxdgy at theta=+pi
    integer :: j

    do j = 1, nrad
       grad_psi0 = grad_psi(1,j)
       gxdgy0 = grad_psi_dot_grad_alpha(1,j)
       gxdgy1 = grad_psi_r(1,j) * grad_alpha_r(mpol,j)  &
            & + grad_psi_z(1,j) * grad_alpha_z(mpol,j)

       ax(:,j) = ax(:,j) + ay(:,j)*(gxdgy0-gxdgy1)/grad_psi0**2
    enddo
  end subroutine derivative_transform_at_theta_cut


  subroutine communicate_field_value_between_neighbour_cells2() !use cylindrical components of the electric field
    use mpi
    use perturbation_field,only: ef_cyl_r_left, ef_cyl_z_left, ef_cyl_phi_left !already known before entering this subroutine
    use perturbation_field,only: ef_cyl_r_right, ef_cyl_z_right, ef_cyl_phi_right !as output
    use magnetic_coordinates,only: m=>mtor, n=>nrad, mpol2
    use domain_decomposition,only: myid, GCLR,Tube_comm,ntube, my_left,my_right

    integer:: status(MPI_STATUS_SIZE),ierr
    !communication between neighbour cells: Every proc is response for one cell which has two boundary grids, only the field on left-boundary-grid is computed by the present proc, the field on the right-boundary is received from the neighbour proc. Note that the definition of cell here is different from the the definition of cell in PIC: the grids are the centers of the cells defined in PIC while the grids are the boundaries of the cells mentioned here.

    call MPI_Sendrecv(ef_cyl_r_left,  (m+1)*n,  MPI_real8, my_left,  4,&
         &            ef_cyl_r_right, (m+1)*n,  MPI_real8, my_right, 4,Tube_COMM,status,ierr)
    call MPI_Sendrecv(ef_cyl_z_left,  (m+1)*n,  MPI_real8, my_left,  5,&
         &            ef_cyl_z_right, (m+1)*n,  MPI_real8, my_right, 5,Tube_COMM,status,ierr)
    call MPI_Sendrecv(ef_cyl_phi_left, (m+1)*n,  MPI_real8, my_left,  6,&
         &            ef_cyl_phi_right,(m+1)*n,  MPI_real8, my_right, 6,Tube_COMM,status,ierr)

    if(GCLR == mpol2-1) then !special treatment at theta cut
       call connect_condition_across_theta_cut(ef_cyl_r_right, direction=1) 
       call connect_condition_across_theta_cut(ef_cyl_z_right, direction=1) 
       call connect_condition_across_theta_cut(ef_cyl_phi_right, direction=1) 
    endif
  end subroutine communicate_field_value_between_neighbour_cells2


  subroutine get_nearby_field_along_field_line(a, a_left, a_right, a_left2, a_right2)
    !get value of field on the two grids that are to the left/right of the present grid
    use constants, only: p_
    use magnetic_coordinates, only: mpol2
    use domain_decomposition, only: GCLR, TUBE_COMM, my_left, my_right, my_left2, my_right2
    use mpi
    implicit none
    real(p_), intent(in) :: a(:,:)
    real(p_), intent(out) :: a_left(:,:), a_right(:,:)
    real(p_), optional, intent(out) :: a_left2(:,:), a_right2(:,:)
    integer :: status(MPI_STATUS_SIZE), ierr, m, n

    m=size(a,1)
    n=size(a,2)

    call MPI_Sendrecv(a,      m*n, MPI_real8, my_right, 1, &
         &            a_left, m*n, MPI_real8, my_left, 1, Tube_COMM, status, ierr)

    call MPI_Sendrecv(a,       m*n, MPI_real8, my_left, 2, &
         &            a_right, m*n, MPI_real8, my_right, 2,Tube_COMM, status,ierr)


    if(GCLR==mpol2-1) then
       call connect_condition_across_theta_cut(a_right, 1)
    endif

    if(GCLR==0) then
       call connect_condition_across_theta_cut(a_left, -1)
    endif

    if(present(a_left2) .neqv. present(a_right2)) stop "a_left2 and a_right2 must be both present or both absent"
    if( .not. present(a_left2)) return
    call MPI_Sendrecv(a,       m*n, MPI_real8, my_right2, 3, &
         &            a_left2, m*n, MPI_real8, my_left2,  3,Tube_COMM,status,ierr)

    call MPI_Sendrecv(a,        m*n, MPI_real8, my_left2,  4, &
         &            a_right2, m*n, MPI_real8, my_right2, 4,Tube_COMM,status,ierr)

    if(GCLR==mpol2-2) then
       call connect_condition_across_theta_cut(a_right2, 1)
    endif

    if(GCLR==mpol2-1) then
       call connect_condition_across_theta_cut(a_right2, 1)
    endif

    if(GCLR==1) then
       call connect_condition_across_theta_cut(a_left2, -1)
    endif
    if(GCLR==0) then
       call connect_condition_across_theta_cut(a_left2, -1)
    endif

  end subroutine get_nearby_field_along_field_line

  
end module communication_connection
module filter_module
  implicit none
  integer, parameter ::  radial_harmonics_included = 25
contains

  subroutine surface_average_of_n0(s, nx)
    ! retain only the magnetic surface average of n=0 harmonic (other poloidal harmonics are filtered out)
    use constants, only: p_, twopi
    use magnetic_coordinates, only: nrad, nz => mpol2, jacobian_av, abs_jacobian
    use domain_decomposition, only: tube_comm, ipol_eq, dtheta2
    use mpi
    implicit none
    integer, intent(in) :: nx !radial
    complex(p_), intent(inout) :: s(nx) ! Amplitude of n=0 Fourier harmonic
    real(p_) :: my_sz(nx), sz(nx, 0:nz-1), sx(nx)
    integer :: j, ierr

    my_sz(:) = real(s(:))*abs_jacobian(ipol_eq, 2:nrad-1)
    call mpi_allgather(my_sz(:), nx, MPI_real8, sz, nx, MPI_real8, tube_comm, ierr)

    do j = 1, nx
       sx(j) = sum(sz(j,:))*dtheta2/(jacobian_av(j+1)*twopi)
    enddo

    s(:) = sx(:) !store the result in the complex array.
  end subroutine 
  
  subroutine mfilter_for_n0(s, n)
    ! The n=0 mode is often numerically unstable if we keep all m harmonics
    ! m filtering (usually retain only m=0) is needed to make the simulation stable
    use constants, only: p_, ii, twopi, pi
    use magnetic_coordinates, only: nz => mpol2
    use domain_decomposition, only: tclr, gclr, ntube, grid_comm, tube_comm
    use mpi
    implicit none
    integer, intent(in) :: n !radial
    complex(p_), intent(inout) :: s(n)
    integer, parameter :: mupp = 0 !poloidal harmonics m in [-mupp:mupp], basis funtions: exp(ii*m*theta)
    logical, save :: is_first = .true.
    integer, save :: my_start, my_end, spacing, my_range
    integer, allocatable, save :: recvcounts(:), displacement(:)
    complex(p_), allocatable, save :: wexp1(:,:), wexp2(:)
    complex(p_), allocatable :: coef(:,:)
    real(p_), allocatable :: my_sz(:), sz(:,:), my_sx(:), sx(:)
    integer :: i, j, k, ierr

    if (is_first .eqv. .true.) then !needs to be genreated only for the first time
       is_first = .false.
       spacing = n/ntube
       my_start = TCLR*spacing + 1
       my_end = my_start + (spacing-1)
       if(TCLR == ntube-1) my_end = n !the last process handles all the remainder if n is not a perfect multiple of ntube
       my_range = my_end - my_start + 1
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:) = spacing
       recvcounts(ntube-1) = recvcounts(ntube-1)+(n-spacing*ntube) !last process contains additional elements
       do i = 0, ntube-1
          displacement(i) = i * spacing
       enddo
       allocate(wexp1(0:nz-1, -mupp:mupp))
       allocate(wexp2(-mupp:mupp))
       do k = -mupp, mupp
          do i = 0, nz-1
             wexp1(i,k) = exp(-ii*k*(-pi+i*twopi/nz))
          enddo
          wexp2(k) = exp(+ii*k*(-pi+gclr*twopi/nz))
       enddo
    endif

    allocate(my_sz(my_start:my_end))
    allocate(sz(my_start:my_end, 0:nz-1))
    allocate(coef(-mupp:mupp, my_start:my_end))
    allocate(my_sx(my_start:my_end))
    allocate(sx(n))

    my_sz(:) = real(s(my_start:my_end))
    call mpi_allgather(my_sz(:), my_range, MPI_real8, sz, my_range, MPI_real8, tube_comm, ierr)

    do j = my_start, my_end
       do k = -mupp, mupp
          coef(k,j) = sum(sz(j,:)*wexp1(:,k)) ! Fourier expansion coefficient
       enddo
       my_sx(j) = real(sum(coef(:,j)*wexp2(:))/nz)  ! Inverse Fourier Transform, and convert the type to real
    enddo

    call mpi_allgatherv(my_sx, my_range, MPI_real8, &
         &     sx, recvcounts, displacement, MPI_real8, grid_comm, ierr)

    s(:) = sx(:) !store the result in the complex array.
  end subroutine mfilter_for_n0

  subroutine mfilter_for_each_n(s, nx, ntor)
    use constants, only: p_, ii, twopi, pi
    use magnetic_coordinates, only: nrad, nz => mpol2
    use radial_module, only: qrad
    use domain_decomposition, only: tclr, gclr, ntube, grid_comm, tube_comm, theta_start
    use mpi
    implicit none
    integer, intent(in) :: nx, ntor !radial, toroidal_harmonic
    complex(p_), intent(inout) :: s(nx)
    integer, parameter :: mupp = 5 !poloidal harmonics m in [-mupp:mupp], bais funtions: exp(ii*m*theta)
    complex(p_) :: coef(-mupp:mupp)
    logical, save :: is_first = .true.
    integer, save :: my_start, my_end, spacing, my_range
    integer, allocatable, save :: recvcounts(:), displacement(:)
    complex(p_), allocatable, save :: wexp1(:,:,:), wexp2(:,:)
    complex(p_), allocatable :: my_sz(:), sz(:,:), my_filtered(:), filtered(:)
    integer :: i, j, k, ierr

    if (is_first .eqv. .true.) then !needs to be genreated only for the first time
       is_first = .false.
       spacing = nx/ntube
       my_start = TCLR*spacing + 1
       my_end = my_start + (spacing-1)
       if(TCLR == ntube-1) my_end = nx !last process handles all the remainder if nx is not a perfect multiple of ntube
       my_range = my_end - my_start + 1
       allocate(recvcounts(0:ntube-1))
       allocate(displacement(0:ntube-1))
       recvcounts(:) = spacing
       recvcounts(ntube-1) = recvcounts(ntube-1)+(nx-spacing*ntube) !last process contains additional elements
       do i = 0, ntube-1
          displacement(i) = i * spacing
       enddo
       allocate(wexp1(0:nz-1, -mupp:mupp, my_start:my_end))
       allocate(wexp2(-mupp:mupp, my_start:my_end))
       do k = -mupp, mupp
          do i = 0, nz-1
             wexp1(i, k, :) = exp(-ii*(k-nint(ntor*qrad(my_start+1:my_end+1)))*(-pi+i*twopi/nz))
          enddo
          wexp2(k, :) = exp(+ii*(k-nint(ntor*qrad(my_start+1:my_end+1)))*(-pi+gclr*twopi/nz))
       enddo
    endif

    allocate(my_sz(my_start:my_end))
    allocate(sz(my_start:my_end, 0:nz-1))
    my_sz(:) = s(my_start:my_end) /exp(ii*ntor*qrad(my_start+1:my_end+1)*(theta_start+pi))
    call mpi_allgather(my_sz(:), my_range, MPI_complex16, sz, my_range, MPI_complex16, tube_comm, ierr)

    allocate(my_filtered(my_start:my_end))
    allocate(filtered(nx))
    do j = my_start, my_end
       do k = -mupp, mupp
          coef(k) = sum(sz(j,:)*wexp1(:,k, j))/nz ! Fourier expansion coefficient
       enddo
       my_filtered(j) = sum(coef(:)*wexp2(:,j))  ! Reconstruction (i.e., Inverse Fourier Transform)
    enddo
    call mpi_allgatherv(my_filtered, my_range, MPI_complex16, &
         &     filtered, recvcounts, displacement, MPI_complex16, grid_comm, ierr)

    s(:) = filtered(:) *exp(ii*ntor*qrad(2:nrad-1)*(theta_start+pi))
  end subroutine 


  subroutine mfilter_for_each_n_tmp(s, nx, ntor) !without raidal parallization, gving the same result as the above
    use constants, only: p_, ii, twopi, pi
    use magnetic_coordinates, only: nrad, nz => mpol2
    use radial_module, only: qrad
    use domain_decomposition, only: gclr, tube_comm, theta_start
    use mpi
    implicit none
    integer, intent(in) :: nx, ntor !radial, toroidal_harmonic
    complex(p_), intent(inout) :: s(nx)
    integer, parameter :: mupp = 5 !poloidal harmonics m in [-mupp:mupp], bais funtions: exp(ii*m*theta)
    logical, save :: is_first = .true.
    complex(p_), allocatable, save :: wexp1(:,:), wexp2(:)
    complex(p_) :: coef(-mupp:mupp,nx)
    complex(p_) :: my_sz(nx), sz(nx, 0:nz-1), sx(nx)
    integer :: i, j, k, ierr

    if (is_first .eqv. .true.) then !needs to be genreated only for the first time
       is_first = .false.
       allocate(wexp1(0:nz-1, -mupp:mupp))
       allocate(wexp2(-mupp:mupp))
       do k = -mupp, mupp
          do i = 0, nz-1
             wexp1(i,k) = exp(-ii*k*(-pi+i*twopi/nz))
          enddo
          wexp2(k) = exp(+ii*k*(-pi+gclr*twopi/nz))
       enddo
    endif

    my_sz(:) = s(:)/exp(ii*ntor*qrad(2:nrad-1)*(theta_start+pi))
    call mpi_allgather(my_sz(:), nx, MPI_complex16, sz, nx, MPI_complex16, tube_comm, ierr)

    do j = 1, nx
       do k = -mupp, mupp
          coef(k,j) = sum(sz(j,:)*wexp1(:,k)) ! Fourier expansion coefficient
       enddo
       sx(j) = sum(coef(:,j)*wexp2(:))/nz  ! Inverse Fourier Transform
    enddo

    s(:) = sx(:) *exp(ii*ntor*qrad(2:nrad-1)*(theta_start+pi))
  end subroutine 



  subroutine toroidal_filter(field,m,n)
    use constants,only:p_
    use magnetic_coordinates,only:nsegment !computational toroidal region is 1/nsegment of the full torus
    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(inout)::   field(0:m-1,0:n-1)
    integer:: nharmonic
    integer:: i,i_negative_wn

    nharmonic=1 !harmonic that needs to be kept in the computational toroidal segment
    !write(*,*) 'nharmonic=',nharmonic
    i_negative_wn=m-nharmonic
    if(nharmonic.eq.0) i_negative_wn=0

    do i=0,m-1
       if((i.ne. nharmonic) .and. (i.ne.i_negative_wn)) field(i,:)=(0._p_,0._p_)
    enddo

!!$    do j=0,n-1
!!$     do i=0,m/2 !scanning over the positive wavenumber
!!$        if(i.eq.nharmonic) then 
!!$           out(i,j)=in(i,j) !positive wavenumber
!!$           i_negative_wn=m-i
!!$           if(i.eq.0) i_negative_wn=0
!!$           out(i_negative_wn,j)=in(i_negative_wn,j) !negative or zero wavenumber
!!$        endif
!!$     enddo
!!$    enddo
  end subroutine toroidal_filter


  subroutine toroidal_reconstruct_v0(s_dft,s,m,n) !manually reconstruct (instead of using inverse DFT) the real space function from the DFT data, the result was verified, which agrees with that given by inverse DFT
    use constants, only: p_, twopi
    use magnetic_coordinates, only: toroidal_range, ygrid, nsegment !computational toroidal region is 1/nsegment of the full torus
    implicit none
    integer, intent(in) :: m,n
    complex(p_), intent(in) :: s_dft(0:m-1,n)
    real(p_), intent(out) :: s(m,n)
    real(p_) :: tor, a, b
    integer :: i, j, i_positive, i_negative
    !complex(p_):: ii=(0._p_,1._p_)
    i_positive=1 !harmonic that needs to be kept in the computational toroidal segment
    i_negative=m-i_positive
    ! write(*,*) 'i_positive=',i_positive,'i_negative=',i_negative
    !s=0._p_
    do j=1,n
       a=2*real(s_dft(i_positive,j))/m
       b=-2*imag(s_dft(i_positive,j))/m
       !a=s_dft(i_positive,j)+s_dft(i_negative,j)
       !b=(s_dft(i_positive,j)-s_dft(i_negative,j))*ii
       do i=1,m
          tor=ygrid(i)
          s(i,j)=a*cos(i_positive*twopi*tor/toroidal_range)+b*sin(i_positive*twopi*tor/toroidal_range)
       enddo
    enddo
  end subroutine toroidal_reconstruct_v0


  subroutine toroidal_reconstruct(s_dft,s,m,n) !manually reconstruct (instead of using inverse DFT) the real space function from the DFT data, the result was verified, which agrees with that given by inverse DFT
    use constants,only:p_
    use constants,only:twopi
    use magnetic_coordinates,only:toroidal_range,ygrid,nsegment !computational toroidal region is 1/nsegment of the full torus
    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(in):: s_dft(0:m-1,n)
    real(p_),intent(out):: s(m,n)
    real(p_):: a(n),b(n),tor(m),cos_tor(m),sin_tor(m)
    integer:: i,j,i_positive,i_negative
    !complex(p_):: ii=(0._p_,1._p_)
    i_positive=1 !harmonic that needs to be kept in the computational toroidal segment
    ! i_negative=m-i_positive
    ! write(*,*) 'i_positive=',i_positive,'i_negative=',i_negative
    !s=0._p_

    a(:)=2*real(s_dft(i_positive,:))/m
    b(:)=-2*imag(s_dft(i_positive,:))/m

    tor(:)=i_positive*twopi*ygrid(1:m)/toroidal_range
    cos_tor(:)=cos(tor(:))
    sin_tor(:)=sin(tor(:))

    !$omp parallel do
    do j=1,n
       do i=1,m
          s(i,j)=a(j)*cos_tor(i)+b(j)*sin_tor(i)
       enddo
    enddo
    !$omp end parallel do
  end subroutine toroidal_reconstruct


  subroutine radial_fourier_filter(field,m,n)
    use constants,only:p_

    implicit none
    integer,intent(in):: m,n
    complex(p_),intent(inout)::   field(0:m-1,0:n-1)
    complex(p_):: out(0:m-1,0:n-1)
    integer:: i,j,j_negative_wn

    out=(0._p_,0._p_) !initialized to zero
    do i=0,m-1
       do j=0,n/2  !scanning over the positive wavenumber
          if(j.le.radial_harmonics_included) then 
             out(i,j)=field(i,j) !positive wavenumber
             j_negative_wn=n-j
             if(j.eq.0) j_negative_wn=0
             out(i,j_negative_wn)=field(i,j_negative_wn) !negative or zero wavenumber
          endif
       enddo
    enddo

    field=out
  end subroutine radial_fourier_filter

  subroutine radial_sine_filter_core(s,m,n)
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(inout):: s(0:m-1,0:n-1)
    integer:: i,j
    do j=0,n-1
       if(j.gt.radial_harmonics_included) s(:,j)=0._p_
    enddo
  end subroutine radial_sine_filter_core

  subroutine radial_sine_high_pass(s, m, n)
    use constants, only: p_
    implicit none
    integer, intent(in) :: m, n
    real(p_), intent(inout) :: s(0:m-1, 0:n-1)
    integer, parameter :: jcut = 2
    integer :: i, j
    do j = 0, n-1
       if(j < jcut) s(:, j) = 0._p_
    enddo
  end subroutine radial_sine_high_pass

  
!!$subroutine radial_sine_filter_em_field()
!!$  use constants,only:p_
!!$  use perturbation_field,only: ex=>ex_left,ey=>ey_left,epar=>epar_left !as input and output
!!$  use perturbation_field,only: mx=>mf_x_left,my=>mf_y_left,mpar=>mf_par_left !as input and output
!!$  use magnetic_coordinates,only: m=>mtor,n=>nrad
!!$  use transform_module
!!$  implicit none
!!$  real(p_)::    epar_dst(m+1,n), ex_dst(m+1,n),  ey_dst(m+1,n)
!!$  real(p_)::    mpar_dst(m+1,n), mx_dst(m+1,n),  my_dst(m+1,n)
!!$
!!$  call oned_sine_transform2(ex,ex_dst,m+1,n) 
!!$  call oned_sine_transform2(ey,ey_dst,m+1,n) 
!!$  call oned_sine_transform2(epar,epar_dst,m+1,n) 
!!$  call radial_sine_filter_core(ex_dst,m+1,n)
!!$  call radial_sine_filter_core(ey_dst,m+1,n)
!!$  call radial_sine_filter_core(epar_dst,m+1,n)
!!$  call oned_inverse_sine_transform2(ex_dst,ex,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(ey_dst,ey,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(epar_dst,epar,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$
!!$  call oned_sine_transform2(mx,mx_dst,m+1,n) 
!!$  call oned_sine_transform2(my,my_dst,m+1,n) 
!!$  call oned_sine_transform2(mpar,mpar_dst,m+1,n) 
!!$  call radial_sine_filter_core(mx_dst,m+1,n)
!!$  call radial_sine_filter_core(my_dst,m+1,n)
!!$  call radial_sine_filter_core(mpar_dst,m+1,n)
!!$  call oned_inverse_sine_transform2(mx_dst,mx,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(my_dst,my,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$  call oned_inverse_sine_transform2(mpar_dst,mpar,m+1,n) !computing 1d inverse DST of s(:,:) along the second dimension
!!$
!!$
!!$end subroutine radial_sine_filter_em_field


  subroutine radial_sine_filter(s)
    use constants,only:p_
    use transform_module
    implicit none
    real(p_),intent(inout):: s(:,:)
    real(p_):: s_dst(size(s,1),size(s,2))
    integer:: m,n

    m=size(s,1)
    n=size(s,2)
    call oned_sine_transform2(s,s_dst,m,n) 
    call radial_sine_filter_core(s_dst,m,n)
    call oned_inverse_sine_transform2(s_dst,s,m,n) !computing 1d inverse DST of s(:,:) along the second dimension
  end subroutine radial_sine_filter


  subroutine radial_sine_reconstruct(s_dst,s,m,n)
    use constants,only:p_
    use constants,only: pi
    use magnetic_coordinates,only: xgrid
    implicit none
    integer,intent(in)::m,n
    real(p_),intent(in):: s_dst(m,0:n-1)
    real(p_),intent(out)::s(m,0:n-1)
    integer:: j,k
    real(p_):: radial_range,x,dx

    !  dx=xgrid(2)-xgrid(1)
    ! radial_range=dx*(n+1)
    s=0._p_
    do j=0,n-1
       !   x=dx*(j+1)
       do k=0,radial_harmonics_included
          s(:,j)=s(:,j)+s_dst(:,k)*sin((k+1)*pi*(j+1)/(n+1)) !see the formula in my notes on Fourier analysis
          !s(:,j)=s(:,j)+s_dst(:,k)*sin((k+1)*pi*x/radial_range) !see the formula in my notes on Fourier analysis
       enddo
    enddo
    s=s/(n+1)
  end subroutine radial_sine_reconstruct

  subroutine toroidal_fourier_radial_sine_filter(s,m,n)
    use transform_module,only: dst_dft
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(inout):: s(0:m-1,0:n-1)
    complex(p_):: s_spectrum(0:m-1,0:n-1)
    real(p_):: s_tmp(0:m-1,0:n-1)

    call dst_dft(s,s_spectrum,m,n)
    call toroidal_reconstruct(s_spectrum,s_tmp,m,n)
    call radial_sine_reconstruct(s_tmp,s,m,n)

  end subroutine toroidal_fourier_radial_sine_filter

  subroutine filter_basic(field)
    use constants,only:p_
    use transform_module
    implicit none
    real(p_),intent(inout):: field(:,:)
    complex(p_):: field_dft(size(field,1),size(field,2))
    real(p_):: field_dst(size(field,1),size(field,2))
    integer:: m,n
    m=size(field,1)
    n=size(field,2)

    call oned_fourier_transform1(field,field_dft,m,n) 
    call toroidal_reconstruct(field_dft,field,m,n)

    call oned_sine_transform2(field,field_dst,m,n) 
    call radial_sine_filter_core(field_dst,m,n)
    call oned_inverse_sine_transform2(field_dst,field,m,n) !computing 1d inverse DST of s(:,:) along the second dimension

  end subroutine filter_basic


  subroutine toroidal_fourier_radial_fourier_filter(s,m,n)
    !use fourn_module,only: twod_fourier_transform_nr,twod_inverse_fourier_transform_nr
    use transform_module
    use constants,only:p_
    implicit none
    integer,intent(in):: m,n
    real(p_),intent(inout):: s(0:m-1,0:n-1)
    complex(p_):: s_spectrum(0:m-1,0:n-1)

    !  call twod_fourier_transform_nr(s,s_spectrum,m,n)
    call twod_fourier_transform(s,s_spectrum,m,n)
    call toroidal_filter(s_spectrum,m,n)
    call radial_fourier_filter(s_spectrum,m,n)
    !  call twod_inverse_fourier_transform_nr(s_spectrum,s,m,n) 
    call twod_inverse_fourier_transform(s_spectrum,s,m,n) 

  end subroutine toroidal_fourier_radial_fourier_filter

end module filter_module


module smoothing_module
  use constants, only: p_, one, two, six
  use communication_connection, only: get_nearby_field_along_field_line
  implicit none
  real(p_), parameter :: weight1 = 0.5_p_, weight2 = -one/six !got to know these values from the GEM code

contains

  subroutine smoothing_along_field_line_core(a)
    real(p_), intent(inout) :: a(:,:)
    real(p_) :: a_left(size(a,1), size(a,2)), a_right(size(a,1),size(a,2))

    call get_nearby_field_along_field_line(a, a_left, a_right) !get value of field on the two grids that are to the left/rightof the present grid

    a=(weight1*a_right + a + weight1*a_left)/(one+two*weight1) !smoothing using weight1

    !the smoothing can be appled multiple times (possibley with different weights):
    call get_nearby_field_along_field_line(a, a_left, a_right)

    a=(weight2*a_right + a +weight2*a_left)/(one+two*weight2) !smoothing using weight2

  end subroutine smoothing_along_field_line_core


  subroutine smoothing_along_field_line_core5(a)
    real(p_), intent(inout) :: a(:,:)
    real(p_) :: a_left(size(a,1),size(a,2)), a_right(size(a,1),size(a,2))
    real(p_) :: a_left2(size(a,1),size(a,2)), a_right2(size(a,1),size(a,2))

    call get_nearby_field_along_field_line(a, a_left, a_right, a_left2, a_right2) 

    a=(-a_left2 + 4*a_left + 10*a + 4*a_right -a_right2)/16.

  end subroutine smoothing_along_field_line_core5

end module smoothing_module
module derivatives_in_xyz
  use constants,only : p_, two
  implicit none

contains
  
subroutine x_derivative0(field, field_x)
    use magnetic_coordinates, only: dradcor
    complex(p_), intent(in) :: field(:,:)
    complex(p_), intent(out) :: field_x(:,:)
    integer :: m, n, i, j

    m=size(field,1)
    n=size(field,2)
    do i = 1, m
       do j = 2, n-1
          field_x(i,j) = (field(i,j+1)-field(i,j-1))/(two*dradcor)
       enddo
       field_x(i,1) = two*field_x(i,2)-field_x(i,3) !linear interpolation to obtain the derivative at the boundary point
       field_x(i,n) = two*field_x(i,n-1)-field_x(i,n-2) !linear interpolation
    enddo
  end subroutine x_derivative0
  
  subroutine x_derivative(field, field_x)
    !calculate the derivative of field with respect to psi in field-line-following-coordinates (psi,theta,alpha)
    use magnetic_coordinates, only: dradcor
    real(p_), intent(in) :: field(:,:)
    real(p_), intent(out) :: field_x(:,:)
    integer :: m, n, i, j

    m=size(field,1)
    n=size(field,2)
    do i=1,m
       do j=2,n-1
          field_x(i,j)= (field(i,j+1)-field(i,j-1))/(two*dradcor)
       enddo
!!$       field_x(i,1)=(field(i,2)-0._p_)/(two*dradcor) !zero boundary condition for the field is used
!!$       field_x(i,n)=(0._p_-field(i,n-1))/(two*dradcor)!zero boundary condition for the field is used
       field_x(i,1)=two*field_x(i,2)-field_x(i,3) !linear interpolation to obtain the derivative at the boundary point
       field_x(i,n)=two*field_x(i,n-1)-field_x(i,n-2) !linear interpolation
    enddo

  end subroutine x_derivative

  subroutine y_derivative(field, field_y)
    !partial derivative of field with respect to y
    use magnetic_coordinates,only: dtor
    real(p_), intent(in) :: field(:,:)
    real(p_), intent(out) :: field_y(:,:)
    integer :: m,n,i,j,i_left,i_right

    m=size(field,1)
    n=size(field,2)
    do i=1,m   
       do j=1,n
          i_left=i-1
          if(i.eq.1) i_left=m !periodic boundary condition
          i_right=i+1
          if(i.eq.m) i_right=1 !periodic boundary condition
          field_y(i,j)=(field(i_right,j)-field(i_left,j))/(two*dtor)
       enddo
    enddo

  end subroutine y_derivative

  subroutine z_derivative(a, a_theta) ! calculating derivative along a magnetic field line
    use domain_decomposition, only: dtheta2, myid, my_right, my_left, Tube_COMM, gclr
    use magnetic_coordinates, only: mpol2
    use communication_connection, only: connect_condition_across_theta_cut
    use mpi
    implicit none
    real(p_), intent(in) :: a(:,:,:)
    real(p_), intent(out) :: a_theta( size(a,1), size(a,2) )
    integer :: m, n, status(MPI_STATUS_SIZE), ierr
    real(p_) :: a_left(size(a,1),size(a,2))

    m = size(a,1)
    n = size(a,2)

    call MPI_Sendrecv(a(:,:,1), m*n, MPI_real8, my_right, 1, &
         &            a_left,   m*n, MPI_real8, my_left,  1, Tube_COMM, status, ierr)

    if(GCLR == 0) call connect_condition_across_theta_cut(a_left, -1)
    a_theta(:,:) = (a(:,:,2) - a_left(:,:))/(two*dtheta2) !centered difference
  end subroutine z_derivative

  subroutine z_derivative0(a, a_theta) ! calculating derivative along a magnetic field line
    use domain_decomposition, only: dtheta2
    implicit none
    real(p_), intent(in) :: a(:, :, :)
    real(p_), intent(out) :: a_theta(:, :)

    a_theta(:,:) = (a(:,:,2) - a(:,:,1))/dtheta2 ! forward difference
  end subroutine z_derivative0

end module derivatives_in_xyz
module gk_polarization
  use constants,only: p_
  implicit none
!!$  real(p_),allocatable:: sigma(:,:)
!!$  complex(p_),dimension(:,:,:),allocatable::  u,  vt
!!$  complex(p_),dimension(:,:,:),allocatable::  ut, v
!!$  real(p_),parameter:: singular_value_threshold=6d-4 !singular value s smaller than this will be revmoved, (1/s be replace by zero in the inverse of sigma matrix)
contains
  subroutine prepare_polarization_matrix(ns, mmm) 
    use constants, only:zero,one,two,pi,twopi,kev,epsilon0, ii
    use control_parameters, only: space_charge_switch, nh_max
    use normalizing, only : tu,qu, nu
    use gk_module, only: mass_gk,charge_gk, gk_flr
    use magnetic_coordinates, only: nrad,mtor,toroidal_range, dradcor, dtheta,&
         & xgrid, xlow, xupp, grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use table_in_mc, only: b_mc
    use domain_decomposition, only: ipol_eq
    use math, only: srcbes, bessi0
    use gk_profile_funcs, only: gkt_func, gkn_func
    complex(p_), intent(out) :: mmm(:,:,0:)
    integer :: nx, j, jp, jeq, m, ns, kn
    real(p_) :: ky, kx1, kx2, kper1_sq, kper2_sq, b1, b2
    real(p_) :: radcor, bval, gx0, gy0, gxdgy0
    real(p_) :: lx, x, ly
    real(p_) :: omega, gyro_radius_sq, gamma1,gamma2, trash !omega, cycontron angular frequency
    complex(p_) :: part1, part2, sum, space_charge1, space_charge2
    real(p_) :: coefficient, debye_length_sq

    mmm = (0._p_,0._p_)
    if(gk_flr(ns) .eqv. .false.) return !no polarization density
    
    debye_length_sq=tu*kev*epsilon0/(nu*qu**2)
    !write(*,*) 'debye_length=',sqrt(debye_length_sq)
    nx = nrad - 1
    lx = xupp - xlow
    ly = toroidal_range

    do j = 1, nx-1
       jeq = j+1
       radcor = xgrid(jeq)
       x = radcor - xlow
       gx0 = grad_psi(ipol_eq, jeq)
       gy0 = grad_alpha(ipol_eq, jeq)
       gxdgy0 = grad_psi_dot_grad_alpha(ipol_eq, jeq)
       bval = b_mc(ipol_eq, jeq)
       omega = bval*abs(charge_gk(ns))/mass_gk(ns)
       gyro_radius_sq = (gkt_func(radcor,ns)*kev/mass_gk(ns))/omega**2
       coefficient = (charge_gk(ns)/qu)**2/(gkt_func(radcor,ns)/tu)*(gkn_func(radcor,ns)/nu)
       do kn = 0, nh_max !toroidal harmonics
          ky = kn*twopi/ly !the toroidal wavenumber
          do jp = 1, nx-1
             sum = 0._p_
             do m = 1, nx-1
                kx1 = m*pi/lx
                kx2 = -kx1
                kper1_sq = kx1**2*gx0**2 + two*kx1*ky*gxdgy0 + ky**2*gy0**2
                kper2_sq = kx2**2*gx0**2 + two*kx2*ky*gxdgy0 + ky**2*gy0**2
                b1 = kper1_sq*gyro_radius_sq
                b2 = kper2_sq*gyro_radius_sq
                call srcbes(b1, gamma1, trash)
                call srcbes(b2, gamma2, trash)
                part1 =  exp(ii*x*kx1)*(one-gamma1)
                part2 = -exp(ii*x*kx2)*(one-gamma2)
                space_charge1 =  exp(ii*x*kx1)*kper1_sq*debye_length_sq
                space_charge2 = -exp(ii*x*kx2)*kper2_sq*debye_length_sq
                !write(*,*) kper1_sq*debye_length_sq,kper2_sq*debye_length_sq
                !----use Pade approximation: gamma=I0(b)e^(-b)=1/(1+b). The results agree with the above more accurate treatment
!!$                part1= exp(ii*x*kx1)*(one-one/(one+b1)) 
!!$                part2=-exp(ii*x*kx2)*(one-one/(one+b2))
                !----pade approximation end----

                !--- bessi0(b1)*exp(-b1), sometime gave NaN, so not reliable, do not use it.
!!$                part1= exp( ii*m*pi*x/lx)*(one-bessi0(b1)*exp(-b1)) 
!!$                part2=-exp(-ii*m*pi*x/lx)*(one-bessi0(b2)*exp(-b2))
                !------
                sum = sum + sin(jp*m*pi/nx)*(coefficient*(part1+part2) &
                     & + space_charge_switch*(space_charge1+space_charge2))
             enddo
             mmm(j, jp, kn) = sum/(two*ii)*two/nx
          enddo
       enddo
    enddo

  end subroutine prepare_polarization_matrix

  subroutine prepare_polarization_matrix2(ns, mmm, nx) 
    ! direct numerical computation of the poloial density matrix,
    !refer to my note "nonlinear_gyrokinetic_equation.tm"
    use constants, only: zero,one,two,pi,twopi,kev,epsilon0, four, ii
    use normalizing, only: bn,ln,tu,qu, nu
    use gk_module, only: gk_flr, mass_gk, charge_gk
    use magnetic_coordinates, only: nrad, mtor, toroidal_range, dradcor, zgrid, dtheta, &
         & xgrid,  r_mc, z_mc, tor_shift_mc, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: bphi_mc
    use magnetic_coordinates, only: grad_alpha_r, grad_alpha_z, grad_psi_r, grad_psi_z
    use gk_profile_funcs, only : gkt_func, gkn_func
    use func_in_mc, only: tor_shift_func
    use control_parameters, only: nh_min, nh_max
    use magnetic_field, only: pfn_func
    use map_to_mc, only: interpolate_from_cylindrical_to_magnetic_coordinates1
    integer, intent(in) :: ns, nx
    complex(p_), intent(out) :: mmm(nx, nx, 0:nh_max)
    integer ::  j, jeq, jp0, jp1, ivper, j1, j2, kn
    real(p_) :: x, bphi
    real(p_) :: dvper,vper_range, vper, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg
    real(p_) :: coefficient, local_contribution, vt, xlow, xupp
    real(p_) :: R0, Z0, R, Z, dr, dz, dxdr, dxdz
    real(p_) :: dydr, dydz, delta_y, tor_shift0, angle1, angle2, c0, c1, xring, theta_tmp
    integer, parameter :: nvper = 100, n1 = 16, n2 = 32

    mmm=(0._p_, 0._p_) !mmm matrix corresponds to -np/nu, where np is the poloarization number density
    if(gk_flr(ns) .eqv. .false.) return

    xlow = xgrid(2) 
    xupp = xgrid(nrad-1)
    vper_range = 3._p_ !in unit of local sqrt(T/m)
    dvper = vper_range/(nvper-1) 

    do j = 1, nx
       jeq = j + 1
       x = xgrid(jeq)
       bphi = bphi_mc(ipol_eq, jeq)
       omega = abs(bphi*charge_gk(ns))/mass_gk(ns)
       coefficient = (charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)*(gkn_func(x,ns)/nu)
       R0 = r_mc(ipol_eq, jeq)
       Z0 = z_mc(ipol_eq, jeq)
       tor_shift0 = tor_shift_mc(ipol_eq, jeq)
       dydr = grad_alpha_r(ipol_eq, jeq)
       dydz = grad_alpha_z(ipol_eq, jeq)
       dxdr = grad_psi_r(ipol_eq, jeq)
       dxdz = grad_psi_z(ipol_eq, jeq)
       vt = sqrt(gkt_func(x,ns)*kev/mass_gk(ns))
       do ivper = 1, nvper
          vper = 0 + dvper*(ivper-1)
          gyro_radius = vper*vt/omega
          do j1 = 1, n1
             angle1 = 0. + twopi/(n1-1)*(j1-1)
             do j2 = 1, n2
                angle2 = 0.0 + twopi/(n2-1)*(j2-1)
                dr = gyro_radius*cos(angle1) + gyro_radius*cos(angle2)
                dz = gyro_radius*sin(angle1) + gyro_radius*sin(angle2)
                xring = pfn_func(R0+dr, Z0+dz)
                !xring = x + dxdr*dr + dxdz*dz
                if((xring >= xupp) .or. (xring < xlow)) cycle !no contribution
                jp0 = 1 + floor((xring-xlow)/dradcor)
                jp1 = jp0 + 1
                jeq = jp0 + 1 !eauilibrium xgrid index corresponding to jp0
                c0 = (xring - xgrid(jeq))/dradcor 
                c1 = one - c0
                !call interpolate_from_cylindrical_to_magnetic_coordinates1(R, Z, theta_tmp) !wrong at theta cut
                !delta_y = tor_shift0 - tor_shift_func(theta_tmp, xring) !wrong at theta cut
                delta_y = dydr*dr + dydz*dz
                do kn = nh_min, nh_max !toroidal harmonics
                   intg = dvper*vper*exp(-vper**2/two)*exp(ii*kn*twopi/toroidal_range*delta_y)
                   mmm(j, jp0, kn) = mmm(j, jp0, kn) + intg*c1
                   mmm(j, jp1, kn) = mmm(j, jp1, kn) + intg*c0
                enddo
             enddo
          enddo
       enddo
       mmm(j,:,:)= -coefficient*mmm(j,:,:)/(n1*n2)       
    enddo

    do j = 1, nx
       jeq = j + 1
       x = xgrid(jeq)
       local_contribution = (gkn_func(x,ns)/nu)*(charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)
       mmm(j, j, :) = mmm(j, j, :) + local_contribution !the non-gyro-averaging contribution
    enddo

  end subroutine prepare_polarization_matrix2


  subroutine prepare_polarization_matrix3(ns, mmm) 
    ! direct numerical computation of the poloial density matrix,
    !refer to my note "nonlinear_gyrokinetic_equation.tm"
    use constants, only: zero,one,two,pi,twopi,kev,epsilon0, four, ii
    use normalizing, only: bn,ln,tu,qu, nu
    use gk_module, only: gk_flr, mass_gk, charge_gk
    use magnetic_coordinates, only: nrad, mtor, toroidal_range, dradcor, zgrid, dtheta, &
         & xgrid,  r_mc, z_mc, tor_shift_mc, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: bphi_mc
    use magnetic_coordinates, only: xlow, xupp, grad_alpha_r, grad_alpha_z, grad_psi_r, grad_psi_z
    use gk_profile_funcs, only : gkt_func, gkn_func
    use func_in_mc, only: tor_shift_func
    use control_parameters, only: nh_min, nh_max
    use magnetic_field, only: pfn_func
    use map_to_mc, only: interpolate_from_cylindrical_to_magnetic_coordinates1
    integer, intent(in) :: ns
    complex(p_), intent(out) :: mmm(:,: , 0:)
    integer ::  j, jeq, m, ivper, j1, j2, kn, k, nx
    real(p_) :: x, bphi
    real(p_) :: dvper, vper_range, vper, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg, sum
    real(p_) :: coefficient, local_contribution, vt, lx, ly
    real(p_) :: R0, Z0, R, Z, dr, dz, dxdr, dxdz
    real(p_) :: dydr, dydz, delta_y, tor_shift0, angle1, angle2, xring
    integer, parameter :: nvper = 50, n1 = 8, n2 = 8

    mmm=(0._p_, 0._p_) !mmm matrix corresponds to -np/nu, where np is the poloarization number density
    if(gk_flr(ns) .eqv. .false.) return

    nx = nrad -1
    lx = xupp - xlow
    ly = toroidal_range
    vper_range = 3._p_ !in unit of local sqrt(T/m)
    dvper = vper_range/(nvper-1) 
    do j = 1, nx-1
       jeq = j + 1
       x = xgrid(jeq)
       bphi = bphi_mc(ipol_eq, jeq)
       omega = abs(bphi*charge_gk(ns))/mass_gk(ns)
       coefficient = (charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)*(gkn_func(x,ns)/nu)
       R0 = r_mc(ipol_eq, jeq)
       Z0 = z_mc(ipol_eq, jeq)
       tor_shift0 = tor_shift_mc(ipol_eq, jeq)
       dydr = grad_alpha_r(ipol_eq, jeq)
       dydz = grad_alpha_z(ipol_eq, jeq)
       dxdr = grad_psi_r(ipol_eq, jeq)
       dxdz = grad_psi_z(ipol_eq, jeq)
       vt = sqrt(gkt_func(x,ns)*kev/mass_gk(ns))
       do ivper = 1, nvper
          vper = 0 + dvper*(ivper-1)
          gyro_radius = vper*vt/omega
          do j1 = 1, n1
             angle1 = 0. + twopi/(n1-1)*(j1-1)
             do j2 = 1, n2
                angle2 = 0. + twopi/(n2-1)*(j2-1)
                dr = gyro_radius*cos(angle1) + gyro_radius*cos(angle2)
                dz = gyro_radius*sin(angle1) + gyro_radius*sin(angle2)
                !xring = pfn_func(R0+dr, Z0+dz) - xlow
                xring = x + dxdr*dr + dxdz*dz - xlow
                if((xring > lx) .or. (xring < 0)) cycle !no contribution
                delta_y = dydr*dr + dydz*dz
                do kn = nh_min, nh_max !toroidal harmonics
                   intg = dvper*vper*exp(-vper**2/two)*exp(ii*kn*twopi/ly*delta_y)
                   do k = 1, nx-1
                      sum = 0
                      do m = 1, nx-1
                         sum = sum + intg*2/nx*sin(k*m*pi/nx)*sin(m*pi/lx*xring)
                      enddo
                      mmm(j, k, kn) = mmm(j, k, kn) + sum
                   enddo
                enddo
             enddo
          enddo
       enddo
       mmm(j,:,:) = -coefficient*mmm(j,:,:)/(n1*n2)       
    enddo

    do j = 1, nx-1
       jeq = j + 1
       x = xgrid(jeq)
       local_contribution = (gkn_func(x,ns)/nu)*(charge_gk(ns)/qu)**2/(gkt_func(x,ns)/tu)
       mmm(j, j, :) = mmm(j, j, :) + local_contribution !the non-gyro-averaging contribution
    enddo

  end subroutine prepare_polarization_matrix3
  

  pure  subroutine prepare_slowing_down_polarization_matrix(mmm, mass, charge, e_cut) 
    ! direct numerical computation of the poloial density matrix for slowing-down distribution
    use constants,only:zero,one,two,pi,twopi, Mev,epsilon0, four, kev
    use normalizing,only:bn,ln,tu,qu, nu
    use magnetic_coordinates,only: nrad,mtor,toroidal_range, dradcor,zgrid,dtheta,&
         & xgrid, xlow, xupp, r_mc, z_mc, &
         & grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use domain_decomposition,only: ipol_eq
    use table_in_mc,only: b_mc
    use func_in_mc, only : tor_shift_func
    use magnetic_field, only : radcor_as_func_of_pfn
    use control_parameters, only : nh_max
    use magnetic_field, only : pfn_func
    use map_to_mc, only : interpolate_from_cylindrical_to_magnetic_coordinates1

    complex(p_), intent(out) :: mmm(:,:,0:)
    real(p_), intent(in) :: mass, charge, e_cut
    complex(p_), parameter :: ii = (0._p_, 1._p_)
    integer ::  nx, jp, j, jeq, m, ierr, ivper, ivpar, j1, j2, ns
    integer :: kn
    real(p_) :: ky, x, bval, gx0, energy, v, v_cut
    real(p_) :: dvper, dvpar, vper, vpar, omega, gyro_radius !omega, cycontron angular frequency
    complex(p_) :: intg
    real(p_) :: local_contribution
    real(p_) :: R0, Z0, R, Z, delta_y, tor_shift0, angle1, angle2, c0, c1, x_tmp, theta_tmp
    integer, parameter :: nvper=50, nvpar=50, n1=8, n2=8

    nx = nrad - 2 !perturbations at the two end points are fixed at zero.
    !allocate(mmm(nx,nx,0:nh_max)) !the matrix corresponds to -np/nu, where np is the poloarization number density

    v_cut = sqrt(2*e_cut/mass)
    dvper = v_cut/(nvper-1)
    dvpar = 2*v_cut/(nvpar-1)

    mmm = (0._p_, 0._p_)
    do j = 1, nx
       jeq = j+1
       x = xgrid(jeq)
       gx0 = grad_psi(ipol_eq,jeq)
       bval = b_mc(ipol_eq,jeq)
       omega = bval*abs(charge)/mass
       R0 = r_mc(ipol_eq, jeq)
       Z0 = z_mc(ipol_eq, jeq)
       tor_shift0 = tor_shift_func(zgrid(ipol_eq), x)

       do ivpar = 1, nvpar
          vpar =  -v_cut + dvpar*(ivpar-1)
          do ivper = 1, nvper
             vper = 0 + dvper*(ivper-1)
             v = sqrt(vpar**2 + vper**2)
             energy = 0.5*mass*v**2
             if(v > v_cut) cycle
             gyro_radius = vper/omega
             do j1 = 1, n1
                angle1 = 0.+twopi/(n1-1)*(j1-1)
                do j2 = 1, n2
                   angle2 = 0.+twopi/(n2-1)*(j2-1)
                   R = R0+gyro_radius*cos(angle1)+gyro_radius*cos(angle2)
                   Z = Z0+gyro_radius*sin(angle1)+gyro_radius*sin(angle2)
                   x_tmp = radcor_as_func_of_pfn(pfn_func(R,Z))
                   if((x_tmp.ge.xupp) .or. (x_tmp.le.xlow)) cycle !no contribution to the polarization
                   call interpolate_from_cylindrical_to_magnetic_coordinates1(R,Z,theta_tmp)
                   delta_y = tor_shift0 - tor_shift_func(theta_tmp, x_tmp)

                   do kn = 0, nh_max !toroidal harmonics
                      ky = twopi*kn/toroidal_range !toroidal mode number
                      !intg = dvper*vper*exp(-vper**2/two)*exp(ii*ky*delta_y)
                      intg = dvpar * dvper * vper * slowing_down_Ed(x, energy) *exp(ii*ky*delta_y)
                      jp = int((x_tmp-xlow)/dradcor) !locating
                      c1=(x_tmp-xgrid(jp+1))/dradcor !linear interpolating coefficient
                      c0=one-c1
                      if(jp == 0) then
                         mmm(j, jp+1, kn)= mmm(j, jp+1, kn) + intg*c1
                      elseif(jp == nx) then
                         mmm(j, jp, kn)= mmm(j,jp,kn) + intg*c0
                      elseif( (jp .gt. 0) .and. (jp .lt. nx)) then
                         mmm(j,jp,  kn)= mmm(j, jp, kn) + intg*c0   
                         mmm(j,jp+1,kn)= mmm(j, jp+1, kn) + intg*c1
                      else
                         !do nothing
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
       mmm(j,:,:)= - mmm(j,:,:)*twopi/(n1*n2)*(charge/qu)**2*tu*kev/nu
    enddo

    do j = 1, nx
       jeq = j+1
       x = xgrid(jeq)
       local_contribution = (charge/qu)**2*(f_Ed_integral(x))*tu*kev/nu
       mmm(j,j, :) = mmm(j,j, :) + 1*local_contribution !add the gyro-averaging-free part to the diagonal elements
    enddo

  end subroutine prepare_slowing_down_polarization_matrix

  
  pure real(p_) function slowing_down(x, E) result(f) !not used
    use constants,only: p_, kev
    use gk_radial_profiles, only : nalpha_object, alpha_normc, alpha_ecrit
    real(p_), intent(in) :: x, E
    real(p_) :: N, Ec

    N = alpha_normc%func(x)
    Ec = alpha_ecrit%func(x)*kev
    f = N/(1+sqrt(E/Ec)**3)

  end function slowing_down

  elemental real(p_) function slowing_down_Ed(x, E) result(f) 
    use constants,only: p_, kev
    use gk_radial_profiles, only : nalpha_object, alpha_normc, alpha_ecrit
    real(p_), intent(in) :: x, E
    real(p_) :: N, Ec

    N = alpha_normc%func(x)
    Ec = alpha_ecrit%func(x)*kev

    f = -N/(1+sqrt(E/Ec)**3)**2*1.5*sqrt(E/Ec)/Ec

  end function slowing_down_Ed

  pure real(p_) function f_Ed_integral(x) result(z)
    use constants,only: p_, kev, Mev, fourpi, atom_mass_unit
    implicit none
    real(p_), intent(in) :: x
    integer, parameter :: n = 100
    real(p_), parameter :: mass = 4*atom_mass_unit
    real(p_) :: vmax, v(n), E(n), dv
    integer :: i

    vmax = sqrt(2*3.5*Mev/mass)
    dv = vmax/(n-1)
    do i = 1, n
       v(i) = 0 + dv*(i-1)
    enddo
    E = 0.5*mass*v**2
    z =sum(slowing_down_Ed(x,E)*v**2)*fourpi*dv

  end function f_Ed_integral


  pure real(p_) function dNdx(x) result (z) 
    use gk_radial_profiles, only : alpha_normc
    real(p_), intent(in) :: x
    real(p_), parameter :: dx = 1d-3

    z = (alpha_normc%func(x+dx) - alpha_normc%func(x-dx))/(2*dx)

  end function dNdx

  pure real(p_) function dEcdx(x) result (z) !SI
    use constants,only: p_, kev
    use gk_radial_profiles, only : alpha_ecrit
    real(p_), intent(in) :: x
    real(p_), parameter :: dx = 1d-3

    z = (alpha_ecrit%func(x+dx) - alpha_ecrit%func(x-dx))*kev/(2*dx)

  end function dEcdx


  pure real(p_) function slowing_down_xd(x, E) result(f) 
    use constants,only: p_, kev
    use gk_radial_profiles, only : alpha_normc, alpha_ecrit
    real(p_), intent(in) :: x, E
    real(p_) :: N, Ec

    N = alpha_normc%func(x)
    Ec = alpha_ecrit%func(x)*kev

    f = dNdx(x)/(1+sqrt(E/Ec)**3) - N/((1+sqrt(E/Ec)**3))**2*1.5*sqrt(E/Ec)*(-E/Ec**2*dEcdx(x))

  end function slowing_down_xd


end module gk_polarization
module gyro_ring_mod
  use constants, only : p_, two, twopi
  
  real(p_), allocatable :: cos_gyro(:,:), sin_gyro(:,:)

contains
  subroutine set_gyro_phase()
    use domain_decomposition, only: myid
    use gk_module, only: gyro_npt, nmmax
    use math, only: random_yj
    implicit none
    integer :: i, j, iseed
    real(p_) :: gyro_angle0, gyro_angle, tmp

    allocate(cos_gyro(gyro_npt,nmmax))
    allocate(sin_gyro(gyro_npt,nmmax))

    iseed = -(1177+myid*3) !set the iseed in different procs
    tmp = random_yj(iseed) !just to trigger the use of the iseed, the generated random number is not used, 
    do i=1, nmmax
       gyro_angle0 = random_yj(0) !0 means using last random number as iseed, first gyro-angle
       gyro_angle0 = gyro_angle0*twopi !scaled to [0:twopi]
       !gyro_angle0=0 !comment in/out this if we do-not/do want the first gyro-angle to be random
       do j = 1, gyro_npt ! gyro-angles are uniform distributed relatively to the first one.
          gyro_angle =  gyro_angle0 + twopi/gyro_npt*(j-1)
          cos_gyro(j,i) = cos(gyro_angle)
          sin_gyro(j,i) = sin(gyro_angle)
       enddo
    enddo
  end subroutine set_gyro_phase

  subroutine gyro_ring(ns, lost, radcor, alpha, theta, mu, x_ring, y_ring, z_ring)
    !output are used by gyro_average and deposit_gk
    use gk_module, only:  gk_flr, nm_gk
    use magnetic_coordinates, only: xlow, xupp
    implicit none
    integer, intent(in) :: ns
    logical, intent(in) :: lost(:)
    real(p_), intent(in) :: theta(:), radcor(:), alpha(:), mu(:)
    real(p_), intent(out) :: x_ring(:, :), y_ring(:, :), z_ring(:, :)
    integer :: nm, k

    nm= nm_gk(ns)

    if(gk_flr(ns) .eqv. .false.) then
       do k = 1, nm
          if (lost(k) .eqv. .true.) cycle
          x_ring(1,k) = radcor(k)
          y_ring(1,k) = alpha(k)
          z_ring(1,k) = theta(k)
       enddo
    else
       do k =1, nm
          if (lost(k) .eqv. .true.) cycle
          call gyro_ring_core(ns, k, radcor(k), alpha(k), theta(k), mu(k), x_ring(:,k), y_ring(:,k), z_ring(:,k))
       enddo
    endif
  end subroutine gyro_ring


pure subroutine gyro_ring_core(ns, k, x0, y0, z0, mu, x, y, z)
    use magnetic_coordinates, only: mpol, nrad, zgrid, xgrid, dtheta, dradcor, &
         & grad_psi, grad_alpha, grad_psi_dot_grad_alpha, grad_psi_dot_grad_theta, &
         & pfn_inner, pfn_bdry, toroidal_range, GSpsi_prime
    use table_in_mc, only : bdgxcgz, b_mc
    use math,only: shift_toroidal
    use gk_module, only:   mass_gk, charge_gk, vn_gk, gyro_npt
    use domain_decomposition, only: theta_start, dtheta2
    use math, only: shift_toroidal
    use interpolate_module, only : linear_2d_interpolate0, locate
    implicit none
    integer, intent(in) :: ns, k
    real(p_), intent(in) :: x0, y0, z0, mu
    real(p_), intent(out) :: x(:), y(:), z(:)
    real(p_) :: bval0, gx0, gxdgy0, gxdgz0, bdgxcgz0
    real(p_) :: vper, gyro_radius, dy1, dy2, dz1, dz2
    integer :: i, j, kp

    call locate(mpol, zgrid, dtheta, z0, i)
    call locate(nrad, xgrid, dradcor,x0, j)
    call linear_2d_interpolate0(mpol,nrad,zgrid,xgrid,dtheta, dradcor, grad_psi,&
         & z0, x0,i,j, gx0)
    call linear_2d_interpolate0(mpol,nrad,zgrid,xgrid,dtheta, dradcor, grad_psi_dot_grad_alpha,&
         & z0, x0,i,j,gxdgy0)
    call linear_2d_interpolate0(mpol,nrad,zgrid,xgrid,dtheta, dradcor, grad_psi_dot_grad_theta,&
         & z0, x0,i,j, gxdgz0)
    call linear_2d_interpolate0(mpol, nrad, zgrid, xgrid, dtheta, dradcor, bdgxcgz,&
         & z0, x0,i,j, bdgxcgz0)
    call linear_2d_interpolate0(mpol, nrad, zgrid, xgrid, dtheta, dradcor, b_mc,&
         & z0, x0,i,j,bval0)

    vper=sqrt(two*bval0*mu)
    gyro_radius=vper*vn_gk(ns)/(bval0*abs(charge_gk(ns))/mass_gk(ns))

    dy1 = gyro_radius*gxdgy0/gx0
    dy2 = gyro_radius*bval0/(GSpsi_prime*gx0)
    dz1 = gyro_radius*gxdgz0/gx0
    dz2 = gyro_radius*bdgxcgz0/(bval0*gx0)

    do kp = 1, gyro_npt
       x(kp)= x0 + gyro_radius*cos_gyro(kp,k)*gx0
       x(kp) = max(x(kp), pfn_inner) !to prevent exceeding the radial computational box
       x(kp) = min(x(kp), pfn_bdry) 
       y(kp) = y0 + cos_gyro(kp,k)*dy1 + sin_gyro(kp,k)*dy2
       call shift_toroidal(y(kp), toroidal_range)
       z(kp) = z0 + cos_gyro(kp,k)*dz1 + sin_gyro(kp,k)*dz2
       z(kp) = max(z(kp), theta_start)
       z(kp) = min(z(kp), theta_start+dtheta2)
    enddo
  end subroutine gyro_ring_core
end module gyro_ring_mod
module gyro_average_mod
contains
pure  subroutine gyro_average0(flr, x_ring, y_ring, z_ring, phix, phix_av)
    use constants, only:  p_
    use gk_module, only: gyro_npt
    implicit none
    real(p_), intent(in) :: x_ring(:), y_ring(:), z_ring(:)
    logical, intent(in)  :: flr
    real(p_), intent(in) :: phix(:,:,:)
    real(p_), intent(out) :: phix_av
    real(p_) :: phixp(gyro_npt)
    integer :: kr, npt
    
    if (flr .eqv. .false.) then
       npt = 1
    else
       npt = gyro_npt
    endif

    do kr = 1, npt
       call field_at_particle0(x_ring(kr), y_ring(kr), z_ring(kr), phix, phixp(kr))
    enddo
    phix_av =sum(phixp(1:npt))/npt

  end subroutine gyro_average0

  pure subroutine field_at_particle0(radcor, alpha, theta, phix, phixp)
    use constants, only: p_, one, zero, twopi
    use magnetic_coordinates, only: mtor, nrad, ygrid, xgrid
    use domain_decomposition, only: dtheta2, theta_start
    implicit none
    real(p_), intent(in) :: radcor, theta, alpha, phix(:,:,:)
    real(p_), intent(out) :: phixp
    real(p_) :: c1, c2, tmp1, tmp2

    c1 = (theta-theta_start)/dtheta2
    c2 = one-c1
    tmp1 = interpolate2d(mtor, nrad, ygrid, xgrid, phix(:,:,1), alpha, radcor) 
    tmp2 = interpolate2d(mtor, nrad, ygrid, xgrid, phix(:,:,2), alpha, radcor) 
    phixp = tmp1*c2 + tmp2*c1

  end subroutine field_at_particle0


 pure subroutine gyro_average(ns, lost, x_ring, y_ring, z_ring, &
       & phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga)
    !output are used in computing drift and pushing weight.
    use constants, only:  p_
    use gk_module, only: gk_flr, gyro_npt, nm_gk
    use magnetic_coordinates, only: xlow, xupp
    implicit none
    integer, intent(in)  :: ns
    logical, intent(in)  :: lost(:)
    real(p_), intent(in) :: x_ring(:,:), y_ring(:,:), z_ring(:,:)
    real(p_), intent(out) :: phix_ga(:), phiy_ga(:), phiz_ga(:)
    real(p_), intent(out) :: ax_ga(:),ay_ga(:),az_ga(:), ahx_ga(:),ahy_ga(:),ahz_ga(:), ah_ga(:)
    real(p_), dimension(gyro_npt) :: phix0, phiy0, phiz0, ax0, ay0, az0, ahx0, ahy0, ahz0, ah0
    integer :: nm, k, kr, npt

    nm = nm_gk(ns)

    if(gk_flr(ns) .eqv. .false.) then
       npt=1
    else
       npt=gyro_npt
    endif

    do k = 1, nm
       if (lost(k) .eqv. .true.) cycle
       do kr = 1, npt
          call field_at_particle(x_ring(kr,k), y_ring(kr,k), z_ring(kr,k), & 
                              & phix0(kr), phiy0(kr), phiz0(kr), ax0(kr), ay0(kr), az0(kr), &
                              & ahx0(kr), ahy0(kr), ahz0(kr), ah0(kr))
       enddo
       phix_ga(k) =sum(phix0(1:npt))/npt
       phiy_ga(k) =sum(phiy0(1:npt))/npt
       phiz_ga(k) =sum(phiz0(1:npt))/npt
       ax_ga(k)   =sum(ax0(1:npt))/npt
       ay_ga(k)   =sum(ay0(1:npt))/npt
       az_ga(k)   =sum(az0(1:npt))/npt
       ahx_ga(k)  =sum(ahx0(1:npt))/npt
       ahy_ga(k)  =sum(ahy0(1:npt))/npt
       ahz_ga(k)  =sum(ahz0(1:npt))/npt
       ah_ga(k)   =sum(ah0(1:npt))/npt
    enddo

  end subroutine gyro_average


  pure subroutine field_at_particle(radcor, alpha, theta, phix0,phiy0,phiz0,ax0,ay0,az0,ahx0,ahy0,ahz0,ah0)
    use constants, only: p_, one, zero, twopi
    use perturbation_field, only: phix, phiy, phiz, ax, ay, az, ahx, ahy, ahz, apara_h
    use magnetic_coordinates, only: mtor, nrad, ygrid, xgrid
    use domain_decomposition, only: dtheta2, theta_start
    implicit none
    real(p_), intent(in) :: radcor, theta, alpha
    real(p_),intent(out) :: phix0, phiy0, phiz0, ax0, ay0, az0, ahx0, ahy0, ahz0, ah0
    real(p_) :: coeff1,coeff2,tmp1,tmp2

       coeff1=(theta-theta_start)/dtheta2
       coeff2=one-coeff1

       tmp1=interpolate2d(mtor, nrad, ygrid, xgrid, phix(:,:,1), alpha, radcor) 
       tmp2=interpolate2d(mtor, nrad, ygrid, xgrid, phix(:,:,2), alpha, radcor) 
       phix0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,phiy(:,:,1) ,alpha,radcor)  
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,phiy(:,:,2) ,alpha,radcor)  
       phiy0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,phiz(:,:,1) ,alpha,radcor)  
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,phiz(:,:,2) ,alpha,radcor)  
       phiz0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ax(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ax(:,:,2) ,alpha,radcor) 
       ax0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ay(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ay(:,:,2) ,alpha,radcor) 
       ay0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,az(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,az(:,:,2) ,alpha,radcor) 
       az0=tmp1*coeff2+tmp2*coeff1


       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ahx(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ahx(:,:,2) ,alpha,radcor) 
       ahx0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ahy(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ahy(:,:,2) ,alpha,radcor) 
       ahy0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,ahz(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,ahz(:,:,2) ,alpha,radcor) 
       ahz0=tmp1*coeff2+tmp2*coeff1

       tmp1=interpolate2d(mtor,nrad,ygrid,xgrid,apara_h(:,:,1) ,alpha,radcor) 
       tmp2=interpolate2d(mtor,nrad,ygrid,xgrid,apara_h(:,:,2) ,alpha,radcor) 
       ah0=tmp1*coeff2+tmp2*coeff1

  end subroutine field_at_particle


  pure function interpolate2d(nx, nz, xarray, zarray, psi,x,z) result(fval)
    use constants,only: p_, one
    implicit none
    integer, intent(in) :: nx, nz
    real(p_), intent(in) :: xarray(nx), zarray(nz), psi(nx,nz)
    real(p_), intent(in) :: x, z
    real(p_) :: fval
    real(p_) :: dx, dz, t1, t2, slope
    integer :: i, j, ii, jj, i_plus1

    dx = xarray(2) - xarray(1)
    i = floor(one+(x-xarray(1))/dx)
    i_plus1 = i+1
    if(i==nx) i_plus1 = 1 !periodic condition

    dz = zarray(2) - zarray(1)
    j = floor(one+(z-zarray(1))/dz)

    if(j>=nz .or. j<=1) then
       fval =0
       return
    endif

    slope = (psi(i_plus1,j)-psi(i,j))/dx
    t1 = psi(i,j)+slope*(x-xarray(i))
    slope = (psi(i_plus1,j+1)-psi(i,j+1))/dx
    t2 = psi(i,j+1)+slope*(x-xarray(i))
    slope = (t2-t1)/dz
    fval = t1+slope*(z-zarray(j))
  end function interpolate2d

end module gyro_average_mod
module poisson
  use constants, only : p_
  save
  complex(p_), allocatable :: poisson_matrix(:,:,:) 
  complex(p_), allocatable :: poisson_matrix0(:,:,:) !pure gk contribution, used in solving n=0 component when using adiabatic electrons, 
  integer, allocatable :: ipiv(:,:), ipiv0(:,:)

contains
  subroutine prepare_poisson_matrix()
    use constants, only: zero, elementary_charge, atom_mass_unit, Mev
    use magnetic_coordinates, only: nrad, xgrid, mtor
    use gk_radial_profiles, only: density_object, temperature_object
    use control_parameters, only: dt_omega_i_axis, nh_max, adiabatic_electrons
    use gk_module, only: nsm
    use gk_polarization
    use normalizing, only: tu,qu, nu
    use domain_decomposition, only: myid, numprocs, tclr, gclr
    use adiabatic_e_profiles, only: initialize_adiabatic_electron, ne_object, te_object
    implicit none
    complex(p_), allocatable :: polarization(:,:,:)
    complex(p_), allocatable :: mmm(:,:,:)
    real(p_) ::  x, tmp
    integer :: nx, i, j, jeq, info, kn, ns, u

    nx = nrad - 2
    allocate(polarization(nx, nx, 0:nh_max), source=(zero,zero))
    allocate(mmm(nx, nx, 0:nh_max))

    do ns = 1, nsm
       if(ns<30) then
          call prepare_polarization_matrix(ns, mmm)
          !call prepare_polarization_matrix2(ns, mmm, nx)
          !call prepare_polarization_matrix3(ns, mmm)
       else
          call prepare_slowing_down_polarization_matrix(mmm, 4*atom_mass_unit, 2*elementary_charge, 3.5*Mev)
       endif
       polarization(:,:,:) = polarization + mmm
    enddo
if(tclr==0) write(*,*) gclr, real(polarization(nx/2,nx/2+1, nh_max)), imag(polarization(nx/2,nx/2+1, nh_max))

    allocate(poisson_matrix(nx, nx, 0:nh_max))
    allocate(ipiv(nx, 0:nh_max))

    poisson_matrix = polarization

    if(adiabatic_electrons .eqv. .true.) then
       allocate(poisson_matrix0(nx, nx, 0:0))
       allocate(ipiv0(nx,0:0))
       poisson_matrix0 = polarization(:,:, 0:0)
       call initialize_adiabatic_electron()
       do j = 1, nx 
          jeq = j+1
          x = xgrid(jeq)
          tmp = (ne_object%func(x)/nu)*(elementary_charge/qu)**2/(te_object%func(x)/tu)
          poisson_matrix(j, j, :) = poisson_matrix(j, j, :) + tmp !add to the diagonal elemets of the matrix
       enddo
       kn = 0
       call ZGETRF(nx, nx, poisson_matrix0(:,:,kn), nx, ipiv0(:,kn), info) 
    endif

    do kn = 0, nh_max !for each toroidal Fourier component, LU factorize the radial coeff matrix
       ! Computes LU factorization of a general matrix using partial pivoting with row interchanges.
       call ZGETRF(nx, nx, poisson_matrix(:,:,kn), nx, ipiv(:,kn), info) 
    enddo

  end subroutine prepare_poisson_matrix


  subroutine solve_poisson(density_left, density_right, potential, phix, phiy, phiz)
    use constants, only: p_, zero, elementary_charge
    use magnetic_coordinates, only: m=>mtor, n=>nrad, jacobian_av, abs_jacobian, &
         & mpol2, xgrid, dradcor, dtor
    use control_parameters, only: fk_switch, filter_radial, dt_omega_i_axis, &
         & ismooth, nh_min, nh_max, adiabatic_electrons
    use gk_radial_profiles, only : density_object, temperature_object
    use normalizing,only: tu,qu, nu
    !use fk_module,only: mass_i,charge_i
    use gk_module, only : w_unit
    use domain_decomposition,only: myid, ipol_eq, tube_comm, multi_eq_cells, dtheta2, dvol
    use math, only : ZGETRS_wrapper, solver_toroidal_mode_number_parallel
    use filter_module
    use smoothing_module
    use derivatives_in_xyz
    use communication_connection
    use transform_module
    use mpi
    implicit none
    real(p_), intent(in) :: density_left(:,:), density_right(:,:)
    real(p_), intent(out), dimension(:,:,:) :: potential, phix, phiy, phiz
    real(p_) :: density(0:m-1,n), source(0:m-1,n-2), source_dst(0:m-1,n-2), phi(0:m-1,n-2)
    complex(p_) :: source_dft(0:m-1, n-2), phi_dft(0:m-1, n-2), phix_dft(0:m-1, n-2)
    real(p_) :: signal_dst(1:1, n-2), signal(1:1, n-2)
    real(p_) :: phi_n0(n-2), sum_phi(n-2), av_phi(n-2)
    integer :: kn, ierr, j,i, it, jeq
    real(p_) :: coeff, x

    call merge_source(density_left, density_right, density)
    do j=1,n
       density(:,j) = density(:,j)*w_unit/dvol(j)/nu
    enddo
    source = density(:,2:n-1)
    
    call oned_DFT_parallel_version(source, source_dft, m, n-2)

    if(nh_min == 0) then !remove harmonics of low kx to suppress numerical instabilities
       call oned_sine_transform2(real(source_dft(0:0,:)), signal_dst, 1, n-2)
       call radial_sine_high_pass(signal_dst, 1, n-2)
       call oned_inverse_sine_transform2(signal_dst, signal, 1, n-2)
       source_dft(0,:) = signal(1,:)
    endif
    
    phi_dft = 0 !only some toroidal harmonics will be solved, others are assumed to be zero

    if((adiabatic_electrons .eqv. .true.) .and. (nh_min==0)) then
       ! Solve for n=0 component while neglecting the electron contribution in the poisson_matrix.
       !The solution is incoorrect, but the surface averaging of it equals to the surface average of the correct solution
       kn=0; call ZGETRS_wrapper(kn, poisson_matrix0, IPIV0, source_dft(kn,:), phi_dft(kn,:))
       ! Convert to real type and include the Jacobian (to perform magnetic surface average):
       phi_n0(:) = real(phi_dft(0,:)) * abs_jacobian(ipol_eq, 2:n-1)

       call MPI_Allreduce(phi_n0, sum_phi, n-2, MPI_Double, MPI_sum, tube_comm, ierr)
       av_phi(:) = (sum_phi(:)/mpol2)/jacobian_av(2:n-1) !magnetic surface averaging
       do j = 1, n-2 !add <phi> to the right-hand side
          x = xgrid(j+1)
          coeff = (density_object(1)%func(x)/nu)*(elementary_charge/qu)**2/(temperature_object(1)%func(x)/tu)
          source_dft(0,j) = source_dft(0,j) + av_phi(j)*coeff
       enddo
    endif

    call solver_toroidal_mode_number_parallel(poisson_matrix, IPIV, source_dft, phi_dft)
    do kn = nh_min, nh_max
       !call mfilter_for_each_n(phi_dft(kn,:), n-2, kn)
    enddo
    do kn = 1, nh_max ! For negative toroidal mode number
       phi_dft(m-kn,:) = conjg(phi_dft(kn,:)) !use the Complex conjugate relation
    enddo
    !if(nh_min == 0) call mfilter_for_n0(phi_dft(0,:), n-2)
    if(nh_min == 0) call surface_average_of_n0(phi_dft(0,:), n-2)

    call oned_backward_DFT_parallel_version(phi_dft, phi, m, n-2) !radial task decomposion
    
    call x_derivative0(phi_dft, phix_dft) !more accurate than taking x derivative in real space
    call oned_backward_DFT_parallel_version(phix_dft, phix(:, 2:n-1, 1), m, n-2) 
    
    potential(:, 2:n-1, 1) = phi(:,:)
    potential(:, 1, 1) = 0._p_ !zero boundary condition
    potential(:, n, 1) = 0._p_  !zero boundary condition
    call update_scalar_at_right_boundary(potential)

    do j = 1, ismooth
       call smoothing_along_field_line_core(potential(:,:,1))
    enddo

    !call x_derivative(potential(:,:,1), phix(:,:,1))
    phix(:, 1, 1) = 0; phix(:, n, 1) = 0
    call y_derivative(potential(:,:,1), phiy(:,:,1))
    call z_derivative(potential(:,:,:), phiz(:,:,1)) 
    call update_derivatives_at_right_boundary(phix, phiy, phiz)

    !call electric_potential_to_field() !this version uses epar, ex, and ey rather than the cylindrical components.
    if(fk_switch==1) call potential_to_cylindrical_field() !this version use cylindrical components Er,Ephi,Ez

  end subroutine solve_poisson

end module poisson

subroutine potential_to_cylindrical_field()
  use mpi
  use constants,only: p_, one,two,one_half
  use magnetic_coordinates,only: mtor,nrad,dradcor,zgrid,dtheta,r_mc
  use magnetic_coordinates,only: grad_psi_r,grad_psi_z,grad_theta_r,grad_theta_z, &
       & grad_alpha_r,grad_alpha_z,grad_alpha_phi

  use perturbation_field,only: potential
  use domain_decomposition,only: theta_start,ipol_eq
  use perturbation_field,only: ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left !as output
  use control_parameters,only: filter_radial
  use constants,only:kev
  use normalizing,only:bn,ln, tu,qu
  use fk_module, only : vn_fk
  use domain_decomposition,only:myid
  use transform_module,only:oned_fourier_transform1, oned_backward_fourier_transform1,&
       &oned_sine_transform2,oned_inverse_sine_transform2
  use filter_module,only: toroidal_filter,radial_sine_filter_core
  use derivatives_in_xyz,only: x_derivative,y_derivative,z_derivative, z_derivative0
  use communication_connection, only: communicate_field_value_between_neighbour_cells2
  implicit none

  real(p_):: potential_psi(mtor,nrad),potential_alpha(mtor,nrad),potential_theta(mtor,nrad)
  real(p_):: ef_cyl_r(mtor,nrad),ef_cyl_z(mtor,nrad),ef_cyl_phi(mtor,nrad)
  !complex(p_):: ef_cyl_r_dft(mtor,nrad),ef_cyl_z_dft(mtor,nrad),ef_cyl_phi_dft(mtor,nrad)
  complex(p_):: potential_dft(mtor,nrad),out(mtor,nrad)
  real(p_):: potential_dst(mtor,nrad) !discrete sine transform
  integer:: i,j,i_left,i_right,j_left,j_right,ierr
  integer:: jeq
  real(p_):: normal
  real(p_):: tmp1(mtor,nrad),tmp2(mtor,nrad),tmp3(mtor,nrad)


  !smoothing is moved outside of this subroutine
!!$  if(filter_toroidal.eqv..true.) then !filter over the toroidal mode number, keeping the perturbation with desired toroidal mode number
!!$     call oned_fourier_transform1(potential,potential_dft,mtor,nrad) !calculating 1d DFT of s(:,:) along the first dimension
!!$     call toroidal_filter(potential_dft,mtor,nrad)
!!$     call oned_backward_fourier_transform1(potential_dft,potential,mtor,nrad)
!!$  endif

!!$  if(filter_radial.eqv..true.) then !filter over the radial mode number, keeping only low-radial-harmonics of the perturbation
!!$     call oned_sine_transform2(potential,potential_dst,mtor,nrad) !calculating 1d DST of s(:,:) along the second dimension
!!$     call radial_sine_filter_core(potential_dst,mtor,nrad)
!!$     call oned_inverse_sine_transform2(potential_dst,potential,mtor,nrad)
!!$  endif

  call x_derivative  (potential(:,:,1),potential_psi)
  call y_derivative(potential(:,:,1),potential_alpha)
  call z_derivative(potential,potential_theta) !derivative along the magnetic field line
  !write(*,*) potential_theta(10,40),tmp1(10,40),'myid=',myid
  !call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  do i=1,mtor
     do j=1,nrad
        jeq=j
        ef_cyl_r(i,j)=-potential_psi(i,j)*grad_psi_r(ipol_eq,jeq)&
             & -potential_theta(i,j)*grad_theta_r(ipol_eq,jeq)-potential_alpha(i,j)*grad_alpha_r(ipol_eq,jeq)

        ef_cyl_z(i,j)=-potential_psi(i,j)*grad_psi_z(ipol_eq,jeq) &
             & -potential_theta(i,j)*grad_theta_z(ipol_eq,jeq)-potential_alpha(i,j)*grad_alpha_z(ipol_eq,jeq)
        !ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mc(ipol_eq,jeq)*grad_alpha_phi(ipol_eq,jeq) !wrong, an additional 1/R factor is wrongly included, a bug found on 2018-Jan.-4 evening
        ef_cyl_phi(i,j)=-potential_alpha(i,j)/r_mc(ipol_eq,jeq) !corrected
     enddo
  enddo


  normal=tu*kev/(ln*bn*vn_fk*qu) !;write(*,*) 'normal=',normal
  !write(*,*) 'normal=',normal
  ef_cyl_r=ef_cyl_r*normal !normalized by the unit used in the code (Bn*vn_fk)
  ef_cyl_z=ef_cyl_z*normal
  ef_cyl_phi=ef_cyl_phi*normal

  do j=1,nrad
     do i=1,mtor
        ef_cyl_r_left(i,j)=ef_cyl_r(i,j) !store the data in the proper arrays
        ef_cyl_z_left(i,j)=ef_cyl_z(i,j)
        ef_cyl_phi_left(i,j)=ef_cyl_phi(i,j)
     enddo
     ef_cyl_r_left(mtor+1,j)=ef_cyl_r_left(1,j)  !peroidic toroidal boundary condition
     ef_cyl_z_left(mtor+1,j)=ef_cyl_z_left(1,j) 
     ef_cyl_phi_left(mtor+1,j)=ef_cyl_phi_left(1,j)
  enddo

  ef_cyl_r_left(:,1)=0._p_
  ef_cyl_z_left(:,1)=0._p_
  ef_cyl_phi_left(:,1)=0._p_

  ef_cyl_r_left(:,nrad)=  0._p_
  ef_cyl_z_left(:,nrad)=  0._p_
  ef_cyl_phi_left(:,nrad)= 0._p_

  call communicate_field_value_between_neighbour_cells2() !use ER,Ephi, EZ
end subroutine potential_to_cylindrical_field



module ampere
  use constants, only: p_
  implicit none
  complex(p_), allocatable :: mpara(:,:,:) !coefficient matrix for parallel Ampere's law
  integer, allocatable :: ipiv(:,:)

contains
  subroutine prepare_ampere_matrix()
    use constants, only: ii, zero, one, two, pi, twopi, kev, epsilon0, c
    use gk_module, only: mass_gk, charge_gk, nsm, gk_flr
    use magnetic_coordinates, only: nrad, nrad, toroidal_range, &
         & xgrid, xlow, xupp, grad_psi, grad_alpha,grad_psi_dot_grad_alpha
    use domain_decomposition, only: ipol_eq
    use gk_profile_funcs, only : gkn_func
    use control_parameters, only : nh_max
    implicit none
    real(p_) :: gy0, gx0, gxdgy0
    integer :: j, jp, m, n, nx, jeq, ns, info
    real(p_) :: part1, lx, ly, x
    real(p_) :: skin(nrad-2), wp2(nrad, nsm)
    complex(p_) :: s, part2

    nx = nrad - 1
    allocate(mpara(nx-1,nx-1,0:nh_max))
    allocate(ipiv(nx-1,0:nh_max))

    do ns = 1, nsm
       do j = 1, nrad
          x = xgrid(j)
          wp2(j,ns) = gkn_func(x,ns)*charge_gk(ns)**2/(mass_gk(ns)*epsilon0) !plasma frequency**2 
       enddo
    enddo

    do j = 1, nrad-2
       skin(j) = -sum(wp2(1+j,:))/c**2  !sum over species, skin current coefficient
    enddo

    lx = xupp - xlow
    ly = toroidal_range
    do n = 0, nh_max !toroidal harmonics
       do j = 1, nx-1
          jeq = j + 1
          gx0 = grad_psi(ipol_eq, jeq)
          gy0 = grad_alpha(ipol_eq, jeq)
          gxdgy0 = grad_psi_dot_grad_alpha(ipol_eq, jeq)
          do jp = 1, nx-1
             s = 0._p_
             do m = 1, nx-1
                part1 = (m*pi/lx)**2*gx0**2 + (n*twopi/ly)**2*gy0**2
                part2 = -ii*n*twopi/ly*(m*pi/lx)*2*gxdgy0
                s = s + sin(jp*m*pi/nx) * (part1*sin(j*m*pi/nx) + part2*cos(j*m*pi/nx))
             enddo
             mpara(j,jp,n) = s*two/nx  
          enddo
       enddo
    enddo

    do j = 1, nx-1
       mpara(j,j,:) = mpara(j,j,:) - skin(j)
    enddo

    do n = 0, nh_max !for each toroidal Fourier component
       call ZGETRF(nx-1, nx-1, mpara(:,:,n), nx-1, ipiv(:,n), info) !LU factorize the radial coeff matrix
    enddo
  end subroutine prepare_ampere_matrix


  subroutine solve_ampere(isolve, jpar_left, jpar_right, apara_s, apara_h, apara, ax, ay, az, ahx, ahy, ahz)
    use constants, only: p_, c, kev, epsilon0
    use normalizing, only: vu, tu, qu, nu
    use magnetic_coordinates, only: m=>mtor, n=>nrad, dtor,dradcor, jacobian
    use gk_module, only: w_unit
    use control_parameters, only: ismooth, nh_min, nh_max
    use derivatives_in_xyz, only: x_derivative0, y_derivative, z_derivative, z_derivative0
    use communication_connection, only: merge_source, update_scalar_at_right_boundary, &
         & update_derivatives_at_right_boundary
    use transform_module
    use domain_decomposition, only: myid, numprocs, dvol
    use filter_module
    use smoothing_module
    use math, only: solver_toroidal_mode_number_parallel
    implicit none
    integer, intent(in) :: isolve
    real(p_), intent(in) :: jpar_left(:,:), jpar_right(:,:)
    real(p_), intent(in) :: apara_s(:,:,:)
    real(p_), intent(out), dimension(:,:,:) :: apara_h
    real(p_), intent(out), dimension(:,:,:) :: apara, ax, ay, az, ahx, ahy, ahz
    real(p_) :: jpar(m,n) !, lap_As(m,n)
    real(p_) :: rhs(m, n-2) !, residual(m, n-2)
    complex(p_) :: rhs_dft(0:m-1, n-2), ah_dft(0:m-1, n-2)
    complex(p_) :: lap_As_dft(m, n), t_dft(0:m-1, n-2), ax_dft(0:m-1, n-2),  ahx_dft(0:m-1, n-2)
    real(p_) :: debye_sq
    integer :: kn, i, j, jeq

    debye_sq = (epsilon0*tu*kev)/(nu*qu**2)

    call merge_source(jpar_left, jpar_right, jpar)
    do j = 1, n  
       jpar(:,j)  = jpar(:,j)*(w_unit/dvol(j)/nu)
    enddo

    !call laplacian(apara_s(:,:,1), lap_As)
    !rhs(:,:) = lap_As(:, 2:n-1) + jpar(:, 2:n-1)/debye_sq*(vu/c)**2
    rhs(:,:) = jpar(:, 2:n-1)/debye_sq*(vu/c)**2

    call oned_DFT_parallel_version(rhs, rhs_dft, m, n-2) !calculating DFT along the first dimension
    call laplacian_dft(apara_s(:, :, 1), lap_As_dft(:, :), m, n)
    rhs_dft = rhs_dft + lap_As_dft(:, 2:n-1)

    ah_dft = 0._p_ !only some harmonics are solved, others are assumed to be zero
    call solver_toroidal_mode_number_parallel(mpara, IPIV, rhs_dft, ah_dft)
    do kn = 1, nh_max ! negative toroidal mode number
       ah_dft(m-kn,:) = conjg(ah_dft(kn,:)) 
    enddo

    do kn = nh_min, nh_max
       !call mfilter_for_each_n(ah_dft(kn,:), n-2, kn)
    enddo
    do kn = 1, nh_max ! negative toroidal mode number
       ah_dft(m-kn,:) = conjg(ah_dft(kn,:)) 
    enddo

    !if(nh_min == 0) call mfilter_for_n0(ah_dft(0,:), n-2)
    if(nh_min == 0) call surface_average_of_n0(ah_dft(0, :), n-2)

    call oned_backward_DFT_parallel_version(ah_dft, apara_h(:, 2:n-1, 1), m, n-2)
    apara_h(:, 1, 1) = 0 ! boundary condition
    apara_h(:, n, 1) = 0 ! boundary condition
    call update_scalar_at_right_boundary(apara_h)

    if(isolve==1) then
       !call x_derivative  (apara_h(:,:,1), ahx(:,:,1))
       call x_derivative0(ah_dft, ahx_dft) !more accurate than taking x derivative in real space
       call oned_backward_DFT_parallel_version(ahx_dft, ahx(:, 2:n-1, 1), m, n-2)

       call y_derivative(apara_h(:,:,1), ahy(:,:,1))
       call z_derivative  (apara_h(:,:,:), ahz(:,:,1))
       call update_derivatives_at_right_boundary(ahx, ahy, ahz)
    endif

    apara = apara_h + apara_s

    call oned_DFT_parallel_version(apara(:, 2:n-1, 1), t_dft, m, n-2) !DFT along the first dimension
    do kn = nh_min, nh_max
       !call mfilter_for_each_n(t_dft(kn,:), n-2, kn)
    enddo
    do kn = 1, nh_max ! negative toroidal mode number
       t_dft(m-kn,:) = conjg(t_dft(kn,:)) 
    enddo

    !if(nh_min == 0) call mfilter_for_n0(t_dft(0,:), n-2)
    if(nh_min == 0) call surface_average_of_n0(t_dft(0, :), n-2)

    call oned_backward_DFT_parallel_version(t_dft, apara(:, 2:n-1, 1), m, n-2)
    call update_scalar_at_right_boundary(apara)

    !call x_derivative(apara(:,:,1), ax(:,:,1))
    call x_derivative0(t_dft, ax_dft) !more accurate than taking x derivative in real space
    call oned_backward_DFT_parallel_version(ax_dft, ax(:, 2:n-1, 1), m, n-2)

    call y_derivative(apara(:,:,1), ay(:,:,1))
    call z_derivative(apara(:,:,:), az(:,:,1))
    call update_derivatives_at_right_boundary(ax, ay, az)
  end subroutine solve_ampere

  
  subroutine apara_s_evolution(aOld, aNew, phiz, dtao)
    use constants, only: p_
    use magnetic_coordinates, only: m=>mtor, n=>nrad
    use control_parameters, only: nh_min, nh_max
    use domain_decomposition, only: ipol_eq
    use table_in_mc, only: w2
    use transform_module
    use filter_module
    implicit none
    real(p_), intent(in) :: aOld(:,:), dtao
    real(p_), intent(in) :: phiz(:,:,:)
    real(p_), intent(out) :: aNew(:,:)
    complex(p_) :: a_dft(0:m-1, n)
    real(p_) :: rate
    integer :: i, j, kn
    ! use dapara_s/dt =-b\cdot\grad_phi
    do i = 1, m
       do j = 2, n-1
          rate = -phiz(i,j,1) * w2(ipol_eq, j)
          aNew(i,j) = aOld(i,j) + rate * dtao
       enddo
    enddo
    aNew(:,1)=0; aNew(:,n)=0 !zero radial boundary condition

    !use dapara_s/dt =0, worse than the above, tested
    !aNew = aOld

  end subroutine apara_s_evolution

  subroutine apara_resplit_and_weight_pullback(w_gk, apara_s, apara_h, ahx, ahy, ahz)
    use constants, only: p_, kev
    use gk_module, only:  nsm, vn_gk, mass_gk, v_gk, charge_gk, xgc, vpar_gk, nm_gk, &
         & zgc, x_ring, y_ring, z_ring, gk_flr, w_unit, ps_vol_gk
    use magnetic_coordinates, only: xlow, xupp
    use normalizing, only: qu, tu, vu
    !use load_gk_mod, only: maxwellian
    use load_gk_mod, only: f0
    use gk_profile_funcs, only: gkt_func
    use gyro_average_mod, only: gyro_average0
    implicit none
    real(p_), intent(inout) :: w_gk(:,:), apara_s(:,:,:), apara_h(:,:,:)
    real(p_), intent(out) :: ahx(:,:,:), ahy(:,:,:), ahz(:,:,:)
    real(p_) :: x, te, ah_av, v, E, factor
    integer :: k, ns
    
    do ns = 1, nsm 
       do k = 1, nm_gk(ns)
          x = xgc(k,ns)
          if ( (x > xupp) .or. (x < xlow)) cycle
          te = gkt_func(x,ns)
          v = v_gk(k,ns)
          !E = 0.5*mass_gk(ns)*(v*vn_gk(ns))**2/kev
          call gyro_average0(gk_flr(ns), x_ring(:,k,ns), y_ring(:,k,ns), z_ring(:,k,ns), &
               &   apara_h, ah_av)
          factor = (charge_gk(ns)/qu)/(te/tu)*(vn_gk(ns)/vu)
          !weight pullback
          !w_gk(k,ns) = w_gk(k,ns) -factor*vpar_gk(k,ns)*ah_av*maxwellian(x,E,ns)* ps_vol_gk(k,ns)*vn_gk(ns)**3/w_unit
          w_gk(k,ns) = w_gk(k,ns) -factor*vpar_gk(k,ns)*ah_av*f0(x,v,ns)*ps_vol_gk(k,ns)/w_unit
       enddo
    enddo
    !Apar resplit
    apara_s = apara_s + apara_h !collect apara_h into apara_s
    apara_h = 0
    ahx = 0
    ahy = 0
    ahz = 0
  end subroutine apara_resplit_and_weight_pullback


  subroutine laplacian(apara, out)
    use constants, only: p_
    use magnetic_coordinates, only: m=>mtor, n=>nrad,  &
         & grad_psi, grad_psi_dot_grad_alpha, grad_alpha
    use domain_decomposition, only: ipol_eq
    use derivatives_in_xyz, only: x_derivative, y_derivative
    implicit none
    real(p_), intent(in) :: apara(m,n)
    real(p_), intent(out) :: out(m,n)
    real(p_), dimension(m,n) :: apara_x, apara_xx, apara_y, apara_yy, apara_xy
    integer :: j, jeq

    call x_derivative(apara, apara_x)
    call x_derivative(apara_x, apara_xx)
    call y_derivative(apara, apara_y)
    call y_derivative(apara_y, apara_yy)
    call y_derivative(apara_x, apara_xy)

    do j = 1, n
       jeq = j
       out(:,j) = grad_psi  (ipol_eq, jeq)**2*apara_xx(:,j)  &
            &   + grad_alpha(ipol_eq, jeq)**2*apara_yy(:,j)  &
            &   + 2*grad_psi_dot_grad_alpha(ipol_eq, jeq)*apara_xy(:,j)

    enddo
  end subroutine laplacian


  subroutine laplacian_dft(As, out, m, n)
    use constants, only: p_, ii, twopi
    use magnetic_coordinates, only:  ly=>toroidal_range, &
         & grad_psi, grad_psi_dot_grad_alpha, grad_alpha
    use domain_decomposition, only: ipol_eq
    use control_parameters, only: nh_min, nh_max
    use derivatives_in_xyz, only: x_derivative0, y_derivative
    use transform_module
    implicit none
    integer, intent(in) :: m, n
    real(p_), intent(in) :: As(0:m-1,n)
    complex(p_), intent(out) :: out(0:m-1,n)
    complex(p_) :: As_dft(0:m-1, n), As_dft_x(0:m-1,n), As_dft_xx(0:m-1,n)
    integer :: j, jeq, nh, kn

    call oned_DFT_parallel_version(As(:, 2:n-1), As_dft(:, 2:n-1), m, n-2) !calculating DFT along the first dimension
    As_dft(:,1)=0; As_dft(:,n)=0
    call x_derivative0(As_dft, As_dft_x)
    call x_derivative0(As_dft_x, As_dft_xx)

    do nh = 0, m-1
       if (nh .le. m/2) then
          kn = nh
       else
          kn = nh - m
       endif
       do j = 1, n
          jeq = j
          out(nh,j) = grad_psi  (ipol_eq, jeq)**2*As_dft_xx(nh,j)  &
               &   - grad_alpha(ipol_eq, jeq)**2*As_dft(nh,j)*(kn*twopi/ly)**2  &
               &   + 2*grad_psi_dot_grad_alpha(ipol_eq, jeq)*As_dft_x(nh,j)*(ii*kn*twopi/ly)
       enddo
    enddo
  end subroutine laplacian_dft
  
  subroutine skin_current_residual(iter, isolve, apara_h, lap_As, jpar_ref, residual)
    use constants,only:zero,one,two,pi,twopi,kev,epsilon0, c, mu0
    use normalizing, only: qu, vu, nu, tu
    use magnetic_coordinates,only: m=>mtor, n=>nrad, xgrid, dtor, dradcor, &
         & ygrid, jacobian, xgrid
    use gk_profile_funcs, only : gkn_func, gkt_func
    use gk_module, only: mass_gk, charge_gk, nsm, vn_gk,  nm_gk, w_unit, ps_vol_gk, &
         & xgc, zgc, ygc, vpar_gk, xgc_mid, zgc_mid, ygc_mid, vpar_gk_mid, v_gk
    use load_gk_mod, only : f0
    use communication_connection, only: merge_source
    use domain_decomposition,only: dtheta2,theta_start, ipol_eq, multi_eq_cells, myid, tclr, dvol
    use gyro_average_mod, only : field_at_particle0
    implicit none
    integer, intent(in) :: iter, isolve
    real(p_), intent(in) :: apara_h(:,:,:), lap_As(:,:), jpar_ref(:,:)
    real(p_), intent(out) :: residual(m,n-2)
    real(p_)  :: debye_sq, wp2, ah_e
    real(p_) :: jpar_left(m,n), jpar_right(m,n), jpar(m,n)
    real(p_) :: x, y, z, te, delta_f, w, vpar
    integer :: ns, i, j, k, i_plus1, j_plus1, jeq
    real(p_) :: cz1, cz2, cy1, cy2, cx1, cx2, kernel

    debye_sq = (epsilon0*tu*kev)/(nu*qu**2)
    do ns=1,nsm
       if(charge_gk(ns)<0) exit !select out electrons
    enddo

    jpar_left = 0; jpar_right = 0;
    do k = 1, nm_gk(ns)     ! Monte-Carlo integration
       if(isolve==1) then
          x = xgc_mid(k,ns)
          y = ygc_mid(k,ns)
          z = zgc_mid(k,ns)
          vpar = vpar_gk_mid(k,ns)
       else
          x = xgc(k,ns)
          y = ygc(k,ns)
          z = zgc(k,ns)
          vpar = vpar_gk(k,ns)
       endif
       te = gkt_func(x,ns)
       call field_at_particle0(x, y, z, apara_h, ah_e)
       delta_f = -charge_gk(ns)/(te*kev)*vpar*vn_gk(ns)*ah_e*(Tu*kev/(qu*vu))*f0(x,v_gk(k,ns),ns)
       w = ps_vol_gk(k,ns)*delta_f / w_unit
       kernel = w*(charge_gk(ns)/qu)*vpar*(vn_gk(ns)/vu)

       cz1 = (z-theta_start)/dtheta2 
       cz2 = one-cz1

       i = floor((y-ygrid(1))/dtor+1) 
       cy1 = (y-ygrid(i))/dtor
       cy2 = one-cy1

       j = floor((x-xgrid(1))/dradcor+1)
       cx1 = (x-xgrid(j))/dradcor
       cx2 = one-cx1

       i_plus1=i+1
       if(i.eq.m) i_plus1=1 !periodic condition
       j_plus1=j+1
       if(j.eq.n) cycle !marker is out of radial computational region

       jpar_left(i,j) = jpar_left(i,j) + kernel*cz2*cy2*cx2
       jpar_left(i_plus1,j) = jpar_left(i_plus1,j) + kernel*cz2*cy1*cx2
       jpar_left(i,j_plus1) = jpar_left(i,j_plus1) + kernel*cz2*cy2*cx1
       jpar_left(i_plus1,j_plus1) = jpar_left(i_plus1,j_plus1)+kernel*cz2*cy1*cx1

       jpar_right(i,j) = jpar_right(i,j) + kernel*cz1*cy2*cx2
       jpar_right(i_plus1,j) = jpar_right(i_plus1,j)+ kernel*cz1*cy1*cx2
       jpar_right(i,j_plus1) = jpar_right(i,j_plus1) + kernel*cz1*cy2*cx1
       jpar_right(i_plus1,j_plus1) = jpar_right(i_plus1,j_plus1)+kernel*cz1*cy1*cx1
    enddo

    call merge_source(jpar_left, jpar_right, jpar)
    do j=1,n  !divided by cell volume (gridpoint is the center of the cell)
       jpar(:,j)  = jpar(:,j)*(w_unit/dvol(j)/nu)
    enddo

    do j =1, n-2
       x = xgrid(j+1)
       wp2= gkn_func(x,ns)*charge_gk(ns)**2/(mass_gk(ns)*epsilon0)
       residual(:,j) = jpar(:,j+1)/debye_sq*(vu/c)**2 -(-wp2/c**2*apara_h(:,j+1, 1))
    enddo
    ! if(tclr==0) then
    !    if(isolve==1 .and. iter==2 ) write(*,'(A40, 2i4, 10ES14.4)') 'ipol_eq, iter, jpar1, jpar2=', iter, ipol_eq, &
    !         & jpar(m/2,n/2)/debye_sq*(vu/c)**2, -wp2/c**2*apara_h(m/2,n/2, 1), lap_As(m/2,n/2), jpar_ref(m/2,n/2)/debye_sq*(vu/c)**2
    ! endif
  end subroutine skin_current_residual
  
end module ampere
module force
contains

!!$subroutine field_perturbation_on_marker(radcor,theta,alpha,active,&
!!$     & epar0,ex0,ey0,mf_par0,mf_x0,mf_y0)



!!$subroutine potential_on_marker(radcor,theta,alpha,active,potential_val)
!!$  use constants,only:p_
!!$  use constants,only: one,zero,twopi
!!$  use perturbation_field,only: potential_left,potential_right
!!$  use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
!!$  use domain_decomposition,only: dtheta2,theta_start
!!$  use interpolate_module,only: linear_2d_interpolate
!!$  implicit none
!!$  real(p_),intent(in)::radcor,theta,alpha
!!$  logical,intent(in)::active
!!$  real(p_),intent(out)::potential_val
!!$  real(p_):: coeff1,coeff2,tmp1,tmp2
!!$
!!$ if(active.eqv..false.) then !force on particles outside the computational region is set to zero
!!$   potential_val=0._p_
!!$  else
!!$
!!$  coeff1=(theta-theta_start)/dtheta2
!!$  coeff2=one-coeff1
!!$  !if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
!!$  call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,potential_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
!!$  call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,potential_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
!!$  potential_val=tmp1*coeff2+tmp2*coeff1
!!$endif
!!$end subroutine potential_on_marker


!!$  subroutine epar_on_marker(radcor,theta,alpha,epar_val)
!!$    use constants,only:p_
!!$    use constants,only: one,zero,twopi
!!$    use perturbation_field,only: epar_left,epar_right
!!$    !use perturbation_field,only: epar_left=>epar_left_half,epar_right=>epar_right_half
!!$    use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
!!$    use domain_decomposition,only: dtheta2,theta_start
!!$    use interpolate_module,only: linear_2d_interpolate
!!$    implicit none
!!$    real(p_),intent(in)::radcor,theta,alpha
!!$    real(p_),intent(out)::epar_val
!!$    real(p_):: coeff1,coeff2,tmp1,tmp2
!!$
!!$    coeff1=(theta-theta_start)/dtheta2
!!$    coeff2=one-coeff1
!!$
!!$    !if(alpha>twopi .or. alpha<zero) write(*,*) 'field_on_marker,alpha=',alpha
!!$    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,epar_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
!!$    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,epar_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
!!$    epar_val=tmp1*coeff2+tmp2*coeff1
!!$
!!$  end subroutine epar_on_marker


  subroutine field_perturbation_on_marker_for_adiabatic_e_model(radcor,theta,alpha,&
       & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val) !for adiabatic electron model
    use constants,only:p_
    use constants,only: one,zero,twopi
    use perturbation_field,only:ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
    use perturbation_field,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right !field at nearby (larger) poloidal location
    use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
    use domain_decomposition,only: dtheta2,theta_start !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_),intent(in)::radcor,theta,alpha
    real(p_),intent(out)::ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
    real(p_):: coeff1,coeff2,tmp1,tmp2

    coeff1=(theta-theta_start)/dtheta2
    coeff2=one-coeff1

    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
    ef_cyl_r_val=tmp1*coeff2+tmp2*coeff1

    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_left,alpha,radcor,tmp1)  
    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_right,alpha,radcor,tmp2)  
    ef_cyl_z_val=tmp1*coeff2+tmp2*coeff1

    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_left,alpha,radcor,tmp1)
    call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_right,alpha,radcor,tmp2)  
    ef_cyl_phi_val=tmp1*coeff2+tmp2*coeff1
  end subroutine field_perturbation_on_marker_for_adiabatic_e_model

  subroutine field_perturbation_on_marker2(radcor,theta,alpha,active,&
       & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val) !using cylindrical components
    use constants,only:p_
    use constants,only: one,zero,twopi
    use perturbation_field,only:ef_cyl_r_left,ef_cyl_z_left,ef_cyl_phi_left
    use perturbation_field,only: ef_cyl_r_right,ef_cyl_z_right,ef_cyl_phi_right !field at nearby (larger) poloidal location
    use magnetic_coordinates,only: mtor,nrad,ygrid,xgrid
    use domain_decomposition,only: dtheta2,theta_start !as input
    use interpolate_module,only: linear_2d_interpolate
    implicit none
    real(p_),intent(in)::radcor,theta,alpha
    logical,intent(in):: active
    real(p_),intent(out)::ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
    real(p_):: coeff1,coeff2,tmp1,tmp2

    if(active.eqv..false.) then !force on particles outside the computational region is set to zero
       ef_cyl_r_val=0._p_
       ef_cyl_z_val=0._p_
       ef_cyl_phi_val=0._p_
    else

       coeff1=(theta-theta_start)/dtheta2
       coeff2=one-coeff1

       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_left,alpha,radcor,tmp1)  !uniform xarray and zarray are assumed
       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_r_right,alpha,radcor,tmp2)  !uniform xarray and zarray are assumed
       ef_cyl_r_val=tmp1*coeff2+tmp2*coeff1

       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_left,alpha,radcor,tmp1)  
       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_z_right,alpha,radcor,tmp2)  
       ef_cyl_z_val=tmp1*coeff2+tmp2*coeff1

       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_left,alpha,radcor,tmp1)
       call linear_2d_interpolate(mtor+1,nrad,ygrid,xgrid,ef_cyl_phi_right,alpha,radcor,tmp2)  
       ef_cyl_phi_val=tmp1*coeff2+tmp2*coeff1
    endif
  end subroutine field_perturbation_on_marker2
end module force



subroutine push_electron_cylindrical_single_particle(dtao,mu,r,z,phi,vpar)
  !input: initial condition of the orbit: mu,r,z,phi,vpar
  !Output: the instanteous value of the orbit after dtao: mu,r,z,phi,vpar
  use constants,only:p_, zero,one,two,one_half,three,six,twopi
  implicit none
  real(p_),intent(in):: dtao,mu
  real(p_),intent(inout):: r,z,phi,vpar  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvpar
  real(p_):: r_e_dot,z_e_dot,phi_e_dot,vpar_e_dot !equations of motion
  real(p_):: kr1,kz1,kvpar1,kphi1,kr2,kz2,kvpar2,kphi2,kr3,kz3,kvpar3,kphi3,kr4,kz4,kvpar4,kphi4 !Runge-Kutta steps

  !4nd order Rung-Kutta method
  kr1=dtao*r_e_dot(r,z,phi,vpar,mu)
  kz1=dtao*z_e_dot(r,z,phi,vpar,mu)
  kvpar1=dtao*vpar_e_dot(r,z,phi,vpar,mu)
  kphi1=dtao*phi_e_dot(r,z,phi,vpar,mu)

  kr2=dtao*r_e_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kz2=dtao*z_e_dot      (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kvpar2=dtao*vpar_e_dot(r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)
  kphi2=dtao*phi_e_dot  (r+kr1*one_half,z+kz1*one_half,phi+kphi1*one_half,vpar+kvpar1*one_half,mu)

  kr3=dtao*r_e_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kz3=dtao*z_e_dot      (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kvpar3=dtao*vpar_e_dot(r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)
  kphi3=dtao*phi_e_dot  (r+kr2*one_half,z+kz2*one_half,phi+kphi2*one_half,vpar+kvpar2*one_half,mu)

  kr4=dtao*r_e_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kz4=dtao*z_e_dot      (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kvpar4=dtao*vpar_e_dot(r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)
  kphi4=dtao*phi_e_dot  (r+kr3,z+kz3,phi+kphi3,vpar+kvpar3,mu)

  dr=kr1/six+kr2/three+kr3/three+kr4/six
  dz=kz1/six+kz2/three+kz3/three+kz4/six
  dphi=kphi1/six+kphi2/three+kphi3/three+kphi4/six
  dvpar=kvpar1/six+kvpar2/three+kvpar3/three+kvpar4/six
  !update
  r=r+      dr
  z=z+      dz
  vpar=vpar+dvpar
  phi=phi+  dphi
  !write(*,*) 'dvpar=',dvpar
end subroutine push_electron_cylindrical_single_particle



function r_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_, twopi
    use gk_module,only: charge_sign_gk  
use magnetic_field, only : b,bphi,b_z,bz,b_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_r,bstar_e_parallel,bstar_e_parallelval
  
  

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=bstar_e_r(r,z,phi,vpar)/bstar_e_parallelval*vpar +&
       & charge_sign_gk(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*&
       & (bphi(r,z)*b_z(r,z)-bz(r,z)*b_phi(r,z)/r)

end function r_e_dot


function z_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_, twopi
    use gk_module,only: charge_sign_gk
use magnetic_field, only : b,bphi,b_z,bz,b_phi, br,b_r
    implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_z,bstar_e_parallel,bstar_e_parallelval
  
 

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval= bstar_e_z(r,z,phi,vpar)/bstar_e_parallelval*vpar + &
       & charge_sign_gk(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*(br(r,z)*b_phi(r,z)/r-bphi(r,z)*b_r(r,z))
end function z_e_dot

function phi_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_
  use constants,only:twopi
  use gk_module,only: charge_sign_gk
  use magnetic_field, only : b,bphi,b_z,bz,b_phi,br,b_r
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_phi,bstar_e_parallel,bstar_e_parallelval
   

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=bstar_e_phi(r,z,phi,vpar)/bstar_e_parallelval*vpar+&
       & charge_sign_gk(1)*mu/(twopi*b(r,z)*bstar_e_parallelval)*(bz(r,z)*b_r(r,z)-br(r,z)*b_z(r,z))
  funcval=funcval/r
end function phi_e_dot

function vpar_e_dot(r,z,phi,vpar,mu) result(funcval)
  use constants,only:p_
  use constants,only:twopi
use magnetic_field, only : b,bphi,b_z,bz,b_phi, b_r
  implicit none
  real(p_):: funcval,r,z,phi,vpar,mu
  real(p_):: bstar_e_r,bstar_e_z,bstar_e_phi,bstar_e_parallel,bstar_e_parallelval
  

  bstar_e_parallelval=bstar_e_parallel(r,z,phi,vpar)
  funcval=-mu*(bstar_e_r(r,z,phi,vpar)/bstar_e_parallelval*b_r(r,z)+bstar_e_z(r,z,phi,vpar)/bstar_e_parallelval*b_z(r,z) &
       & +bstar_e_phi(r,z,phi,vpar)/bstar_e_parallelval*b_phi(r,z)/r )
end function vpar_e_dot


function bstar_e_r(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_gk
use magnetic_field, only : br,unitbphi_z,unitbz_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  

  funcval=br(r,z)+charge_sign_gk(1)*vpar/twopi*(unitbz_phi(r,z)/r-unitbphi_z(r,z))

end function bstar_e_r


function bstar_e_z(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only:charge_sign_gk
use magnetic_field, only : b,bphi,bz,unitbphi_r,unitbr_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_):: unitbphi
  

  unitbphi=bphi(r,z)/b(r,z)
  funcval=bz(r,z)+charge_sign_gk(1)*vpar/twopi*(unitbphi_r(r,z)+unitbphi/r-unitbr_phi(r,z)/r)

end function bstar_e_z

function bstar_e_phi(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_gk
use magnetic_field, only : bphi,unitbr_z,unitbz_r
    implicit none
  real(p_)::funcval,r,z,phi,vpar
  

  funcval=bphi(r,z)+charge_sign_gk(1)*vpar/twopi*(unitbr_z(r,z)-unitbz_r(r,z))

end function bstar_e_phi


function bstar_e_parallel(r,z,phi,vpar) result(funcval)
  use constants,only:p_
  use constants,only:twopi,one
    use gk_module,only: charge_sign_gk
use magnetic_field, only : b,br,bz,bphi, unitbr_z,unitbz_r,unitbphi_r,unitbphi_z, unitbr_phi,unitbz_phi
    implicit none
  real(p_):: funcval,r,z,phi,vpar
  real(p_)::unitb_dot_curl_unitb
  real(p_):: unitbr,unitbphi,unitbz,bval
  
  bval=b(r,z)
  unitbr=br(r,z)/bval
  unitbz=bz(r,z)/bval
  unitbphi=bphi(r,z)/bval

  unitb_dot_curl_unitb=unitbr*(unitbz_phi(r,z)/r-unitbphi_z(r,z))&
       & +unitbphi*(unitbr_z(r,z)-unitbz_r(r,z))&
       & +unitbz*(unitbphi_r(r,z)+unitbphi/r-unitbr_phi(r,z)/r)
  funcval=bval+charge_sign_gk(1)*vpar/twopi*unitb_dot_curl_unitb
end function bstar_e_parallel

module boris
contains

subroutine push_full_orbit_cylindrical_boris(charge_sign,dtao,r,phi,z,vr,vphi,vz)
  !input: initial condition of the orbit: r,phi,z,vr,vphi,vz
  !Output: the instanteous value of the orbit after dtao: r,phi,z,vr,vphi,vz
!in essence, the algorithm uses Cartesian coordinates, and then takes into account the rotation of the the basis vectors due to the change of particle's location
  use constants,only:p_
  use constants,only:zero,one,two,twopi
  use math,only: cross_product_in_cartesian
  use magnetic_field, only : br,bz,bphi 
  implicit none
  real(p_),intent(in):: charge_sign, dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of orbit

!  real(p_):: er,ez,ephi !function names
  real(p_):: tx,ty,tz,t,factor,sx,sy,sz
  real(p_):: vx_minus,vy_minus,vz_minus
  real(p_):: vx_prime,vy_prime,vz_prime
  real(p_):: vx_plus,vy_plus,vz_plus
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: alpha

  vx_minus=vr  +  er  (r,z,phi)*dtao/two*(twopi*charge_sign)
  vy_minus=vphi+  ephi(r,z,phi)*dtao/two*(twopi*charge_sign)
  vz_minus=vz  +  ez  (r,z,phi)*dtao/two*(twopi*charge_sign)

  tx=br(r,z)*dtao/two*(twopi*charge_sign)
  ty=bphi(r,z)*dtao/two*(twopi*charge_sign)
  tz=bz(r,z)*dtao/two*(twopi*charge_sign)

  call cross_product_in_cartesian(vx_minus,vy_minus,vz_minus,tx,ty,tz,cx,cy,cz)

  vx_prime=  vx_minus+ cx
  vy_prime=  vy_minus+ cy
  vz_prime=  vz_minus+ cz

  t=sqrt(tx*tx+tz*tz+ty*ty)
  factor=two/(one+t*t)
  sx=tx*factor
  sy=ty*factor
  sz=tz*factor

  call cross_product_in_cartesian(vx_prime,vy_prime,vz_prime,sx,sy,sz,cx,cy,cz)

  vx_plus =  vx_minus+ cx
  vy_plus =  vy_minus+ cy
  vz_plus =  vz_minus+ cz

  vx=vx_plus+  er  (r,z,phi)*dtao/two*(twopi*charge_sign)
  vy=vy_plus+  ephi(r,z,phi)*dtao/two*(twopi*charge_sign)
  vz=vz_plus+  ez  (r,z,phi)*dtao/two*(twopi*charge_sign)


  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  alpha=asin(y/r)

  phi=phi+alpha

  vr=cos(alpha)*vx+sin(alpha)*vy
  vphi=-sin(alpha)*vx+cos(alpha)*vy
end subroutine push_full_orbit_cylindrical_boris


subroutine push_full_orbit_cylindrical_boris2(dtao,radcor,theta,alpha, active,&
     & r,phi,z,vr,vphi,vz,vr_mid,vphi_mid,vz_mid)
  !derived from "push_full_orbit_cylindrical_boris", the modification is that: additional input phi_mid, which is value of phi at t_{n+1/2},
  !and additional output: vr_mid,vphi_mid,vz_mid, which are projection of velocity at t_{n+1/2} onto the basis vector at t_{n+1/2}, instead of t_{n+1}. Here, for general cases, n can also be an half-integer
  use constants,only:p_, zero,one,two,twopi
  use fk_module,only: fk_nonlinear
  use magnetic_field, only : br,bz,bphi 
  use math,only: cross_product_in_cartesian
  use force
  implicit none
  real(p_),intent(in):: dtao
  real(p_),intent(inout):: r,phi,z,vr,vphi,vz  !instantaneous value of (r,phi,z)  at t_{n} (as input) and t_{n+1} (as output), (vr,vphi,vz) is the projection of velocity at t_{n-1/2} at basis vector at t_{n} (as input); and the projection of velocity at t_{n+1/2} on the basis vector at t_{n+1} (as output}
  !  real(p_),intent(in):: phi_mid !value of phi at t_{n+1/2}, not used now, instead this value is calculated direclty in this subroutine
  real(p_),intent(in):: radcor,theta,alpha !perturbed field is interpolated in magnetic coordinates, so we need the magnetic coordinates of particles
  logical,intent(in):: active
  real(p_),intent(out):: vr_mid,vphi_mid,vz_mid !projection velocity at t_{n+1/2} on the basis vector at t_{n+1/2}
  !real(p_):: er,ez,ephi !function names
  real(p_):: er_val,ez_val,ephi_val
  real(p_):: tx,ty,tz,t,factor,sx,sy,sz
  real(p_):: vx_minus,vy_minus,vz_minus
  real(p_):: vx_prime,vy_prime,vz_prime
  real(p_):: vx_plus,vy_plus,vz_plus
  real(p_):: cx,cy,cz
  real(p_):: x,y,vx,vy
  real(p_):: dphi,phi0


  if(fk_nonlinear.eq.1) then
     call  field_perturbation_on_marker2(radcor,theta,alpha, active, er_val,ez_val,ephi_val) 
  else
     er_val=0._p_
     ez_val=0._p_
     ephi_val=0._p_
  endif
  phi0=phi !record the old value of phi
  vx_minus=vr   +  er_val*dtao/two*(twopi)
  vy_minus=vphi +  ephi_val*dtao/two*(twopi)
  vz_minus=vz   +  ez_val*dtao/two*(twopi)

  tx=br(r,z)*dtao/two*(twopi)
  ty=bphi(r,z)*dtao/two*(twopi)
  tz=bz(r,z)*dtao/two*(twopi)

  call cross_product_in_cartesian(vx_minus,vy_minus,vz_minus,tx,ty,tz,cx,cy,cz)

  vx_prime=  vx_minus+ cx
  vy_prime=  vy_minus+ cy
  vz_prime=  vz_minus+ cz

  t=sqrt(tx*tx+tz*tz+ty*ty)
  factor=two/(one+t*t)
  sx=tx*factor
  sy=ty*factor
  sz=tz*factor

  call cross_product_in_cartesian(vx_prime,vy_prime,vz_prime,sx,sy,sz,cx,cy,cz)

  vx_plus =  vx_minus+ cx
  vy_plus =  vy_minus+ cy
  vz_plus =  vz_minus+ cz

  vx=vx_plus +  er_val*dtao/two*(twopi)
  vy=vy_plus +  ephi_val*dtao/two*(twopi)
  vz=vz_plus +  ez_val*dtao/two*(twopi)

  x=r+vx*dtao
  y=vy*dtao
  z=z+vz*dtao


  r=sqrt(x*x+y*y)

  dphi=asin(y/r)

  phi=phi+dphi

  vr=cos(dphi)*vx+sin(dphi)*vy !projected to basis vectors at the t_{n+1} particle location
  vphi=-sin(dphi)*vx+cos(dphi)*vy !projected to basis vectors at the t_{n+1} particle location

  !  dphi=phi_mid-phi0 !turn out not reliable because the toroidal angular momentum conservation is not as good as the following simple treatment
  dphi=(phi-phi0)/two !equivalent to dphi=dphi/two
  vr_mid=cos(dphi)*vx+sin(dphi)*vy !the projection of v_{n+1/2} on basis vectors at t_{n+1/2}
  vphi_mid=-sin(dphi)*vx+cos(dphi)*vy
  vz_mid=vz

end subroutine push_full_orbit_cylindrical_boris2


function er(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: er,r,z,phi

  er=0._p_

end function er


function ez(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: ez,r,z,phi

  ez=0._p_

end function ez


function ephi(r,z,phi) !R component of electric field
  use constants,only:p_
  implicit none
  real(p_):: ephi,r,z,phi

  ephi=0._p_

end function ephi





subroutine normalize_full_orbit_variables(nmarker,r,phi,z, vr,vphi,vz)
 use constants,only:p_
 use normalizing,only: Ln
 use fk_module, only: vn_fk

  implicit none
  integer,intent(in):: nmarker
  real(p_),intent(inout):: r(nmarker),phi(nmarker),z(nmarker)
  real(p_),intent(inout):: vr(nmarker),vphi(nmarker),vz(nmarker)

 !normalization
  r=r/Ln !convert to unit Ln
  z=z/Ln !convert to unit Ln
!  phi=phi
  vr=vr/vn_fk !normalized by vn_fk
  vphi=vphi/vn_fk
  vz=vz/vn_fk
end subroutine normalize_full_orbit_variables


subroutine particle_variables_to_guiding_center_variables(mass, charge, vn, r,phi,z,vr0,vphi0,vz0,rg,phig,zg, mu, vpar)
  use constants,only:p_,two
  use normalizing,only: bn
  use magnetic_field, only : br,bz,bphi 
  implicit none
  real(p_),intent(in):: mass, charge, vn
  real(p_),intent(in):: r,phi,z,vr0,vphi0,vz0 !particle variables, vr0, vphi0, vz0 in unit of vn
  real(p_),intent(out):: rg,phig,zg, mu, vpar !guiding-center variables

  real(p_):: brval,bphival,bzval,bval,factor
  real(p_):: vr,vz,vphi,v !in SI unit
  
  brval=br(r,z)
  bzval=bz(r,z)
  bphival=bphi(r,z)
  bval=sqrt(brval**2+bphival**2+bzval**2)
  vr=vr0*vn !to SI unit
  vz=vz0*vn !to SI unit
  vphi=vphi0*vn !to SI unit

  v=sqrt(vr**2+vz**2+vphi**2)

factor=mass/(bval**2*charge)

  !  rg=r+mass/(bval**2*charge)*(-vz*bphival+vphi*bzval) ! wrong!
  rg=sqrt((r+factor*(vphi*bzval-vz*bphival))**2+(factor*(vz*brval-vr*bzval))**2)
  !  phig=phi+atan(mass/(bval**2*charge)*(-vr*bzval+vz*brval)/r) !wrong!
  phig=phi+asin(factor*(-vr*bzval+vz*brval)/rg) !range of asin function is [-pi/2,pi/2]
  !zg=z+mass/(bval**2*charge)*(-vr*bphival-vphi*brval) !wrong, pointed by Yingfeng Xu
  zg=z+factor*(vr*bphival-vphi*brval) !corrected.

vpar=(vr*brval+vz*bzval+vphi*bphival)/bval
mu=(v**2-vpar**2)/(two*bval)
vpar=vpar/vn
mu=mu/(vn**2/bn)
end subroutine particle_variables_to_guiding_center_variables


subroutine guiding_center_variables_to_particle_variables(mass, charge, r_gc,z_gc,phi_gc, mu, vpar, r,z,phi, vr,vz,vphi) !see my notes for the formula
    use constants,only:p_
    use constants,only: twopi,two,one,four
    use normalizing,only:ln,bn
    use math,only: cross_product_in_cartesian
    use domain_decomposition,only: myid
    use magnetic_field, only : br,bz,bphi
    implicit none
    real(p_),intent(in):: mass, charge
    real(p_),intent(in):: r_gc,z_gc,phi_gc, mu, vpar
    real(p_),intent(out):: r,z,phi, vr,vz,vphi
    real(p_):: vn
    real(p_)::br_val,bz_val,bphi_val,bval,omega, gyro_radius
    real(p_):: mu_si, vpar_si,vper_si  !*_si indicates quanties in S.I. units
    real(p_):: root1,root2,a,b,c !coefficients of the quadratic equation
    real(p_):: vx,vy, vz2
    real(p_):: shiftx, shifty, shiftz
    integer,parameter:: niteration=2 !one iteration is usually enought
    integer:: k
    vn=ln/(twopi/(bn*charge/mass))
    mu_si=mu*(mass*vn**2/bn)
    vpar_si=vpar*vn

    r=r_gc !initial guess
    z=z_gc
    phi=phi_gc 
    do k=1, niteration
       br_val=br(r,z)
       bz_val=bz(r,z)
       bphi_val=bphi(r,z)
       bval=sqrt(br_val**2+bz_val**2+bphi_val**2)

       omega=bval*abs(charge)/mass
       vper_si=sqrt(mu_si*two*bval/mass)
!!$    gyro_radius=vper_si/omega
!!$    r=r_gc+gyro_radius
!!$    z=z_gc
!!$    phi=phi_gc

       vr=0._p_
       a=(bz_val/bphi_val)**2+one
       b=-two*bval*bz_val/bphi_val**2*vpar_si
       c=(bval*vpar_si/bphi_val)**2-two*bval*mu_si/mass-vpar_si**2
       root1=(-b+sqrt(b**2-four*a*c))/(two*a)
       root2=(-b-sqrt(b**2-four*a*c))/(two*a)
!!$       if((bphi_val*charge).gt.0) then
!!$          vz=root1/vn
!!$       else
!!$          vz=root2/vn
!!$       endif


       !          vz=root1/vn
       vz=root2/vn

       vphi=(bval*vpar-vz*bz_val)/bphi_val

       vx=vr*vn
       vy=vphi*vn
       vz2=vz*vn

       call cross_product_in_cartesian(vx,vy,vz2,br_val,bphi_val,bz_val,shiftx,shifty,shiftz)

       shiftx=-shiftx/(bval**2*charge/mass)
       shifty=-shifty/(bval**2*charge/mass)
       shiftz=-shiftz/(bval**2*charge/mass)

       r=r_gc + shiftx
       z=z_gc + shiftz
       phi=phi_gc+asin(shifty/r)
       !if(myid.eq.0)write(*,*) 'iteration=', k, 'r, z, phi=', r, z, phi
    enddo

!!$    write(*,*) 'v_squared=',vphi**2+vz**2, vpar**2+mu_si*two*bval/mass/vn**2 !to very they agree with each other
!!$    yj: block
!!$      real(p_):: rg2,phig2,zg2, mu2, vpar2
!!$      call particle_variables_to_guiding_center_variables(mass, charge, vn, r,phi,z,vr,vphi,vz,rg2,phig2,zg2, mu2, vpar2)
!!$     write(*,*) 'rg,phig,zg, mu, vpar=', rg2,r_gc, phig2, phi_gc, zg2, z_gc, mu2, mu, vpar2, vpar !to very they agree with each other
!!$    end block yj
  end subroutine guiding_center_variables_to_particle_variables


end module boris



module initial_half_step_for_boris
implicit none
contains
subroutine backward_half_step_for_boris(charge_sign,dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm
  !input:  r,z,phi,vr,vz,vphi, at the same time t_{0}
  !Output: the projection of velocity at t_{-1/2} onto the basis vectors at t_{0}
  use constants,only:p_
  use constants,only:one_half
  use constants,only:zero,one,two,one_half,three,six,twopi,kev
  implicit none
  real(p_),intent(in):: charge_sign,dtao
  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
  real(p_):: dr,dz,dphi,dvr,dvz,dvphi
!  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot
!  real(p_):: vz_fo_dot !equations of motion
!  real(p_):: vr_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
  real(p_):: step,vr_new,vphi_new,r0,z0,phi0,dt
  integer,parameter:: m=100
  integer:: i

!write(*,*)  "energy_before evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
r0=r
phi0=phi
z0=z

step=-0.5_p_*dtao
dt=step/m
do i=1,m
  kr1=    one_half*dt*r_fo_dot(r,z,phi,vr,vz,vphi)
  kz1=    one_half*dt*z_fo_dot(r,z,phi,vr,vz,vphi)
  kphi1=  one_half*dt*phi_fo_dot(r,z,phi,vr,vz,vphi)
  kvr1=   one_half*dt*vr_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  kvphi1= one_half*dt*vphi_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  kvz1=   one_half*dt*vz_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)

  dr=   dt*r_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dz=   dt*z_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dphi= dt*phi_fo_dot  (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvr=  dt*vr_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvphi=dt*vphi_fo_dot (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
  dvz=  dt*vz_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)

  z=z+dz
  r=r+dr
  phi=phi+dphi
  vz=vz+dvz
  vr=vr+dvr
  vphi=vphi+dvphi
enddo

dphi=phi-phi0
vr_new=vr
vphi_new=vphi

!write(*,*)  "energy_before projection=", 0.5_p_*mass*(vr_new**2+vz**2+vphi_new**2)*vn**2/kev
!projection of the new velocity onto the old basis vectors
vr=  -vphi_new*sin(dphi)+vr_new*cos(dphi) 
vphi=vphi_new*cos(dphi)+vr_new*sin(dphi) 
!write(*,*)  "energy_after projection=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
r=r0 !go back to the original value
z=z0 !go back to the original value
phi=phi0 !go back to the original value
!write(*,*) 'dvphi=',dvphi
end subroutine backward_half_step_for_boris


subroutine forward_half_step_for_boris(charge_sign,dtao,r0,z0,phi0,vr0,vz0,vphi0,r1,z1,phi1,vr1,vz1,vphi1) 
  !input: initial condition of the orbit: r0,z0,phi0,vr0,vz0,vphi0, at the same time t_{0}
  !Output: the instanteous value of (r,z,phi) after half_dtao, i.e., t_{1/2}, and the projection of velocity at t_{0} onto the basis vectors at t_{1/2}
  use constants,only:p_
  use constants,only:zero,one,two,one_half,twopi
  use math,only: shift_toroidal
  use magnetic_coordinates, only : toroidal_range
  implicit none
  real(p_),intent(in):: charge_sign,dtao
  real(p_),intent(in):: r0,z0,phi0,vr0,vz0,vphi0  !instantaneous value of orbit at t0
  real(p_),intent(out):: r1,z1,phi1 !location at t0+half_dtao
  real(p_),intent(out):: vr1,vz1,vphi1   !projection of the old velocity (t=t0) on the new basis vector (determined by the new (r,z,phi))
  real(p_):: step,dt,dr,dz,dphi,dvr,dvz,dvphi
!  real(p_):: r_fo_dot,z_fo_dot,phi_fo_dot,vr_fo_dot,vz_fo_dot,vphi_fo_dot !equations of motion
  real(p_):: kr1,kz1,kphi1,kvr1,kvz1,kvphi1 !Runge-Kutta steps
!  real(p_):: kr2,kz2,kphi2,kvr2,kvz2,kvphi2 !Runge-Kutta steps
  integer,parameter::m=100 ! if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple steps
  integer:: k
  real(p_):: r,z,phi,vr,vz,vphi !working variables

  r=r0
  z=z0
  phi=phi0
  vr=vr0
  vz=vz0
  vphi=vphi0

  step=0.5_p_*dtao
  dt=step/m
  do k=1,m
     !2nd order Rung-Kuta method
     kr1=     one_half*dt*r_fo_dot(r,z,phi,vr,vz,vphi)
     kz1=     one_half*dt*z_fo_dot(r,z,phi,vr,vz,vphi)
     kphi1=   one_half*dt*phi_fo_dot(r,z,phi,vr,vz,vphi)
     kvr1=    one_half*dt*vr_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
     kvz1=    one_half*dt*vz_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
     kvphi1=  one_half*dt*vphi_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)

     dr=    dt*r_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dz=    dt*z_fo_dot    (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dphi=  dt*phi_fo_dot  (r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvr=   dt*vr_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvz=   dt*vz_fo_dot   (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)
     dvphi= dt*vphi_fo_dot (charge_sign,r+kr1,z+kz1,phi+kphi1,vr+kvr1,vz+kvz1,vphi+kvphi1)

     !update
     r=r+      dr
     z=z+      dz
     phi=phi+  dphi
     vr=vr+dvr
     vz=vz+dvz
     vphi=vphi+dvphi
     !write(*,*) 'dvphi=',dvphi
  enddo

r1=r
z1=z
phi1=phi

dphi=phi1-phi0

vz1=vz0
vr1=  vphi0*sin(dphi)+vr0*cos(dphi) !projection of the old velocity on the new basis vector 
vphi1=vphi0*cos(dphi)-vr0*sin(dphi) !projection of the old velocity on the new basis vector 

!!$phi1=phi1-int(phi1/twopi)*twopi !shift into the range [0:twopi]
!!$if(phi1.lt.0) phi1=phi1+twopi !shift into the range [0:twopi]
!   call shift_to_zero_twopi_range(phi1)
   call shift_toroidal(phi1, toroidal_range)
end subroutine forward_half_step_for_boris



function r_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: r_fo_dot,r,z,phi,vr,vz,vphi

  r_fo_dot=vr
end function r_fo_dot

function phi_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: phi_fo_dot,r,z,phi,vr,vz,vphi

  phi_fo_dot=vphi/r
end function phi_fo_dot

function z_fo_dot(r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  implicit none
  real(p_):: z_fo_dot,r,z,phi,vr,vz,vphi

  z_fo_dot=vz
end function z_fo_dot

function vr_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  use constants,only:p_
  use constants,only:twopi
  use magnetic_field, only : bphi,bz
  implicit none
  real(p_):: charge_sign,vr_fo_dot,r,z,phi,vr,vz,vphi
 
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi)) !wrong, term vphi**2/r**2 is missing
!  vr_fo_dot=twopi*(-vz*bphi(r,z,phi)+vphi*bz(r,z,phi))+ vphi**2/r**2 !2017-Aug.12, a bug found, where r**2 should be replaced by r. kinetic energy is not well conserved by the buggy code (compared with results given by the Cartesian version of orbit integrator), which forced me to examine the codes to find possible bugs, and I finally found this bug. After correcting this code, the conservation of kinetic energy given by this code is as good as that of the Cartesian version.
  vr_fo_dot=twopi*charge_sign*(-vz*bphi(r,z)+vphi*bz(r,z))+ vphi**2/r !correct
end function vr_fo_dot



function vphi_fo_dot(charge_sign,r,z,phi,vr,vz,vphi) 
  use constants,only:p_
  use constants,only:twopi
  use magnetic_field, only : br,bz
  implicit none
  real(p_):: charge_sign,vphi_fo_dot,r,z,phi,vr,vz,vphi

  vphi_fo_dot=twopi*charge_sign*(-vr*bz(r,z)+vz*br(r,z))-vphi*vr/r !the last term is inertial force

end function vphi_fo_dot

function vz_fo_dot(charge_sign,r,z,phi,vr,vz,vphi)
  use constants,only:p_,twopi
  use magnetic_field, only : br,bphi
  implicit none
  real(p_):: charge_sign,vz_fo_dot,r,z,phi,vr,vz,vphi
  vz_fo_dot=twopi*charge_sign*(vr*bphi(r,z)-vphi*br(r,z))
end function vz_fo_dot
end module initial_half_step_for_boris



subroutine test_full_orbit(dtao0)
  use constants,only:p_
  use constants,only:zero,kev,pi,twopi
  use fk_module,only: vn=>vn_fk,tn=>tn_fk
  use fk_module,only:mass=>mass_i,charge=>charge_i
  use boris
  use initial_half_step_for_boris
  use magnetic_field, only : psi_func, br, bz,bphi,b
  implicit none
  real(p_),intent(in):: dtao0
  real(p_),parameter::ln=1.0_p_
  real(p_)::r,z,phi,vr,vz,vphi,t !,r0,z0,phi0
  integer:: kk
  integer,parameter:: maxstep=29000
  real(p_):: kin_eng,pphi
  real(p_):: bval,brval,bzval,bphival
  real(p_):: v,energy,sai,pitch_angle
  real(p_):: rg,zg,phig,rg1,zg1,phig1
  real(p_):: b_dot_v,omega_local,dtao
  integer,parameter:: n_tor_period=1
  logical,parameter:: check_boundary_loss=.true.
  character(100), parameter::orbit_file="fo_go.txt"
  real(p_):: r0,z0,phi0,vr0,vphi0,vz0
  r= 2.1_p_
  z= 0._p_
  phi=0._p_
  vr=1.0d6
  vz=1.0d6
  vphi=5d5
!  dtao=1._p_

dtao=1.0*dtao0

  !--guiding-center orbit, to roughly verify the full orbit--
  bval=b(r,z)
  brval=br(r,z)
  bzval=bz(r,z)
  bphival=bphi(r,z)

!!$  call particle_to_guiding_center_location(ns,r,phi,z,vr,vphi,vz,brval,bphival,bzval,rg,phig,zg)
!!$
!!$  v=sqrt(vr**2+vz**2+vphi**2)
!!$  energy=0.5_p_*mass*v*v/kev
!!$  write(*,*) 'kinetic energy (kev)=', energy
!!$
!!$  b_dot_v=brval*vr+bphival*vphi+bzval*vz
!!$  sai=b_dot_v/(bval*v)
!!$  pitch_angle=acos(sai)/pi*180._p_
!!$
!!$  call orbit(mass,charge,energy,pitch_angle,phig,rg,zg,dtao,n_tor_period,check_boundary_loss,orbit_file)

  !----
  !then calculate full orbit

  r=r/ln
  z=z/ln
  vr=vr/vn
  vz=vz/vn
  vphi=vphi/vn
!r=1.6182173678310701;z=-0.17532304528354634;phi=0.51458392767352690
!vr= 2.2153551525963856E-002; vz= -3.2974526832309299E-003; vphi=-3.1659760504767286E-002

!   omega_local=b(r*ln,z*ln)*charge/mass
!   dtao=twopi/omega_local/tn/8_p_ !the time-step is chosen as in terms of the local gyro-period
t=0._p_
!write(*,*) 'local gyro-period=',twopi/omega_local
 
     write(*,*) 'Using Boris algorithm to push full orbit'
     open(163,file='full_orbit_boris.txt')
     !call backward_half_step_cartesian(dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm
     r0=r;z0=z;phi0=phi; vr0=vr; vz0=vz; vphi0=vphi
     call forward_half_step_for_boris(sign(1._p_,charge),dtao,r0,z0,phi0,vr0,vz0,vphi0,r,z,phi,vr,vz,vphi) !for testing

!     do kk=1,maxstep
     do kk=1,200
        call push_full_orbit_cylindrical_boris(1.d0,dtao,r,phi,z,vr,vphi,vz)

        t=t+dtao
        kin_eng=0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
        pphi=mass*r*ln*vphi*vn+charge*psi_func(r*ln,z*ln)
        !write(163,*) t, t*tn,r*ln,z*ln,phi,vr,vz,vphi,kin_eng, (vr**2+vz**2)/b(r,z),pphi   
        write(163,*) t+dtao/2, (t+dtao/2)*tn,r*ln,z*ln,phi,vr,vz,vphi,kin_eng, (vr**2+vz**2)/b(r,z),pphi !for the case of using forward initialization

     enddo
     write(*,*) 'kk=',kk
     close(163)

end subroutine test_full_orbit



!!$subroutine backward_half_step_cartesian(dtao,r,z,phi,vr,vz,vphi) !push only velocity, to set initial condition for the first step of boris algorithm, !using multi-steps, instead of one step, considering dtao used in Boris may be comparable to the gyro-period
!!$!actually working in Cartesian coordinates (i.e., constant basis vectors)
!!$  use constants,only:p_
!!$  use constants,only:zero,one,two,one_half,three,six,twopi,kev
!!$  use normalizing,only: vn=>vn_fk
!!$  use fk_module,only:mass=>mass_i
!!$  implicit none
!!$  real(p_),intent(in):: dtao
!!$  real(p_),intent(inout):: r,z,phi,vr,vz,vphi  !instantaneous value of orbit
!!$  real(p_):: dvx,dvy,dvz
!!$  real(p_):: x_fo_dot,y_fo_dot,z_cartesian_fo_dot
!!$  real(p_):: vx_fo_dot,vy_fo_dot,vz_cartesian_fo_dot
!!$  real(p_):: kx1,ky1,kz1,kvx1,kvy1,kvz1 !Runge-Kutta steps
!!$!  real(p_):: kx2,ky2,kz2,kvx2,kvy2,kvz2 !Runge-Kutta steps
!!$  real(p_)::vx,vy,x,y,dx,dy,dz,z0,dt
!!$  real(p_):: step
!!$  integer,parameter::m=100 !if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple rk steps
!!$  integer::k
!!$
!!$write(*,*)  "energy calculated in Cartesian ,before evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
!!$  x=r
!!$  y=0._p_
!!$  z0=z
!!$  vx=vr
!!$  vy=vphi
!!$
!!$step=-0.5_p_*dtao
!!$dt=step/m
!!$do k=1,m
!!$  kx1=one_half*dt*x_fo_dot(x,y,z,vx,vy,vz)
!!$  ky1=one_half*dt*y_fo_dot(x,y,z,vx,vy,vz)
!!$  kz1=one_half*dt*z_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  kvx1=   one_half*dt*vx_fo_dot(x,y,z,vx,vy,vz)
!!$  kvy1=   one_half*dt*vy_fo_dot(x,y,z,vx,vy,vz)
!!$  kvz1=   one_half*dt*vz_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$
!!$  dx=dt*x_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dy=dt*y_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dz=dt*z_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvx=   dt*vx_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvy=   dt*vy_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvz=   dt*vz_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$
!!$  !update
!!$  x=x+dx
!!$  y=y+dy
!!$  z=z+dz
!!$  vx=vx+dvx
!!$  vy=vy+dvy
!!$  vz=vz+dvz
!!$enddo
!!$
!!$z=z0 !resume to the original value
!!$vr=vx
!!$vphi=vy
!!$
!!$write(*,*)  "energy calculated in Cartesian after evolution=", 0.5_p_*mass*(vr**2+vz**2+vphi**2)*vn**2/kev
!!$!write(*,*) 'dvphi=',dvphi
!!$end subroutine backward_half_step_cartesian
!!$
!!$
!!$function x_fo_dot(x,y,z,vx,vy,vz)
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: x_fo_dot,x,y,z,vx,vy,vz
!!$
!!$  x_fo_dot=vx
!!$end function 
!!$
!!$
!!$function y_fo_dot(x,y,z,vx,vy,vz)
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: y_fo_dot,x,y,z,vx,vy,vz
!!$
!!$  y_fo_dot=vy
!!$end function 
!!$
!!$
!!$function z_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: z_cartesian_fo_dot,x,y,z,vx,vy,vz
!!$
!!$  z_cartesian_fo_dot=vz
!!$end function
!!$
!!$ 
!!$function vx_fo_dot(x,y,z,vx,vy,vz) 
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vx_fo_dot,x,y,z,vx,vy,vz
!!$  real(p_):: bphi,bz,br !function names
!!$  real(p_):: brval,bphival,bzval,by
!!$  real(p_)::r,phi
!!$  r=sqrt(x*x+y*y)
!!$  phi=acos(x/r)
!!$  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
!!$  brval=br(r,z,phi)
!!$  bphival=bphi(r,z,phi)
!!$  bzval=bz(r,z,phi)
!!$  by=brval*sin(phi)+bphival*cos(phi)
!!$  vx_fo_dot=twopi*(-vz*by+vy*bzval) 
!!$end function
!!$
!!$
!!$function vy_fo_dot(x,y,z,vx,vy,vz) !without inertial force
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vy_fo_dot,x,y,z,vx,vy,vz
!!$  real(p_):: br,bz,bphi  !function names
!!$  real(p_):: brval,bphival,bzval,bx
!!$  real(p_)::r,phi
!!$  r=sqrt(x*x+y*y)
!!$  phi=acos(x/r)
!!$  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
!!$  brval=br(r,z,phi)
!!$  bphival=bphi(r,z,phi)
!!$  bzval=bz(r,z,phi)
!!$  
!!$  bx=brval*cos(phi)-bphival*sin(phi)
!!$
!!$  vy_fo_dot=twopi*(-vx*bzval+vz*bx)
!!$
!!$end function
!!$
!!$
!!$function vz_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  use constants,only:p_
!!$  use constants,only:twopi
!!$  implicit none
!!$  real(p_):: vz_cartesian_fo_dot,x,y,z,vx,vy,vz
!!$  real(p_):: bphi,br ,bz !function names
!!$  real(p_):: brval,bphival,bzval,bx,by
!!$  real(p_)::r,phi
!!$  r=sqrt(x*x+y*y)
!!$  phi=acos(x/r)
!!$  if(y<0) phi=-phi  !phi is assumed in the range [-pi,pi]
!!$  brval=br(r,z,phi)
!!$  bphival=bphi(r,z,phi)
!!$  bzval=bz(r,z,phi)
!!$
!!$  bx=brval*cos(phi)-bphival*sin(phi)
!!$  by=brval*sin(phi)+bphival*cos(phi)
!!$ vz_cartesian_fo_dot=twopi*(vx*by-vy*bx)
!!$
!!$end function
!!$
!!$
!!$
!!$
!!$
!!$subroutine forward_half_step_cartesian(dtao,r0,z0,phi0,vr0,vz0,vphi0,r1,z1,phi1,vr1,vz1,vphi1) 
!!$  use constants,only:p_
!!$  use constants,only:zero,one,two,one_half,three,six,twopi,kev
!!$  use normalizing,only: vn=>vn_fk
!!$  use fk_module,only:mass=>mass_i
!!$  implicit none
!!$  real(p_),intent(in):: dtao
!!$  real(p_),intent(in):: r0,z0,phi0,vr0,vz0,vphi0
!!$  real(p_),intent(out):: r1,z1,phi1,vr1,vz1,vphi1
!!$  real(p_):: x_fo_dot,y_fo_dot,z_cartesian_fo_dot
!!$  real(p_):: vx_fo_dot,vy_fo_dot,vz_cartesian_fo_dot
!!$  real(p_):: kx1,ky1,kz1,kvx1,kvy1,kvz1 !Runge-Kutta steps
!!$!  real(p_):: kx2,ky2,kz2,kvx2,kvy2,kvz2 !Runge-Kutta steps
!!$  real(p_):: x,y,z,vx,vy,vz,dx,dy,dz,dvx,dvy,dvz !working variables
!!$  real(p_):: step,dt,alpha
!!$  integer,parameter::m=100 !if dtao is comparable or larger than the gyro-period, then the first backward half step needs to be finished with multiple rk steps
!!$  integer::k
!!$
!!$write(*,*)  "energy calculated in Cartesian ,before evolution=", 0.5_p_*mass*(vr0**2+vz0**2+vphi0**2)*vn**2/kev
!!$  x=r0
!!$  y=0._p_
!!$  z=z0
!!$  vx=vr0
!!$  vy=vphi0
!!$  vz=vz0
!!$
!!$step=0.5_p_*dtao
!!$dt=step/m
!!$do k=1,m
!!$  kx1=one_half*dt*x_fo_dot(x,y,z,vx,vy,vz)
!!$  ky1=one_half*dt*y_fo_dot(x,y,z,vx,vy,vz)
!!$  kz1=one_half*dt*z_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$  kvx1=   one_half*dt*vx_fo_dot(x,y,z,vx,vy,vz)
!!$  kvy1=   one_half*dt*vy_fo_dot(x,y,z,vx,vy,vz)
!!$  kvz1=   one_half*dt*vz_cartesian_fo_dot(x,y,z,vx,vy,vz)
!!$
!!$  dx=dt*x_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dy=dt*y_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dz=dt*z_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvx=   dt*vx_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvy=   dt*vy_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$  dvz=   dt*vz_cartesian_fo_dot(x+kx1,y+ky1,z+kz1,vx+kvx1,vy+kvy1,vz+kvz1)
!!$
!!$  !update
!!$  x=x+dx
!!$  y=y+dy
!!$  z=z+dz
!!$  vx=vx+dvx
!!$  vy=vy+dvy
!!$  vz=vz+dvz
!!$enddo
!!$
!!$ r1=sqrt(x*x+y*y)
!!$
!!$  alpha=asin(y/r1)
!!$
!!$  phi1=phi0+alpha
!!$z1=z
!!$
!!$
!!$  vr1=cos(alpha)*vr0+sin(alpha)*vphi0
!!$  vphi1=-sin(alpha)*vr0+cos(alpha)*vphi0
!!$vz1=vz0
!!$write(*,*)  "energy calculated in Cartesian after evolution=", 0.5_p_*mass*(vr1**2+vz1**2+vphi1**2)*vn**2/kev
!!$!write(*,*) 'dvphi=',dvphi
!!$end subroutine forward_half_step_cartesian
module fk_particle_coordinates_transform_module
contains
subroutine compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !get the radial coordinates of markers, relocate markers if they are outside of the computational box, then calculate the (theta,alpha) coordinates
  use constants,only:p_
  use constants,only:pi,twopi
  use poloidal_flux_2d,only:xarray,zarray,nx,nz
  use magnetic_coordinates,only:xlow,xupp
  use magnetic_coordinates,only:mpol,nrad,zgrid, xgrid, tor_shift_mc, toroidal_range
  use magnetic_field,only:radcor_as_func_of_pfn
  use magnetic_field, only : pfn_func
    use map_to_mc, only : interpolate_from_cylindrical_to_magnetic_coordinates1

  use interpolate_module
    use math,only: shift_toroidal
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(inout):: r_i(nmarker_i),phi_i(nmarker_i),z_i(nmarker_i)
  real(p_),intent(out):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  logical,intent(out):: active_i(nmarker_i),touch_bdry_i(nmarker_i)
  integer:: k,flag,next_seed
  real(p_):: theta,rannum1,rannum2
  logical:: outside_box
  real(p_):: tor_shift

  do k=1,nmarker_i
     outside_box=r_i(k).ge.xarray(nx) .or.r_i(k).le.xarray(1) .or. z_i(k).ge.zarray(nz) .or.  z_i(k).le.zarray(1) 
     if(outside_box.eqv..true.) then
        touch_bdry_i(k)=.true. !marker is lost forever, will never be re-introduced
        active_i(k)=.false.
     else
        radcor_i(k)=radcor_as_func_of_pfn(pfn_func(r_i(k),z_i(k))) !calculate the radial coordinate
        if(radcor_i(k).ge.xgrid(nrad) .or. radcor_i(k).le.xgrid(1)) then 
           touch_bdry_i(k)=.true. !marker is lost forever, will never be re-introduced
           active_i(k)=.false.
        else if(radcor_i(k).lt.xlow .or. radcor_i(k).gt.xupp) then
           z_i(k)=-z_i(k) !relocate the marker by up-down reversing
           radcor_i(k)=radcor_as_func_of_pfn(pfn_func(r_i(k),z_i(k))) !re-calculate the radial coordinate, which is different from the original value if the flux-surface is not up-down symmetric
           active_i(k)=.false.
           touch_bdry_i(k)=.false.
        else
           touch_bdry_i(k)=.false.
           active_i(k)=.true.
        endif
     endif
  enddo

!calculate (theta,alpha) coordinates
 do k=1,nmarker_i
     !if(active_i(k).eqv..false.) cycle !wrong, inactive markers' theta must be computed so that they can be sorted by the sorting subroutine
     if(touch_bdry_i(k).eqv..true.) cycle
     call interpolate_from_cylindrical_to_magnetic_coordinates1(r_i(k),z_i(k),theta_i(k))
     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,tor_shift_mc,&
          & theta_i(k),radcor_i(k),tor_shift) !interpolating in magnetic coordinates to get tor_shift
     alpha_i(k)=phi_i(k)-tor_shift !generalized toroidal angle
     call shift_toroidal(alpha_i(k),toroidal_range)
     call shift_toroidal(phi_i(k),toroidal_range) !phi_i is needed in deposition, so should be shifted to the desired range
  enddo

end subroutine compute_particle_magnetic_coordinates


subroutine count_lost_markers_fk() 
  use fk_module
  use domain_decomposition,only: myid
  implicit none
  integer:: k, nlost
  nlost=0
  do k=1,nmarker_i
     if(touch_bdry_i(k).eqv..true.) then
        nlost=nlost+1
     endif
  enddo
write(*,*) 'number of lost fk markers=', nlost, 'myid=',myid, ',total number=', nmarker_i
end subroutine count_lost_markers_fk


subroutine clean_up_lost_markers_fk() !for fully kinetic markers, drop lost markers so that we we have smaller particle arrays
  use fk_module,only:nmarker_i,w_i, ps_vol_i, active_i, touch_bdry_i, radcor_i, theta_i, alpha_i,&
       & r_i,phi_i,z_i, vr_i,vphi_i,vz_i, vr_i_mid, vz_i_mid, vphi_i_mid
  implicit none
  integer:: k, kclean

  kclean=0
  do k=1,nmarker_i
     if(touch_bdry_i(k).eqv..false.) then

        kclean=kclean+1

        r_i(kclean)=r_i(k)
        phi_i(kclean)=phi_i(k)
        z_i(kclean) =z_i (k)
        vr_i(kclean)=vr_i(k)
        vphi_i(kclean)=vphi_i(k)
        vz_i(kclean)=vz_i(k)
        w_i(kclean)=w_i(k)
        ps_vol_i(kclean)=ps_vol_i(k)
        active_i(kclean)=active_i(k)
        touch_bdry_i(kclean)=touch_bdry_i(k)
        radcor_i(kclean)=radcor_i(k)
        theta_i(kclean)=theta_i(k)
        alpha_i(kclean)=alpha_i(k)
        vr_i_mid(kclean)=vr_i_mid(k)
        vz_i_mid(kclean)= vz_i_mid(k)
        vphi_i_mid(kclean)= vphi_i_mid(k)
     endif
  enddo
  nmarker_i=kclean !update the number of markers
end subroutine clean_up_lost_markers_fk


end module fk_particle_coordinates_transform_module
module drift
contains
  pure subroutine compute_drift(ns, lost, x, z, mu, vpar, &
       & phix, phiy, phiz, ax, ay, az, &
       & xdrift0, zdrift0, ydrift0, mirror_force, zdrift00, xdrift1, zdrift1, ydrift1)
    !output are used by the tajectory pusher and the weight pusher
    use constants, only: p_, one, twopi, kev
    use table_in_mc, only: w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, bdgxcgz, bdgycgz, bdgxcgy, b_mc
    use magnetic_coordinates, only: mpol, nrad, zgrid, xgrid, xlow, xupp
    use interpolate_module, only: linear_2d_interpolate
    use gk_module, only: nm_gk, charge_sign_gk, vn_gk
    use normalizing, only: vu, bn, ln, qu,tu
    implicit none
    integer,  intent(in)  :: ns ! species index
    logical, intent(in) :: lost(:)
    real(p_), intent(in)  :: x(:), z(:), mu(:), vpar(:)
    real(p_), intent(in)  :: phix(:), phiy(:), phiz(:), ax(:), ay(:), az(:)
    real(p_), intent(out) :: xdrift0(:), ydrift0(:), zdrift0(:), zdrift00(:)
    real(p_), intent(out) :: mirror_force(:)
    real(p_), intent(out) :: xdrift1(:), ydrift1(:), zdrift1(:)
    real(p_) :: w1v, w2v, w3v, w4v, w5v, w6v, w7v, w8v, w9v, w10v
    real(p_) :: factor1, cs, norm1, norm2, bval, bsq, bdgxcgy0, bdgxcgz0, bdgycgz0
    integer :: k

    cs = charge_sign_gk(ns)
    norm1 = tu*kev/(qu*vn_gk(ns)*bn*ln)
    norm2 = tu*kev/(qu*vu*bn*ln)
    do k = 1, nm_gk(ns) !number of MC markers
       if (lost(k) .eqv. .true.) cycle
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w1,z(k),x(k),w1v) 
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w2,z(k),x(k),w2v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w3,z(k),x(k),w3v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w4,z(k),x(k),w4v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w5,z(k),x(k),w5v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w6,z(k),x(k),w6v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w7,z(k),x(k),w7v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w8,z(k),x(k),w8v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w9,z(k),x(k),w9v)
       !call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,w10,z(k),x(k),w10v)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bdgxcgz,z(k),x(k),bdgxcgz0)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bdgycgz,z(k),x(k),bdgycgz0)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,bdgxcgy,z(k),x(k),bdgxcgy0)
       call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,b_mc, z(k),x(k), bval)
       bsq = bval**2
       factor1 = one + cs*vpar(k)/twopi*w1v
       !factor1=one !for testing, no difference
       xdrift0(k) = cs*vpar(k)**2/twopi*w3v/factor1 +cs*mu(k)/(twopi*factor1)*w6v 
       ydrift0(k) = cs*vpar(k)**2/(twopi*factor1)*w5v +cs*mu(k)/(twopi*factor1)*w8v
       zdrift0(k) = vpar(k)*(w2v+cs*vpar(k)/twopi*w4v)/factor1 &
            &  + cs*mu(k)/(twopi*factor1)*w7v
       !zdrift00(k) = zdrift0(k) - vpar(k)*w2v !removing the parallel streaming
       zdrift00(k) = vpar(k)*(cs*vpar(k)/twopi*w4v)/factor1 &
            &  + cs*mu(k)/(twopi*factor1)*w7v

       !mirror_force(k)=-mu(k)/factor1*(w9v+cs*vpar(k)/twopi*w10v)
       mirror_force(k) = - mu(k) * w9v !tested, agree with the result using the previous line
       !mirror_force2(k) = mu(k)/twopi*vpar(k)*w10v

       xdrift1(k) =  (-phiy(k)*bdgxcgy0 - phiz(k)*bdgxcgz0)/bsq*norm1 &
            &  + vpar(k)/bsq*(ay(k)*bdgxcgy0 + az(k)*bdgxcgz0)*norm2

       ydrift1(k)  = (phix(k)*bdgxcgy0 - phiz(k)*bdgycgz0)/bsq*norm1 & 
            &  + vpar(k)/bsq*(ax(k)*(-bdgxcgy0) + az(k)*bdgycgz0)*norm2

       zdrift1(k)  = (phix(k)*bdgxcgz0 + phiy(k)*bdgycgz0)/bsq*norm1 &
            &  + vpar(k)/bsq*(ax(k)*(-bdgxcgz0) + ay(k)*(-bdgycgz0))*norm2
    enddo

  end subroutine compute_drift
end module drift
module gk_trajectory_pusher
contains
  pure subroutine push_gc(ns, lost, dt, xdrift0, zdrift0, ydrift0, mirror_force,&
       & xdrift1, zdrift1, ydrift1, x_old, z_old, y_old, vpar_old, weight, x_new, z_new, y_new, vpar_new)
    use constants, only: p_, pi, one_half
    use gk_module, only: nm_gk, gk_nonlinear
    use math, only: shift_toroidal
    use magnetic_coordinates, only : toroidal_range, xlow, xupp
    use magnetic_field, only: qfunc
    implicit none
    integer, intent(in) :: ns
    logical, intent(in) :: lost(:)
    real(p_), intent(in) :: dt
    real(p_), intent(in) :: xdrift0(:), zdrift0(:), ydrift0(:), mirror_force(:)
    real(p_), intent(in) :: xdrift1(:), zdrift1(:), ydrift1(:)
    real(p_), intent(in) :: x_old(:), y_old(:), z_old(:), vpar_old(:)
    real(p_), intent(inout) :: weight(:)
    real(p_), intent(out) :: x_new(:), y_new(:), z_new(:), vpar_new(:)
    real(p_) :: xdrift, zdrift, ydrift
    integer :: nm, k

    nm = nm_gk(ns)
    do k = 1, nm
       if(lost(k) .eqv. .true.) then
          z_new(k) = z_old(k)
          x_new(k) = x_old(k)
          cycle
       endif

       if(gk_nonlinear(ns)==1) then
          xdrift = xdrift0(k) + xdrift1(k)
          zdrift = zdrift0(k) + zdrift1(k)
          ydrift = ydrift0(k) + ydrift1(k)
       else
          xdrift = xdrift0(k)
          zdrift = zdrift0(k)
          ydrift = ydrift0(k)
       endif

       x_new(k) = x_old(k) + xdrift*dt
       z_new(k) = z_old(k) + zdrift*dt
       y_new(k) = y_old(k) + ydrift*dt
       vpar_new(k)  = vpar_old(k) + mirror_force(k)*dt

       ! Comment out the following if-structure if you do not want particle refilling:
       if ((x_new(k) > xupp) .or. (x_new(k) < xlow) ) then
          x_new(k) = x_old(k)
          z_new(k) = - z_old(k)
          y_new(k) = y_old(k) + 2*qfunc(x_old(k))*z_old(k)
          weight(k) = 0
        endif

       if ((z_new(k) .ge. pi) .or. (z_new(k) < -pi)) &
            & call shift_gc_theta_then_alpha(x_new(k), z_new(k), y_new(k))
       call shift_toroidal(y_new(k), toroidal_range)
    enddo

  end subroutine push_gc

  pure  subroutine shift_gc_theta_then_alpha(radcor, theta, alpha)
    !shift theta, and then shift alpha so that phi is not changed, i.e., keep the particle in the same sptial location
    use constants, only: p_, twopi, pi
    use magnetic_field, only: qfunc
    implicit none
    real(p_), intent(in) :: radcor
    real(p_), intent(inout) :: theta, alpha
    integer :: ishift

    ishift = floor(theta/twopi)
    theta = theta - ishift * twopi
    if(theta .ge. pi) then
       theta = theta - twopi
       ishift = ishift + 1
    endif

    !alpha shift is needed to make the cylindrical toroidal angle phi be the same as the value before the theta-shifting,
    !i.e., make the particle lie at the same sptaial location when doing the theta-shifting
    alpha = alpha + ishift*twopi*qfunc(radcor) 
  end subroutine shift_gc_theta_then_alpha

end module gk_trajectory_pusher
module gk_weight_pusher
contains
  pure  subroutine push_gk_weight(ns, lost, dtao, radcor, vpar, v, &
       & xdrift0, ydrift0, zdrift0, zdrift00, mirror_force, xdrift1, ydrift1, zdrift1, &
       & phix, phiy, phiz, ahx, ahy, ahz, ah, w, w_new) 
    use constants, only: p_, two, kev
    use normalizing, only: vu, qu, tu
    use gk_module, only: charge_gk, mass_gk, vn_gk, ps_vol_gk, w_unit, gk_nonlinear, nm_gk
    use gk_profile_funcs, only: gkt_func, gkn_func, gkdndx_func, gkdtdx_func
    use load_gk_mod, only: f0 !equilibrium distribution function
    use gk_polarization, only: slowing_down_xd, slowing_down_Ed
    use magnetic_coordinates, only: xlow, xupp
    implicit none
    integer, intent(in) :: ns
    logical, intent(in) :: lost(:)
    real(p_), intent(in) :: dtao
    real(p_), intent(in) :: radcor(:), vpar(:), v(:)
    real(p_), intent(in) :: xdrift0(:), ydrift0(:), zdrift0(:), zdrift00(:)
    real(p_), intent(in) :: mirror_force(:)
    real(p_), intent(in) :: xdrift1(:), ydrift1(:), zdrift1(:)
    real(p_), intent(in) :: phix(:), phiy(:), phiz(:), ahx(:), ahy(:), ahz(:), ah(:)
    real(p_), intent(in) :: w(:) !weight at t_n
    real(p_), intent(out) :: w_new(:) !weight at t_n+dtao
    real(p_) :: gradient, t0, n0, norm, norm2, rate, x, dfdx, dfdE, E
    integer :: k, nm

    !    if(ns==3) goto 1000
    nm = nm_gk(ns)
    do k = 1, nm  !for Maxwellian distribution
       if (lost(k) .eqv. .true.) cycle
       x = radcor(k)
       !E = 0.5*mass_gk(ns)*(v(k)*vn_gk(ns))**2/kev
       t0 = gkt_func(x, ns)
       n0 = gkn_func(x, ns)
       gradient = -gkdndx_func(x,ns)/n0 - gkdtdx_func(x,ns)/t0 &
            & *(v(k)**2/(two*t0*kev/(mass_gk(ns)*vn_gk(ns)**2))-1.5_p_)
       norm = (charge_gk(ns)/qu)/(t0/tu)
       norm2 = norm*vn_gk(ns)/vu

       rate = xdrift1(k)*gradient &
            & - norm*(phix(k)*xdrift0(k) +phiy(k)*ydrift0(k) +phiz(k)*zdrift00(k))  &
            & - norm*(phix(k)*xdrift1(k) +phiy(k)*ydrift1(k) +phiz(k)*zdrift1(k))*gk_nonlinear(ns)  & 
            & + norm2*vpar(k)*(ahx(k)*xdrift0(k)+ ahy(k)*ydrift0(k) + ahz(k)*zdrift0(k)) &
            & + norm2*vpar(k)*(ahx(k)*xdrift1(k)+ ahy(k)*ydrift1(k) + ahz(k)*zdrift1(k))*gk_nonlinear(ns) & 
            & - norm2*(-mirror_force(k))*ah(k)

       w_new(k) = w(k) + rate* f0(x,v(k),ns)*ps_vol_gk(k,ns)/w_unit*dtao
    enddo

    !w_new(1:nm) = sum(w_new(1:nm))/nm !weight smoothing, for testing
    return

1000 do k = 1, nm ! for slowing-down distribution
       if (lost(k) .eqv. .true.) cycle
       x = radcor(k)
       E = 0.5*mass_gk(ns)*v(k)**2*vn_gk(ns)**2
       dfdx = slowing_down_xd(x, E)
       dfdE = slowing_down_Ed(x, E)
       norm = (charge_gk(ns)/qu)*tu*kev
       norm2 = norm*vn_gk(ns)/vu

       rate = - xdrift1(k) * dfdx &
            & + norm*(phix(k)*xdrift0(k)+phiy(k)*ydrift0(k)+phiz(k)*zdrift0(k))*dfdE  &
            & + norm*(phix(k)*xdrift1(k)+phiy(k)*ydrift1(k)+phiz(k)*zdrift1(k))*dfdE*gk_nonlinear(ns)  & 
            & - norm2*vpar(k)*(ahx(k)*xdrift0(k)+ ahy(k)*ydrift0(k) + ahz(k)*zdrift0(k))*dfdE &
            & - norm2*vpar(k)*(ahx(k)*xdrift1(k)+ ahy(k)*ydrift1(k) + ahz(k)*zdrift1(k))*dfdE*gk_nonlinear(ns) & 
            & + norm2*(-mirror_force(k))*ah(k)*dfdE
       w_new(k) = w(k) + rate*ps_vol_gk(k,ns)*vn_gk(ns)**3/w_unit*dtao

    enddo

  end subroutine push_gk_weight
end module gk_weight_pusher
module density_temperature_funcs
  use constants,only: p_
  implicit none
  real(p_),parameter :: a=0.6_p_ !meter, minor radius
  real(p_),parameter :: r0=a*0.5_p_ !meter, radial center of the simulation box
  real(p_),parameter :: ti_scale=0.3_p_ !dimension-less
  real(p_),parameter :: ni_scale=0.3_p_ !dimension-less
contains
  
 function ti_func(radcor) result (z) !in the unit of kev
    use constants,only: one
    use fk_module,only: ti0,kappa_ti
    use func_in_mc, only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=ti0*exp(-kappa_ti*a*ti_scale*tanh((r-r0)/(a*ti_scale)))
  end function ti_func

  
function ni_func(radcor) result (z) !in the SI unit
    use constants,only: one
    use fk_module,only: ni0,kappa_ni
    use func_in_mc, only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=ni0*exp(-kappa_ni*a*ni_scale*tanh((r-r0)/(a*ni_scale)))
  end function ni_func



 function kappa_ti_func(radcor) result (z)
    use constants,only: one
    use fk_module,only: kappa_ti
    use func_in_mc,only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=kappa_ti*(one-tanh((r-r0)/(a*ti_scale))**2)
  end function kappa_ti_func


  function kappa_ni_func(radcor) result (z)
    use constants,only: one
    use fk_module,only : kappa_ni
    use func_in_mc, only : minor_r_radcor !function
    real(p_),intent(in) :: radcor
    real(p_) :: r, z

    r=minor_r_radcor(radcor)
    z=kappa_ni*(one-tanh((r-r0)/(a*ni_scale))**2)
  end function kappa_ni_func

end module density_temperature_funcs


module push_ion_weight_module
contains
subroutine push_ion_weight(dtao,nmarker_i,active_i,radcor_i,theta_i,alpha_i, &
     & grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i,&
     & v_i, vpar_i,vx_i,vy_i,w_i,w_i_star)
  use constants,only:p_
  use constants,only: two,twopi,one_half,kev
  use normalizing,only: Ln
  use fk_module,only: vn_fk, ni0,mass_i,ps_vol_i,normalizing_factor
!  use force,only: field_perturbation_on_marker
  use func_in_mc,only: minor_r_prime
  use magnetic_coordinates, only : GSpsi_prime
  use density_temperature_funcs,only: kappa_ti_func, kappa_ni_func, ti_func
  implicit none
  integer,intent(in):: nmarker_i
  logical,intent(in):: active_i(nmarker_i)
  real(p_),intent(in):: dtao, radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i)
  real(p_),intent(in):: grad_psi_i(nmarker_i),grad_alpha_i(nmarker_i),grad_psi_dot_grad_alpha_i(nmarker_i),&
       & bval_i(nmarker_i), v_i(nmarker_i), vpar_i(nmarker_i),&
       & vx_i(nmarker_i),vy_i(nmarker_i),w_i(nmarker_i)
  real(p_),intent(out):: w_i_star(nmarker_i)
  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val
  real(p_):: tmp,gradient_term
  integer:: k

  !tmp=ti0*kev/(mass_i*vn_fk**2)
  !eq_particle_number=ni0*sqrt((mass_i/(twopi*ti0*kev))**3)*Ln**3*vn_fk**3/normalizing_factor !explicitly cancel the exp(-v^2/vt^2) dependence, valid only for the case that both physical equilibrium distribution and marker distribution are Maxwellian.

!$omp parallel do private(gradient_term,ef_par_val,ef_x_val,ef_y_val,mf_par_val,mf_x_val,mf_y_val,wiprime)
  do k=1,nmarker_i
     if( active_i(k).eqv..false.) then
        w_i_star(k)=0._p_
     else
        tmp=ti_func(radcor_i(k))*kev/(mass_i*vn_fk**2)
        eq_particle_number=ps_vol_i(k)*fi0(radcor_i(k),v_i(k))
        gradient_term=kappa_ni_func(radcor_i(k))+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti_func(radcor_i(k))
 !       call field_perturbation_on_marker(radcor_i(k),theta_i(k),alpha_i(k),active_i(k),&
        !            & ef_par_val,ef_x_val,ef_y_val)
        ef_x_val=0 !testing
        ef_y_val=0
        ef_par_val=0
        wiprime=eq_particle_number*twopi/tmp*(ef_par_val*vpar_i(k)+ef_x_val*vx_i(k)+ef_y_val*vy_i(k) ) & 
             &+eq_particle_number*gradient_term*GSpsi_prime*minor_r_prime(radcor_i(k))/bval_i(k)**2 &
             !& *(grad_psi_i(k)**2*grad_alpha_i(k)**2-grad_psi_dot_grad_alpha_i(k)**2)*ef_y_val &
             & *(bval_i(k)**2/GSpsi_prime**2)*ef_y_val !&
!em             &-eq_particle_number*gradient_term*minor_r_prime(radcor_i(k))/bval_i(k)*(vx_i(k)*mf_par_val&
!em             & -vpar_i(k)*(mf_x_val*grad_psi_i(k)**2+mf_y_val*grad_psi_dot_grad_alpha_i(k)))

        w_i_star(k)=w_i(k)+wiprime*dtao
        !if (isnan(wiprime)) write(*,'(a20,30(1pe12.1))') 'warning**nan appear', w_i_star(k),wiprime
          !write(*,*) ef_par_val
     endif
  enddo
!$omp end parallel do
!w_i_star=0._p_ !for testing
        !write(*,*) 'max wi=',maxval(abs(w_i)),'max wi_star', maxval(abs(w_i_star))
end subroutine push_ion_weight

subroutine push_ion_weight_using_electric_cylindrical_components(dtao,radcor_i,theta_i,alpha_i,&
     & r_i,z_i,phi_i,vr_i,vz_i,vphi_i,active_i,w_i,w_i_star,nmarker_i)
  use constants,only: p_, two,twopi,one_half,kev
  use normalizing,only: ln
  use fk_module, only : vn_fk
  use fk_module,only:ps_vol_i,mass_i,v_i,ni0 ,normalizing_factor
  use magnetic_coordinates,only:mpol,nrad,zgrid,xgrid, &
       & grad_psi_r,grad_psi_z
  use magnetic_field, only : br,bz,bphi
  use force
  use interpolate_module
  use func_in_mc,only: minor_r_prime
  use density_temperature_funcs,only: kappa_ti_func, kappa_ni_func, ti_func
  implicit none
  real(p_),intent(in):: dtao
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: radcor_i(nmarker_i),theta_i(nmarker_i),alpha_i(nmarker_i),w_i(nmarker_i)
  real(p_),intent(in):: r_i(nmarker_i),z_i(nmarker_i),phi_i(nmarker_i)
  real(p_),intent(in):: vr_i(nmarker_i),vz_i(nmarker_i),vphi_i(nmarker_i)
  logical,intent(in):: active_i(nmarker_i)
  real(p_),intent(out):: w_i_star(nmarker_i)

  real(p_):: eq_particle_number,wiprime
  real(p_):: ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val
  real(p_):: tmp,tmp2

  real(p_):: br_val,bz_val,bphi_val, grad_psi_r0,grad_psi_z0,bval
  integer:: k,ierr

  do k=1,nmarker_i
     v_i(k)=sqrt(vr_i(k)**2+vz_i(k)**2+vphi_i(k)**2)
  enddo
  !tmp=ti0*kev/(mass_i*vn_fk**2)

  !eq_particle_number=ni0*sqrt((mass_i/(twopi*ti0*kev))**3)*ln**3*vn_fk**3/normalizing_factor !explicitly cancel the exp(-v^2/vt^2) dependence, valid only for the case that both physical equilibrium distribution and marker distribution are Maxwellian.
  do k=1,nmarker_i
     if( active_i(k).eqv..false.) then
        w_i_star(k)=0._p_
     else
        tmp=ti_func(radcor_i(k))*kev/(mass_i*vn_fk**2)
        eq_particle_number=ps_vol_i(k)*fi0(radcor_i(k),v_i(k)) !general case
        tmp2=kappa_ni_func(radcor_i(k))+(v_i(k)**2/(two*tmp)-1.5_p_)*kappa_ti_func(radcor_i(k))

        bz_val=bz(r_i(k),z_i(k))
        br_val=br(r_i(k),z_i(k))
        bphi_val=bphi(r_i(k),z_i(k))
        bval=sqrt(bz_val*bz_val+br_val*br_val+bphi_val*bphi_val)

        call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
             & grad_psi_r,theta_i(k),radcor_i(k),grad_psi_r0)  !in magnetic coordinates
        call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
             & grad_psi_z, theta_i(k),radcor_i(k),grad_psi_z0)  !in magnetic coordinates
        call field_perturbation_on_marker2(radcor_i(k),theta_i(k),alpha_i(k), active_i(k),&
             & ef_cyl_r_val,ef_cyl_z_val,ef_cyl_phi_val)

        wiprime=eq_particle_number*twopi/tmp*(ef_cyl_r_val*vr_i(k)+ef_cyl_z_val*vz_i(k)+ef_cyl_phi_val*vphi_i(k)) &
             &+eq_particle_number*tmp2*minor_r_prime(radcor_i(k))/bval**2*&
             (grad_psi_z0*ef_cyl_r_val*bphi_val-grad_psi_r0*ef_cyl_z_val*bphi_val&
             & -grad_psi_z0*ef_cyl_phi_val*br_val+grad_psi_r0*ef_cyl_phi_val*bz_val)

        w_i_star(k)=w_i(k)+wiprime*dtao
     endif
  enddo
end subroutine push_ion_weight_using_electric_cylindrical_components


subroutine ion_velocity_and_metric_in_mc(nmarker_i, r_i , z_i,  &
     & radcor_i , theta_i , active_i,  &
     & vr_i , vphi_i , vz_i,  &
     & grad_psi_i , grad_alpha_i , grad_psi_dot_grad_alpha_i , bval_i,  & 
     & v_i, vpar_i , vx_i , vy_i ) !evaluate variables appearing on the rhs of the ion weight evolution equation, for linear case, these variables are independent of perturbations, and thus can be calculated once and stored, to be used multiple times in the iterations involved in the implicit delta-f method
  use constants,only: p_
  use constants,only: one
  use magnetic_coordinates,only: mpol,nrad,xgrid,zgrid,r_mc,z_mc, &
  & grad_psi,grad_alpha_r,grad_alpha_z,grad_alpha,grad_psi_dot_grad_alpha, &
  & grad_psi_r,grad_psi_z
use magnetic_field, only : br,bphi,bz
  use interpolate_module
  
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: r_i (nmarker_i), z_i (nmarker_i)
  real(p_),intent(in):: radcor_i (nmarker_i), theta_i (nmarker_i)
  logical,intent(in)::  active_i (nmarker_i)
  real(p_),intent(in):: vr_i (nmarker_i), vphi_i (nmarker_i), vz_i (nmarker_i)
  real(p_),intent(out):: grad_psi_i (nmarker_i), grad_alpha_i (nmarker_i), &
       & grad_psi_dot_grad_alpha_i (nmarker_i), bval_i (nmarker_i)
  real(p_),intent(out):: v_i (nmarker_i), vpar_i (nmarker_i), vx_i (nmarker_i), vy_i (nmarker_i)


  real(p_):: br_val,bphi_val,bz_val
  real(p_):: grad_psi_r_val,grad_psi_z_val,grad_alpha_r_val,grad_alpha_z_val
  integer:: k

  v_i (1:nmarker_i)=sqrt(vr_i (1:nmarker_i)**2+vz_i (1:nmarker_i)**2+vphi_i (1:nmarker_i)**2)

  do k=1, nmarker_i
     if(active_i (k).eqv..false.) cycle
     br_val=br(r_i (k),z_i (k))
     bphi_val=bphi(r_i (k),z_i (k))
     bz_val=bz(r_i (k),z_i (k))
     bval_i (k)=sqrt(br_val**2+bphi_val**2+bz_val**2)

     vpar_i (k)=(br_val*vr_i (k)+bz_val*vz_i (k)+bphi_val*vphi_i (k))/bval_i (k)

     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
          & grad_psi_r,theta_i (k),radcor_i (k),grad_psi_r_val)  !in magnetic coordinates

     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
          & grad_psi_z,theta_i (k),radcor_i (k),grad_psi_z_val)  !in magnetic coordinates

     vx_i (k)=(vr_i (k)*grad_psi_r_val+vz_i (k)*grad_psi_z_val)

     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
          & grad_alpha_r,theta_i (k),radcor_i (k),grad_alpha_r_val)  !in magnetic coordinates
     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
          & grad_alpha_z,theta_i (k),radcor_i (k),grad_alpha_z_val)  !in magnetic coordinates

     vy_i (k)=vphi_i (k)/r_i (k)+vr_i (k)*grad_alpha_r_val+vz_i (k)*grad_alpha_z_val

     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
          & grad_psi,theta_i (k),radcor_i (k),grad_psi_i (k))  !in magnetic coordinates

     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
          & grad_alpha,theta_i (k),radcor_i (k),grad_alpha_i (k))  !in magnetic coordinates

     call linear_2d_interpolate(mpol,nrad,zgrid,xgrid,&
          & grad_psi_dot_grad_alpha,theta_i (k),radcor_i (k),grad_psi_dot_grad_alpha_i (k))  !in magnetic coordinates
  enddo

  !for testing
!!$  do k=1,nmarker_i
!!$     if(active_i (k).eqv..true.) then
!!$        cos_theta=grad_psi_dot_grad_alpha_i(k)/(grad_psi_i(k)*grad_alpha_i(k))
!!$        sin_theta=sqrt(1-cos_theta**2) !sign?
!!$        theta1=atan((cos_theta-(vy_i(k)/vx_i(k)*grad_psi_i(k)/grad_alpha_i(k)))/sin_theta)
!!$        vper_i(k)=vx_i(k)/grad_psi_i(k)/cos(theta1)
!!$      write(*,*) sqrt(vpar_i(k)**2+vper_i(k)**2),v_i(k)  !tested, good agreement
!!$    endif
!!$  enddo

end subroutine ion_velocity_and_metric_in_mc

subroutine fk_push_first_step()
  use constants,only:two, one_half
  use fk_module,only: nmarker_i,dtao_fk, w_i, w_i_mid
  use fk_module,only: r_i, z_i, phi_i, r_i_old,z_i_old,phi_i_old, r_i_mid, z_i_mid, phi_i_mid
  use fk_module,only: vr_i_old,vz_i_old,vphi_i_old
  use fk_module,only: vr_i_integer, vz_i_integer, vphi_i_integer
  use fk_module,only: vr_i_mid, vz_i_mid, vphi_i_mid
  use fk_module,only: touch_bdry_i,active_i,radcor_i, theta_i, alpha_i
  use fk_module,only: touch_bdry_i_mid,active_i_mid,radcor_i_mid, theta_i_mid, alpha_i_mid
  use fk_module,only: grad_psi_i,grad_alpha_i,grad_psi_dot_grad_alpha_i,bval_i, v_i,vpar_i,vx_i,vy_i
  use fk_particle_coordinates_transform_module,only:compute_particle_magnetic_coordinates
  use sort_ions, only:sort_ions_according_to_poloidal_location
  implicit none
  
  r_i_old=r_i; z_i_old=z_i; phi_i_old=phi_i
  vr_i_old=vr_i_mid; vz_i_old=vz_i_mid;vphi_i_old=vphi_i_mid
  call push_ion_orbit_first_step(dtao_fk)
  vr_i_integer=(vr_i_old+vr_i_mid)/two; vz_i_integer=(vz_i_old+vz_i_mid)/two; vphi_i_integer=(vphi_i_old+vphi_i_mid)/two
!!$  call ion_velocity_and_metric_in_mc(nmarker_i, r_i_old, z_i_old,  &
!!$       & radcor_i, theta_i, active_i,  vr_i_integer, vphi_i_integer, vz_i_integer,  &
!!$       & grad_psi_i , grad_alpha_i , grad_psi_dot_grad_alpha_i , bval_i,  & 
!!$       & v_i, vpar_i , vx_i , vy_i) !output is used by ion weight pusher
!!$
!!$  call push_ion_weight(one_half*dtao_fk,nmarker_i,active_i,radcor_i,theta_i,alpha_i, &
!!$       & grad_psi_i , grad_alpha_i, grad_psi_dot_grad_alpha_i , bval_i,  & 
!!$       & v_i, vpar_i , vx_i, vy_i, w_i, w_i_mid) 
    call push_ion_weight_using_electric_cylindrical_components(one_half*dtao_fk,radcor_i,theta_i,alpha_i,&
             r_i_old,z_i_old,phi_i_old,vr_i_integer,vz_i_integer,vphi_i_integer,active_i,w_i,w_i_mid,nmarker_i)

  r_i_mid=(r_i+r_i_old)/two;  z_i_mid=(z_i+z_i_old)/two; phi_i_mid=(phi_i+phi_i_old)/two
  !prepare ion magnetic coordinates so that markers can be sorted and deposited in the coordinates
  call compute_particle_magnetic_coordinates(nmarker_i, r_i_mid, phi_i_mid, z_i_mid, &
       & radcor_i_mid, active_i_mid, touch_bdry_i_mid, theta_i_mid, alpha_i_mid)
  call sort_ions_according_to_poloidal_location(theta_i_mid) !prepare for ion weight pusher and deposition
end subroutine fk_push_first_step


subroutine fk_push_second_step()
  use fk_module,only: dtao_fk, nmarker_i, r_i_mid, z_i_mid, phi_i_mid,&
       & radcor_i_mid , theta_i_mid , alpha_i_mid, active_i_mid,  vr_i_mid , vphi_i_mid , vz_i_mid,  &
       & grad_psi_i_mid , grad_alpha_i_mid , grad_psi_dot_grad_alpha_i_mid , bval_i_mid,  & 
       & v_i_mid, vpar_i_mid , vx_i_mid , vy_i_mid
  use fk_module,only: w_i, w_i_star !output
  use fk_module,only: r_i,phi_i,z_i, radcor_i,theta_i,alpha_i,active_i,touch_bdry_i
  use fk_particle_coordinates_transform_module,only:compute_particle_magnetic_coordinates
  use sort_ions, only:sort_ions_according_to_poloidal_location
  implicit none

!!$  call ion_velocity_and_metric_in_mc(nmarker_i, r_i_mid, z_i_mid,  &
!!$       & radcor_i_mid , theta_i_mid , active_i_mid,  vr_i_mid , vphi_i_mid , vz_i_mid,  &
!!$       & grad_psi_i_mid , grad_alpha_i_mid , grad_psi_dot_grad_alpha_i_mid , bval_i_mid,  & 
!!$       & v_i_mid, vpar_i_mid , vx_i_mid , vy_i_mid) !output is used by ion weight pusher
!!$  call push_ion_weight(dtao_fk,nmarker_i,active_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid, &
!!$       & grad_psi_i_mid , grad_alpha_i_mid, grad_psi_dot_grad_alpha_i_mid , bval_i_mid,  & 
!!$       & v_i_mid, vpar_i_mid , vx_i_mid , vy_i_mid,w_i,w_i_star) 
   call push_ion_weight_using_electric_cylindrical_components(dtao_fk,radcor_i_mid,theta_i_mid,alpha_i_mid,&
             r_i_mid,z_i_mid,phi_i_mid,vr_i_mid,vz_i_mid,vphi_i_mid,active_i_mid,w_i,w_i_star,nmarker_i)

  w_i(1:nmarker_i)=w_i_star(1:nmarker_i)
  call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,&
       & radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !prepare ion magnetic coordinates so that markers can be sorted and deposited in the magnetic coordinates
  call sort_ions_according_to_poloidal_location(theta_i) !prepare for ion weight pusher and deposition

end subroutine fk_push_second_step

subroutine push_ion_orbit_first_step(dtao_fk) !wrapper of "push_full_orbit_cylindrical_boris_with_additional_output" subroutine
  use constants,only:p_
  use constants,only: twopi
  use fk_module,only: nmarker_i,touch_bdry_i
  use fk_module,only: r_i,phi_i,z_i,vr_i,vphi_i,vz_i !as input and output
!  use fk_module,only: phi_i_mid
  use fk_module,only: vr_i_mid,vphi_i_mid,vz_i_mid !as output
  use fk_module,only: radcor_i,theta_i,alpha_i,active_i
  use boris
  implicit none
  real(p_),intent(in):: dtao_fk
  integer:: k

  !$omp parallel do
  do k=1,nmarker_i
     if(touch_bdry_i(k).eqv..true.) cycle
     !call push_full_orbit_cylindrical_boris(dtao_fk,r_i(k),phi_i(k),z_i(k),vr_i(k),vphi_i(k),vz_i(k))
     call push_full_orbit_cylindrical_boris2(dtao_fk,radcor_i(k),theta_i(k),alpha_i(k), active_i(k), &
          & r_i(k),phi_i(k),z_i(k),vr_i(k),vphi_i(k),vz_i(k),&
          & vr_i_mid(k),vphi_i_mid(k),vz_i_mid(k)) !obtain velocity at t_{n+1/2}, outputting the projections of this velocity onto both the basis vectors at t_{n+1/2} and those at t_{n+1}
!!$     !{r,z,phi}_i_mid is known before entering this subroutine ({r,z,phi}_i_mid is ethier the initial codintion of the second-pusher or the output of it.)
  enddo
  !$omp end parallel do
end subroutine push_ion_orbit_first_step

!!$subroutine push_ion_orbit_second_step(dtao_fk) 
!!$  !calculate the velocity at t_{n+1}, giving the projections of this velocity onto (1) the basis vectors at t_{n+1} (to prepare for the deposition process) and (2) the basis vectors at t_{n+3/2} (to prepare input for the next secod-pusher)
!!$  !output the location at t_{n+3/2} (to prepare input for the next secod-pusher)
!!$  !(the spatial location at t_{n+1} is already computed by the first boris-pusher)
!!$  use constants,only:p_
!!$  use constants,only: twopi
!!$  use fk_module,only: nmarker_i,touch_bdry_i_mid
!!$  use fk_module,only: r_i_mid,phi_i_mid,z_i_mid !as input (t_{n+1/2}) and output (t_{n+3/2})
!!$  use fk_module,only: vr_i_integer_mid,vphi_i_integer_mid,vz_i_integer_mid !input (projection of v at t_{n} onto basis vector at t_{n+1/2}) and output (projection of v at t_{n+1} onto basis vector at t_{n+3/2})
!!$  use fk_module,only: vr_i_integer,vphi_i_integer,vz_i_integer !as output, projection of velocity at t_{n+1} onto the basis vectors at t_{n+1}
!!$  use fk_module,only: phi_i !as input, location at t_{n+1}
!!$  use domain_decomposition,only: myid
!!$  use boris
!!$  implicit none
!!$  real(p_),intent(in):: dtao_fk
!!$  integer:: k
!!$
!!$  !$omp parallel do
!!$  do k=1,nmarker_i
!!$     if(touch_bdry_i_mid(k).eqv..true.) cycle
!!$     call push_full_orbit_cylindrical_boris2(dtao_fk,r_i_mid(k),phi_i_mid(k),z_i_mid(k),&
!!$          & vr_i_integer_mid(k),vphi_i_integer_mid(k),vz_i_integer_mid(k),&
!!$          & phi_i(k),vr_i_integer(k),vphi_i_integer(k),vz_i_integer(k))
!!$ enddo
!!$  !$omp end parallel do
!!$!if(myid.eq.0) write(*,*) z_i_mid(20), dtao_fk*vz_i_integer(20)
!!$end subroutine push_ion_orbit_second_step

function fi0(radcor,v) result (z) ! v in unit of vn, fi0 in unit 1/(Ln**3*vn**3)
  use constants,only: p_
  use constants,only: two,twopi,kev
  use normalizing,only: Ln
  use fk_module,only: mass_i, vn_fk
  use density_temperature_funcs,only: ni_func, ti_func
  implicit none
  real(p_),intent(in)  :: radcor,v
  real(p_)  :: z
  real(p_)  :: v_si, ni, ti !local variables
  v_si=v*vn_fk
  ti=ti_func(radcor)
  ni=ni_func(radcor)
  z=ni*sqrt((mass_i/(twopi*ti*kev))**3)*exp(-mass_i*v_si**2/(two*ti*kev))
  z=z*(vn_fk**3*Ln**3)
end function fi0

end module push_ion_weight_module
module diagnosis_mod
contains

  subroutine compute_particle_and_heat_flux(time, ns, lost, mu_gk, vpar_gk, w_gk,&
       &  xgc, zgc, xdrift1, xdrift0)
    use constants, only: p_, two, kev, three
    use magnetic_coordinates, only: radial_width, dradcor, dtheta, mpol, toroidal_range, xgrid, &
         & jacobian, vol, grad_psi
    use func_in_mc, only: b_mc_func, grad_psi_func
    use gk_module,only: mass_gk, charge_gk, nm_gk, vn_gk, w_unit, nsm
    use gk_profile_funcs, only : gkt_func, gkn_func, gkdtdx_func
    use radial_module, only : baxis, minor_a
    use domain_decomposition, only: myid
    use mpi
    implicit none
    real(p_), intent(in) :: time
    integer, intent(in) :: ns
    logical, intent(in) :: lost(:)
    real(p_), intent(in) :: mu_gk(:), vpar_gk(:), w_gk(:), xgc(:), zgc(:), xdrift1(:), xdrift0(:)
    integer, parameter :: nsub = 20
    real(p_) :: tmp_xa(nsub), dv(nsub), myheat_flux(nsub), heat_flux(nsub), myheat_flux0(nsub), heat_flux0(nsub)
    real(p_) :: myptcl_flux(nsub), ptcl_flux(nsub), mydensity_perturbation(nsub), density_perturbation(nsub)
    real(p_) :: myptcl_flux0(nsub), ptcl_flux0(nsub)
    real(p_) :: diffusivity(nsub), gyro_bohm(nsub)
    real(p_) :: x, dx, kinetic, gradient, trash
    integer :: nm, i, j, k, jeq, ierr
    character(len=64) :: fn1, fn2,  fn1b, fn2b, fn3
    logical, save :: is_first = .true.
    integer, allocatable, save :: u1(:), u2(:), u1b(:), u2b(:), u3(:)

    nm = nm_gk(ns)

    if ((is_first .eqv. .true.) .and. (myid==0)) then 
       is_first=.false.
       allocate(u1(nsm), u1b(nsm))
       allocate(u2(nsm), u2b(nsm))
       allocate(u3(nsm))       
       do i = 1, nsm
          fn1 = 'heat_flux_nsx.txt'; fn1b = 'heat_flux0_nsx.txt'
          fn2 = 'ptcl_flux_nsx.txt'; fn2b = 'ptcl_flux0_nsx.txt'
          fn3 = 'dens_pert_nsx.txt'
          write(fn1(13:13),'(i1.1)') i
          write(fn2(13:13),'(i1.1)') i
          write(fn1b(14:14),'(i1.1)') i
          write(fn2b(14:14),'(i1.1)') i
          write(fn3(13:13),'(i1.1)') i
          open(newunit=u1(i), file=fn1, position="append")
          open(newunit=u2(i), file=fn2, position="append")
          open(newunit=u1b(i), file=fn1b, position="append")
          open(newunit=u2b(i), file=fn2b, position="append")
          open(newunit=u3(i), file=fn3, position="append")
       enddo
    endif

    dx = radial_width/nsub
    do j = 1, nsub
       tmp_xa(j) = xgrid(1) + dx*(j-1)
    enddo

    myheat_flux = 0._p_; myheat_flux0 = 0._p_
    myptcl_flux = 0._p_; myptcl_flux0 = 0._p_
    density_perturbation = 0._p_
    do k = 1, nm
       if(lost(k) .eqv. .true.) cycle
       j = floor((xgc(k)- tmp_xa(1))/dx)+1
       j = min(j, nsub)
       j = max(j, 1)
       kinetic = vpar_gk(k)**2/two + mu_gk(k)*b_mc_func(zgc(k), xgc(k))
       myheat_flux(j) = myheat_flux(j) + w_gk(k)*kinetic*xdrift1(k)/grad_psi_func(zgc(k), xgc(k))
       myptcl_flux(j) = myptcl_flux(j) + w_gk(k)        *xdrift1(k)/grad_psi_func(zgc(k), xgc(k))
       myheat_flux0(j) = myheat_flux0(j) + w_gk(k)*kinetic*xdrift0(k)/grad_psi_func(zgc(k), xgc(k))
       myptcl_flux0(j) = myptcl_flux0(j) + w_gk(k)        *xdrift0(k)/grad_psi_func(zgc(k), xgc(k))
       mydensity_perturbation(j) = mydensity_perturbation(j) + w_gk(k)
    enddo

    call MPI_Reduce(myheat_flux, heat_flux, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(myptcl_flux, ptcl_flux, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(myheat_flux0, heat_flux0, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(myptcl_flux0, ptcl_flux0, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(mydensity_perturbation, density_perturbation, nsub, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)    

    if(myid==0) then
       do j = 1, nsub
          jeq = 1 + (j-1)*int(dx/dradcor)
          x =  tmp_xa(j)
          dv(j) = sum(abs(jacobian(1:mpol-1, jeq)))*dx*dtheta*toroidal_range
          ptcl_flux(j) = ptcl_flux(j)*w_unit            *vn_gk(ns)   /dv(j) 
          heat_flux(j) = heat_flux(j)*w_unit*mass_gk(ns)*vn_gk(ns)**3/dv(j) !in unit m^2/s
          ptcl_flux0(j) = ptcl_flux0(j)*w_unit            *vn_gk(ns)   /dv(j) 
          heat_flux0(j) = heat_flux0(j)*w_unit*mass_gk(ns)*vn_gk(ns)**3/dv(j) !in unit m^2/s

          density_perturbation(j) = density_perturbation(j)*w_unit/dv(j) 
          gyro_bohm(j) = sqrt(mass_gk(ns))*sqrt(gkt_func(x,ns)*kev)**3/(minor_a*Baxis**2*charge_gk(ns)**2) !in unit m^2/s
          gradient = -gkn_func(x,ns)*gkdtdx_func(x,ns)*kev*sum(grad_psi(:,jeq))/mpol !assume constant density
          diffusivity(j) = 2./3.0*heat_flux(j)/gradient
       enddo
       write(u1(ns),'(50(1pe20.8))') time, sum(diffusivity)/nsub, sum(diffusivity/gyro_bohm)/nsub, &
            & sum(heat_flux)/nsub, (heat_flux(j), j=1, nsub)
       write(u2(ns),'(50(1pe20.8))') time, (ptcl_flux(j), j=1, nsub)
       write(u1b(ns),'(50(1pe20.8))') time, trash, trash, &
            & sum(heat_flux0)/nsub, (heat_flux0(j), j=1, nsub)
       write(u2b(ns),'(50(1pe20.8))') time, (ptcl_flux0(j), j=1, nsub)

       write(u3(ns),'(50(1pe20.8))') time, (density_perturbation(j), j=1, nsub)
    endif

  end subroutine compute_particle_and_heat_flux
 
  
subroutine count_lost_markers_gk(ns, radcor)  !diagnosis
  use constants, only : p_
    use gk_module, only: nm_gk
    use domain_decomposition, only: myid
    use magnetic_coordinates, only: xlow, xupp
    implicit none
    integer, intent(in) :: ns
    real(p_), intent(in) :: radcor(:)
    integer :: nlost, k

    nlost = 0
    do k = 1, nm_gk(ns)
       if( (radcor(k) >= xupp) .or.  (radcor(k) <= xlow) ) nlost = nlost + 1
    enddo
    write(*,'(A50, 20i15)') 'myid, species, total_makers, lost_markers =', myid, ns, nm_gk(ns), nlost
  end subroutine count_lost_markers_gk

  
subroutine parallel_current_diagnostic(niter,iter)
use fk_module,only: ni0,charge_i, vn_fk
use gk_module,only: ngk0
!use perturbation_field,only: jpar_i_old, jpar_e_old

integer,intent(in):: niter,iter

!if(iter.eq.niter) write(*,*) iter, jpar_e_old(5,5)*vn_gk*ngk0*(-charge_i),jpar_i_old(5,5)*vn_fk*ni0*charge_i,&
!     & jpar_e_old(5,5)*vn_gk*ngk0*(-charge_i)/(jpar_i_old(5,5)*vn_fk*ni0*charge_i)

 !write(*,*) jpar_e_old(15,15)*vn_gk*ngk0*(-charge_i),jpar_i_old(15,15)*vn_fk*ni0*charge_i,&
  !   & jpar_e_old(15,15)*vn_gk*ngk0*(-charge_i)/(jpar_i_old(15,15)*vn_fk*ni0*charge_i)

end subroutine parallel_current_diagnostic

  
subroutine check_rms_ion_weight(w_i_star,nmarker_i)
  use constants,only: p_
  use domain_decomposition,only: myid,numprocs
  use mpi
  implicit none
  integer,intent(in):: nmarker_i
  real(p_),intent(in):: w_i_star(nmarker_i)
  real(p_):: my_sum,all_sum
  integer:: ierr
  all_sum=0.
  my_sum=sum(w_i_star(1:nmarker_i)**2)/nmarker_i
  call MPI_Reduce(my_sum, all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD,ierr);
  if(myid.eq.0) write(*,*) 'rms_of_weight_averaged_over_all_particles=',sqrt(all_sum)/numprocs !assume that number of parcles in each process is equal to each other
end subroutine check_rms_ion_weight


subroutine draw_grids_on_theta_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !on theta=constant surface, the subroutine name is wrong, not necessarily top view
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nrad
  real(p_),intent(in)::tor_shift_mc(mpol,nrad),r_mc(mpol,nrad),z_mc(mpol,nrad)
  real(p_):: phi,alpha,dalpha
  integer:: i,j,itor,mtor

  i=10
  i=20
  i=60
  i=1
  !  i=mpol
  mtor=20
  dalpha=twopi/mtor

  open(113,file='grids_on_theta_isosurface.txt')
  do itor=1,mtor+1
     !  do itor=1,2
     alpha=dalpha*(itor-1)
     !     do j=1,nrad,10
     !    do j=1,nrad,1
     do j=1,1
        phi=alpha+tor_shift_mc(i,j) !phi is changing due to the radial dependene of tor_shift
        write(113,*) phi,r_mc(i,j),z_mc(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

  open(113,file='grids_on_theta_isosurface2.txt')
  do j=1,nrad,10
     do itor=1,mtor+1
        alpha=dalpha*(itor-1)
        phi=alpha+tor_shift_mc(i,j)
        write(113,*) phi,r_mc(i,j),z_mc(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)
end subroutine draw_grids_on_theta_isosurface


subroutine draw_alpha_isosurface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !on a nonzero radial range
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nrad
  real(p_),intent(in)::tor_shift_mc(mpol,nrad),r_mc(mpol,nrad),z_mc(mpol,nrad)
  real(p_):: phi,alpha0
  integer:: i,j

  alpha0=twopi/8._p_
  open(113,file='alpha_isosurface.txt')
  do j=1,101,5
     do i=1,mpol
        phi=alpha0+tor_shift_mc(i,j) 
        !phi=alpha0
        write(113,*) phi,r_mc(i,j),z_mc(i,j)
     enddo
     write(113,*)
     write(113,*)
  enddo
  close(113)

end subroutine draw_alpha_isosurface


subroutine draw_alpha_contours_on_a_magnetic_surface(mpol,nrad,tor_shift_mc,r_mc,z_mc) !field lines on a magnetic surface
  use constants,only: p_
  use constants,only: two,twopi
  !  use magnetic_coordinates,only:mtor
  implicit none
  integer,intent(in):: mpol,nrad
  real(p_),intent(in)::tor_shift_mc(mpol,nrad),r_mc(mpol,nrad),z_mc(mpol,nrad)
  real(p_):: phi,alpha0,tor_range
  integer:: i,j,ialpha,ishift,u,iphi
  integer,parameter::nalpha=10,nphi=20

  open(newunit=u,file='alpha_contours_on_magnetic_surface.txt')
  !do j=1,30,2
  j=nrad/2 !select a radial location, i.e., a magnetic surface 
  tor_range=twopi/4
  do ialpha=1,nalpha
     !alpha0=0._p_+tor_range/(nalpha-1)*(ialpha-1)
     alpha0=0._p_
     do i=1,mpol
        phi=alpha0+tor_shift_mc(i,j)
!!$        if(phi<0. .or. phi>tor_range) then !shift into the range [0:tor_range]
!!$           ishift=floor(phi/tor_range)
!!$           phi=phi-ishift*tor_range 
!!$           alpha0=alpha0-ishift*tor_range
!!$           write(u,*)
!!$           write(u,*)
!!$        endif
        !phi=alpha0
        write(u,*) phi,z_mc(i,j),r_mc(i,j)
     enddo
     write(u,*)
     write(u,*)
  enddo
  close(u)
  open(newunit=u,file='ref_surface')
  do iphi=1,nphi
     phi=0.+twopi/(nphi-1)*(iphi-1)
     do i=1,mpol
        write(u,*) phi,z_mc(i,j),r_mc(i,j)
     enddo
     write(u,*)
     write(u,*)
  enddo
  close(u)
phi=alpha0+tor_shift_mc(1,j)

call field_line_tracing_simplified(r_mc(1,j),z_mc(1,j),phi)
end subroutine draw_alpha_contours_on_a_magnetic_surface


subroutine calculate_possibility_density(v,total_number,sample_number,starting_value,ending_value)
  use constants,only:p_
  use  constants,only:one,two,four,five,twopi,eight
  implicit none

  integer,intent(in):: total_number,sample_number
  real(p_),intent(in):: v(total_number)
  real(p_),intent(in):: starting_value,ending_value
  real(p_):: xcenter(sample_number-1)
  real(p_):: possibility_density(sample_number-1)
  real(p_):: bin_npt(sample_number-1)
  real(p_):: interval
  real(p_):: x(sample_number)
  integer:: i,j

  interval=(ending_value-starting_value)/(sample_number-1)
  do i=1,sample_number
     x(i)=starting_value+interval*(i-1)
  enddo

  bin_npt=0
  do i=1,sample_number-1
     do j=1,total_number
        if (v(j).ge.x(i) .and. v(j).le.x(i+1)) bin_npt(i)=bin_npt(i)+1
     enddo
  enddo

  do i=1,sample_number-1
     xcenter(i)=(x(i)+x(i+1))/two
     possibility_density(i)=bin_npt(i)/(total_number)/interval
     write(*,*) xcenter(i), possibility_density(i)
  enddo

end subroutine calculate_possibility_density


subroutine report(t)
  use magnetic_coordinates,only: mtor,nrad
!  use perturbation_field, only:ex=>ex_left,ey=>ey_left,epar=>epar_left

!  use perturbation_field,only:source_e1,source_e2 !,jper_x_i_left,jper_y_i_left
!  use perturbation_field,only:source1,source2,source3,source_faraday1, source_faraday2,source_faraday3
  use constants,only:p_

  real(p_),intent(in):: t
  integer:: i,j

i=mtor/2
j=nrad/2

!write(*,'(7(1pe14.4))') t,ex(i,j),ey(i,j),epar(i,j),mf_x(i,j),mf_y(i,j),mf_par(i,j)
!write(*,*) jper_x_i_left(mtor/2,nrad/2),jper_y_i_left(mtor/2,nrad/2)!,source_e1(mtor/2,nrad/2),source_e2(mtor/2,nrad/2)
!write(*,*) source1(mtor/3,nrad/3),source2(mtor/3,nrad/3),source3(mtor/3,nrad/3)
!write(*,*) source_faraday1(mtor/3,nrad/3),source_faraday2(mtor/3,nrad/3),source_faraday3(mtor/3,nrad/3)
end subroutine report


subroutine bperp_perturbation(t, u)
  use constants, only: p_, kev
  use normalizing, only: tu, qu, vu
  use perturbation_field, only: apara
  use magnetic_coordinates, only: nrad, mtor, grad_psi_dot_grad_alpha, grad_psi, grad_alpha
  use domain_decomposition, only: ipol_eq
  use derivatives_in_xyz, only: x_derivative, y_derivative
  implicit none
  real(p_), intent(in) :: t
  integer, intent(in) :: u
  real(p_) :: bperp(nrad), ax(mtor, nrad), ay(mtor, nrad), gx(nrad), gy(nrad), gxdgy(nrad)
  integer :: i
  
  call x_derivative(apara(:,:,1), ax(:,:))
  call y_derivative(apara(:,:,1), ay(:,:))
  ax = ax*tu*kev/(qu*vu) ! change to SI unit
  ay = ay*tu*kev/(qu*vu)
  gx(:) = grad_psi(ipol_eq, :)
  gy(:) = grad_alpha(ipol_eq, :)
  gxdgy(:) = grad_psi_dot_grad_alpha(ipol_eq, :)
  i = mtor/2
  bperp = sqrt(ax(i,:)**2*gx(:)**2 + ay(i,:)**2*gy(:)**2 + 2*ax(i,:)*ay(i,:)*gxdgy(:))
  write(u, '(3000(1pe18.4))') t, bperp(1:nrad:2)
end subroutine bperp_perturbation


subroutine mode_evolution(t, s, m, n, file_unit) 
  use constants, only: p_
  implicit none
  real(p_), intent(in) :: t
  integer, intent(in) :: m, n
  real(p_), intent(in) :: s(0:m-1, n)
  integer, intent(in) :: file_unit
  integer :: i0, j0, j

  i0 = m/2 !toroidal location
  j0 = n/2 !radial location
  !write(file_unit,'(3000(1pe18.4))') t, (s(i0, j),  j = 1, n, 2)
  write(file_unit,'(3000(1pe18.4))') t, s(i0, j0)

end subroutine mode_evolution

subroutine mode_evolution2(t) 
  use constants,only: one
  use constants,only:p_
  use magnetic_coordinates,only: m=>mtor,n=>nrad
  use perturbation_field,only: ef_cyl_phi_left,ef_cyl_r_left,ef_cyl_z_left
  use perturbation_field,only: ef_cyl_phi_right,ef_cyl_r_right,ef_cyl_z_right
  use transform_module,only: twod_fourier_transform
  implicit none
  real(p_),intent(in):: t
  real(p_):: a(0:m-1,0:n-1)
  complex(p_):: a_fft(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative

  do i=0,m-1
     do j=0,n-1
        a(i,j)=ef_cyl_phi_left(i+1,j+1)
     enddo
  enddo
  call twod_fourier_transform(a,a_fft,m,n)

  ipositive=1
  inegative=m-ipositive

  write(*,'(20(1pe14.4))') t,(real(a_fft(ipositive,j)),imag(a_fft(ipositive,j)),j=0,4),&
       & ef_cyl_phi_left(m/2,n/2)

end subroutine mode_evolution2


subroutine mode_evolution3(t,a,m,n,file_unit) 
  use constants,only: one
  use constants,only:p_
  use transform_module,only: twod_fourier_transform,dst_dft
  use fourn_module,only: twod_fourier_transform_nr
  implicit none
  integer,intent(in):: m,n
  real(p_),intent(in):: a(0:m-1,0:n-1)
  real(p_),intent(in):: t
  integer,intent(in):: file_unit
  complex(p_):: a_spectrum(0:m-1,0:n-1)
  integer:: i,j
  integer:: ipositive,inegative

!m=size(a,1)
!n=size(a,2)
!call dst_dft(a,a_spectrum,m,n) !using DST for radial direction and DFT for toroidal direction
!a_spectrum=a_spectrum/((n+1)*m) !the corresponding expansion coeficient is the dst_dft devided by (n+1)*m, see my notes on Fourier analysis
call twod_fourier_transform(a,a_spectrum,m,n) 
!call twod_fourier_transform_nr(a,a_spectrum,m,n) 
a_spectrum=a_spectrum/(n*m) 

  ipositive=1
  !inegative=m-ipositive

  write(file_unit,'(20(1pe20.8))') t,(real(a_spectrum(ipositive,j)),imag(a_spectrum(ipositive,j)),j=0,3)

end subroutine mode_evolution3


subroutine mode_evolution4(t, s, m, n, file_unit) 
  use constants,only: one, p_
  use transform_module
  use control_parameters, only : nh_min
  !  use fourn_module,only: twod_fourier_transform_nr
  implicit none
  integer, intent(in) :: m,n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  real(p_), intent(in) :: t
  integer, intent(in) :: file_unit
  complex(p_) :: spectrum(0:m-1, 0:n-1)
  integer :: ntor, kr

  !  call oned_sine_transform2(s,a_dst,m,n) 
  !  call oned_fourier_transform1(a_dst,a_spectrum,m,n) 

  call dst_dft(s, spectrum, m, n) !using DST for radial direction and DFT for toroidal direction
  spectrum = spectrum/(2*(n+1)*m) !the expansion coeficient is the dst_dft devided by (n+1)*m, see my notes on Fourier analysis

  ntor = nh_min
  write(file_unit,'(20ES18.4)') t, (real(spectrum(ntor,kr)), imag(spectrum(ntor,kr)), kr=0,3)

end subroutine mode_evolution4

subroutine mode_evolution5(t, s, m, n, file_unit)
  use constants,only:p_
  use, intrinsic :: iso_c_binding
  use my_FFTW3, only: plan_toroidal, in1, out1
  use control_parameters, only : nh_min, nh_max
  implicit none
  include 'fftw3.f03'
  real(p_), intent(in) :: t
  integer, intent(in) :: m,n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  integer, intent(in) :: file_unit
  integer :: i, jrad

!  call oned_fourier_transform1(s,s_spectrum,m,n)                          

  jrad=154

  in1(:) = s(:,jrad) !copy in, meanwhile convert real array to complex array                                      
  call fftw_execute_dft(plan_toroidal, in1(:), out1(:))

  write(file_unit,'(2000(1pe20.8))') t, (real(out1(i))/m, imag(out1(i))/m, i = nh_min, nh_max)
end subroutine mode_evolution5


subroutine mode_evolution6(t, s, m, n, file_unit) 
  use constants, only: p_
  use, intrinsic :: iso_c_binding
  use my_FFTW3, only: plan_toroidal, in1, out1
  use control_parameters, only: nh_min, nh_max
  implicit none
  include 'fftw3.f03'
  real(p_), intent(in) :: t
  integer, intent(in) :: m, n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  integer, intent(in) :: file_unit
  integer :: i, j0
  complex(p_) :: spectrum(0:m-1)

  j0 = n/2 !choose a radial location

  in1(:) = s(:,j0) !copy in, meanwhile convert real array to complex array
  call fftw_execute_dft(plan_toroidal, in1(:), out1(:))
  spectrum(:) = out1(:)/m

  write(file_unit,'(900(1pe18.4))') t, (real(spectrum(i)), imag(spectrum(i)), i=nh_min, nh_max)

end subroutine mode_evolution6

subroutine mode_evolution7(t, s, m, n, file_unit) 
  use constants, only: p_
  use, intrinsic :: iso_c_binding
  use my_FFTW3, only: plan_toroidal, in1, out1
  use control_parameters, only: nh_min, nh_max
  implicit none
  include 'fftw3.f03'
  real(p_), intent(in) :: t
  integer, intent(in) :: m, n
  real(p_), intent(in) :: s(0:m-1,0:n-1)
  integer, intent(in) :: file_unit
  integer :: i, j
  complex(p_) :: spectrum(0:m-1)

  j = n/2
  in1(:) = s(:,j) !copy in, meanwhile convert real array to complex array
  call fftw_execute_dft(plan_toroidal, in1(:), out1(:))
  spectrum(:) = out1(:)/m

  write(file_unit,'(2000(1pe18.4))') t, (real(spectrum(i)), imag(spectrum(i)), i = nh_min, nh_max)

end subroutine mode_evolution7



end module diagnosis_mod


module mode_structure
  implicit none
contains
  subroutine mode_structure_in_xy_plane(kt,GCLR,a,partial_file_name)
    use constants,only:p_
    use magnetic_coordinates,only:xgrid,ygrid
!    use domain_decomposition,only:myid
    integer,intent(in)::kt,GCLR
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in)::partial_file_name
    character(100):: full_file_name
    integer:: i,j,m,n,u

    full_file_name='ms/xy_polxxx_txxxxxx'//partial_file_name
    write(full_file_name(10:12),'(i3.3)') GCLR
    write(full_file_name(15:20),'(i6.6)') kt

    open(newunit=u,file=full_file_name)
    m=size(a,1)
    n=size(a,2)
    do j=1,n
       do i=1,m
          write(u,*) xgrid(j), ygrid(i), a(i,j)
       enddo
       write(u,*)
    enddo
    close(u)
  end subroutine mode_structure_in_xy_plane

  subroutine mode_structure_in_xz_plane(kt, a, partial_file_name)
    use constants,only:p_
    use constants,only:pi
    use magnetic_coordinates,only: xgrid,zgrid,dtheta,mtor
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR,dtheta2
    use constants,only:twopi
    use mpi
    implicit none
    integer,intent(in)::kt
    real(p_),intent(in):: a(:,:)
    character(len=*),intent(in)::partial_file_name
    integer:: itor,ipol,j,ierr,m,n 
    real(p_):: my_a_xz_plane(size(a,2)),a_xz_plane(size(a,2),0:numprocs/ntube-1)
    real(p_)::my_theta
    character(100)::full_file_name
    integer:: u !file unit number

    m=size(a,1)
    n=size(a,2)

    itor=mtor/2 !choose a alpha (i.e., y) grid
    do j=1,n
       my_a_xz_plane(j)=a(itor,j)
    enddo
    call MPI_gather(my_a_xz_plane, n, MPI_real8, &
                  & a_xz_plane,    n, MPI_real8, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       full_file_name='ms/xz_txxxxxx'//partial_file_name
       write(full_file_name(8:13),'(i6.6)') kt
       open(newunit=u,file=full_file_name)
       do ipol = 0, numprocs/ntube-1 !poloidal direction
          my_theta=-pi+ipol*dtheta2
          do j=1,n !radial direction
             write(u,*) xgrid(j), my_theta,a_xz_plane(j,ipol)
          enddo
          write(u,*)
       enddo
       close(u)
    endif
  end subroutine mode_structure_in_xz_plane


  subroutine mode_structure_in_yz_plane(kt, ayxz, partial_file_name)
    use constants, only: p_, pi, twopi
    use magnetic_coordinates, only: mtor, ygrid, zgrid, nrad, mpol2, &
         & r_mc, z_mc, tor_shift_mc
    use domain_decomposition, only: myid, tube_comm, multi_eq_cells, gclr
    use mpi
    implicit none
    integer, intent(in) :: kt
    real(p_), intent(in) :: ayxz(:, :, :)
    character(len=*), intent(in) :: partial_file_name
    real(p_), allocatable :: ay(:), ayz(:,:)
    real(p_) :: phi
    character(100) :: full_file_name
    integer :: itor, ipol, ipol_eq, jrad, ierr, m, n, u
    INTEGER :: status(MPI_STATUS_SIZE)
    
    m = size(ayxz, 1) !toroidal
    n = size(ayxz, 2) !radial
    allocate(ay(m))
    allocate(ayz(m, 0:mpol2))
    jrad = nrad/2 !choose a radial index
    
    ay(:) = ayxz(:, jrad, 1)
    call MPI_gather(ay,  m, MPI_real8, &
         &          ayz, m, MPI_real8, 0, tube_COMM, ierr)

    ! z boundary
    if(gclr==mpol2-1) call mpi_send(ayxz(:,jrad, 2), m, MPI_real8, 0, 123, tube_comm, ierr)
    if(gclr==0) call mpi_recv(ayz(:, mpol2), m, MPI_real8, mpol2-1, 123, tube_comm, status, ierr)

    if(myid==0) then
       full_file_name='ms/yz_txxxxxx'//partial_file_name
       write(full_file_name(8:13),'(i6.6)') kt
       open(newunit=u, file=full_file_name)
       do ipol = 0, mpol2 !z : poloidal
          ipol_eq=multi_eq_cells*ipol + 1 !index in the equilibrium grids
          do itor =1, mtor !y : toroidal 
             phi = ygrid(itor) + tor_shift_mc(ipol_eq, jrad)
             write(u,'(19ES16.5E3)') ygrid(itor), zgrid(ipol_eq), ayz(itor,ipol), &
                  & r_mc(ipol_eq,jrad), z_mc(ipol_eq,jrad), phi
          enddo
          write(u,*) 
          !if(ipol.eq.GCLR_cut) write(u,*) !to inform gnuplot that this is a new data block, to aviud line connection between the following data and previous data
       enddo
       close(u)
    endif
  end subroutine mode_structure_in_yz_plane


  subroutine mode_structure_in_poloidal_plane(kt, a, str)
    use constants, only: p_, one
    use domain_decomposition, only: theta_start, dtheta2
    implicit none
    integer, intent(in) :: kt
    character(*), intent(in) :: str
    real(p_), intent(in) :: a(:, :, :)
    integer, parameter :: nz = 40 !along field line
    real(p_) :: theta(0:nz-1), c
    real(p_) :: perturb(size(a,1), size(a,2), 0:nz-1)
    real(p_) :: a0(size(a,1), size(a,2), 0:nz-1), a1(size(a,1),size(a,2), 0:nz-1)
    character(len=100) :: file_name
    integer :: m, n, j, i

    m = size(a,1) !toroidal
    n = size(a,2) !radial

    do i = 0, nz-1 !interpolate to more z gridpoints 
       theta(i) = theta_start + dtheta2/nz*i
       c = (theta(i) - theta_start)/dtheta2
       perturb(:,:,i) = a(:,:,1)*(one-c) + a(:,:,2)*c
    enddo

    do i = 0, nz-1
       do j = 1, n
          a0(:, j, i) = sum(perturb(:,j, i))/m
       enddo
    enddo

    file_name='ms/poloidal_plane_txxxxxx'//str//'_neq0'
    write(file_name(20:25),'(i6.6)') kt
    call to_poloidal_plane(a0, m, n, nz, theta, file_name)

    a1 = perturb - a0
    file_name='ms/poloidal_plane_txxxxxx'//str//'_nneq0'
    write(file_name(20:25),'(i6.6)') kt
    call to_poloidal_plane(a1, m, n, nz, theta, file_name)

  end subroutine mode_structure_in_poloidal_plane

  subroutine to_poloidal_plane(yxz, ny, nx, nz, theta, file_name)
    use constants, only: p_, pi, twopi
    use magnetic_coordinates, only: ygrid, r_mc, z_mc, tor_shift_mc, mpol, mpol2, nsegment, &
         & zgrid, pfn, tfn, toroidal_range
    use domain_decomposition, only: myid, tube_comm, dtheta2
    use interpolate_module
    use math, only: shift_toroidal
    use mpi
    implicit none
    integer, intent(in) :: ny, nx, nz
    real(p_), intent(in) :: yxz(ny, nx, 0:nz-1), theta(0:nz-1)
    character(len=*), intent(in) :: file_name
    real(p_), allocatable :: field(:,:), field0(:,:)
    real(p_) :: tmp(ny+1)
    real(p_) :: phi, my_theta, alpha, r, z, th, tor_shift
    integer :: j, i, iz,  u, ierr

    allocate(field(nx, 0:nz-1))
    allocate(field0(nx, 0:nz*mpol2-1))

    !choose a cylindrical toroidal angle, for which the mode structure is computed
    phi = 0.5_p_*twopi/nsegment

    do j = 1, nx !radial
       do i = 0, nz-1 !local theta range
          call linear_1d_interpolate(mpol, zgrid, tor_shift_mc(:,j), theta(i), tor_shift)
          alpha = phi - tor_shift ! alpha of the fixed phi
          call shift_toroidal(alpha, toroidal_range)
          !Interpolate in alpha, to the same cylindrical toroidal angle:
          tmp(1:ny) = yxz(:, j, i)
          tmp(ny+1) = yxz(1, j, i) !periodic boundary condition
          call linear_1d_interpolate(ny+1, ygrid, tmp, alpha, field(j, i)) 
       enddo
    enddo

    call MPI_gather(field, nx*nz, MPI_real8, field0, nx*nz, MPI_real8, 0, tube_COMM, ierr)

    if(myid == 0) then
       open(newunit=u, file=file_name)
       do j = 1, nx
          do i = 0, mpol2 - 1
             my_theta = -pi + i*dtheta2
             do iz = 0, nz-1
                th = my_theta + dtheta2/nz*iz
                call linear_1d_interpolate(mpol, zgrid, r_mc(:,j), th, r)
                call linear_1d_interpolate(mpol, zgrid, z_mc(:,j), th, z)
                write(u,'(19ES16.5E3)') r, z, field0(j, i*nz+iz), pfn(j), th, tfn(j)
             enddo
          enddo
       enddo
       close(u)
    endif
  end subroutine to_poloidal_plane

end module mode_structure


module spectrum_diagnostic
  use constants,only:p_
  implicit none
  save
  real(p_):: spectrum_amplitude_old
contains
  subroutine spectrum_diagnostic_routine(iter,t,a,file_unit)
    use constants,only:pi
    use magnetic_coordinates,only: xgrid,zgrid,dtheta,mtor
    use domain_decomposition,only:numprocs,myid,ntube,tube_comm,GCLR
    use transform_module,only: twod_fourier_transform_xz,oned_fourier_transform1
    use constants,only:twopi
    use mpi
    integer,intent(in):: iter
    real(p_),intent(in):: a(:,:),t
    integer,intent(in):: file_unit
    complex(p_):: a_dft(size(a,1),size(a,2))
    complex(p_):: my_a_xz_plane(size(a,2)),a_xz_plane(size(a,2),numprocs/ntube)
    complex(p_):: a_xz_plane_fft(size(a,2),numprocs/ntube)
    integer:: itor,ipol,j,ierr,m,n,mz
    real(p_):: spectrum_amplitude

    m=size(a,1)
    n=size(a,2)
    mz=numprocs/ntube
!!$    itor=mtor/2 !choose a y grid
!!$    do j=1,n
!!$       my_a_xz_plane(j)=a(itor,j)
!!$    enddo
    call oned_fourier_transform1(a,a_dft,m,n) !Fourier transform along toroidal direction
    do j=1,n
       my_a_xz_plane(j)=a_dft(2,j) !select the fundamental harmonic of the toroidal expansion
    enddo

    call MPI_gather(my_a_xz_plane, n, MPI_complex16, &
         & a_xz_plane,    n, MPI_complex16, 0, tube_COMM, ierr)
    if(myid.eq.0) then
       call twod_fourier_transform_xz(a_xz_plane,a_xz_plane_fft,n,mz)
       write(file_unit,'(20(1pe20.8))') t,(real(a_xz_plane(j,j)),imag(a_xz_plane(j,j)),j=1,3)
!!$       spectrum_amplitude=abs(a_xz_plane_fft(4,4))
!!$       write(*,*) iter,'relative_error_in_low_kx_low_kz_spectrum_amplitude=',&
!!$            & abs(spectrum_amplitude_old-spectrum_amplitude)/spectrum_amplitude
!!$       spectrum_amplitude_old=spectrum_amplitude !prepare for the next convergence checking
    endif

  end subroutine spectrum_diagnostic_routine
end module spectrum_diagnostic


subroutine visualize_grid()
  use constants, only : p_, twopi, zero
  use magnetic_coordinates, only : r_mc, z_mc, mpol,nrad, mtor, ygrid, tor_shift_mc, toroidal_range
  use magnetic_field, only : psi_func, qfunc0
  use math, only : shift_toroidal
  implicit none
  integer, parameter :: npt=2000
  integer :: u,i,j,k
  real(p_) :: phi, phi0
  real(p_) :: rf(npt), zf(npt), phif(npt), dphi, qval

  open(newunit=u,file='grid.txt')
  do k=1, mtor
     do i=1, mpol
        do j=1,nrad
           phi = ygrid(k) + tor_shift_mc(i,j)
           !call shift_toroidal(phi,toroidal_range) 
           write(u,*) r_mc(i,j), z_mc(i,j), phi
        enddo
     enddo
  enddo
  close(u)

  i = (mpol+1)/2
  j = nrad/2
  !phi0 = ygrid(1) + tor_shift_mc(i,j)
  phi0 = 0
  qval = qfunc0(psi_func(r_mc(i,j), z_mc(i,j)))
  dphi = twopi/2*abs(qval)/(npt-1)
  call field_line_tracing0(r_mc(i,j), z_mc(i,j), phi0, npt, dphi, rf,zf,phif)
  open(newunit=u,file='grid_z.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  call field_line_tracing0(r_mc(i,j),z_mc(i,j), phi0, npt, -dphi, rf,zf,phif)
  open(newunit=u,file='grid_z2.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)

  !-------------------------------
  j = 1
  phi0=ygrid(1) + tor_shift_mc(i,j)
  qval = qfunc0(psi_func(r_mc(i,j), z_mc(i,j)))
  dphi = twopi*abs(qval)/(npt-1)/2
  call field_line_tracing0(r_mc(i,j), z_mc(i,j), phi0, npt, dphi, rf,zf,phif)
  open(newunit=u,file='grid_z3.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  call field_line_tracing0(r_mc(i,j),z_mc(i,j), phi0, npt, -dphi, rf,zf,phif)
  open(newunit=u,file='grid_z4.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  !-------------------------------
  j = nrad
  phi0=ygrid(1) + tor_shift_mc(i,j)
  qval = qfunc0(psi_func(r_mc(i,j), z_mc(i,j)))
  dphi = twopi*abs(qval)/(npt-1)/2
  call field_line_tracing0(r_mc(i,j), z_mc(i,j), phi0, npt, dphi, rf,zf,phif)
  open(newunit=u,file='grid_z5.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  call field_line_tracing0(r_mc(i,j),z_mc(i,j), phi0, npt, -dphi, rf,zf,phif)
  open(newunit=u,file='grid_z6.txt')
  do k=1,npt
     write(u,*) rf(k), zf(k), phif(k)
  enddo
  close(u)
  
  open(newunit=u,file='inner_bdry.txt')
  do i =1, mpol
     write(u,*) r_mc(i, 1), z_mc(i, 1)
  enddo
  close(u)

  open(newunit=u,file='outer_bdry.txt')
  do i =1, mpol
     write(u,*)  r_mc(i, nrad), z_mc(i, nrad)
  enddo
  close(u)

end subroutine visualize_grid


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
program main
  use constants, only: p_, twopi, pi, kev, two, mu0, one_half, zero
  use control_parameters, only: kstart, kend, dt_omega_i_axis, iplot_mode_structure, &
       &  store_restart_data, fk_switch, diagnosis, adiabatic_electrons
  use normalizing, only: bn,ln, dtao_main, qu, tu, vu, nu
  use magnetic_coordinates, only: mpol,nrad, xlow, xupp, nsegment, dtheta, vol, grad_psi, &
       & mpol2, mtor, zgrid, tor_shift_mc, toroidal_range, GSpsi_prime
  use magnetic_field, only: qfunc, pfn_func
  use func_in_mc, only: minor_r_radcor
  use table_in_mc, only: prepare_table_in_mc
  use radial_module, only: baxis, r_axis, minor_a, radcor_fixed, j_fixed, qpsi
  use fk_module, only: mass_i,charge_i,dtao_fk, nmarker_i, ni0,ti0,kappa_ti, &
       & nmarker_i_per_cell, r_i,z_i,phi_i, radcor_i,theta_i,alpha_i,tor_shift_i,ps_vol_i,&
       & touch_bdry_i,active_i, touch_bdry_i_mid,active_i_mid, &
       & radcor_i_mid,theta_i_mid,alpha_i_mid, &
       & vr_i,vz_i,vphi_i, ntouch_bdry_i, total_ntouch_bdry_i, &
       & vr_i_integer,vz_i_integer,vphi_i_integer,&
       & vr_i_mid,vphi_i_mid,vz_i_mid, &
       & r_i_mid,z_i_mid,phi_i_mid,vr_i_integer_mid,vz_i_integer_mid,vphi_i_integer_mid, &
       & grad_psi_i_mid,grad_alpha_i_mid,grad_psi_dot_grad_alpha_i_mid,bval_i_mid, &
       & v_i_mid,vpar_i_mid,vx_i_mid,vy_i_mid, &
       & w_i,w_i_star,w_i_mid, initialize_fk, tn_fk, vn_fk, omegan_fk
  use gk_module, only: nsm, nmmax, nm_gk, mass_gk, charge_gk, tgk0, ngk0, &
       & dtao_gk, vn_gk, w_unit, xgc, zgc, ygc, vpar_gk, w_gk, w_gk_mid, &
       & xgc_mid, zgc_mid, ygc_mid, vpar_gk_mid, mu_gk, v_gk, ps_vol_gk, lost_gc, &
       & x_ring, y_ring, z_ring, initialize_gk, sort_gk_markers
  use load_gk_mod, only: load_gk
  use gyro_ring_mod, only: set_gyro_phase, gyro_ring
  use gyro_average_mod, only: gyro_average
  use drift, only: compute_drift
  use gk_trajectory_pusher, only: push_gc !, count_lost_markers_gk
  use gk_weight_pusher, only: push_gk_weight
  use perturbation_field, only: allocate_field_matrix, potential, phix, phiy, phiz, &
       &  apara, apara_h, apara_s, apara_s_old, ax, ay, az, ahx, ahy, ahz
  use domain_decomposition, only: myid,numprocs, nvp, tube_comm,grid_comm,ntube,gclr,tclr, &
       & dtheta2,theta_start,my_right,my_left, my_right2, my_left2, multi_eq_cells, ipol_eq, dvol
  use misc, only: calculate_dvol

  use fk_particle_coordinates_transform_module, only: compute_particle_magnetic_coordinates, &
       & clean_up_lost_markers_fk, count_lost_markers_fk
  use mode_structure
  use diagnosis_mod, only: report, mode_evolution, mode_evolution6, bperp_perturbation, &
       & compute_particle_and_heat_flux
  use deposit_fk_module, only: deposit_fk
  use deposit_gk_module, only: deposit_gk
  use poisson, only: prepare_poisson_matrix, solve_poisson
  use ampere, only: prepare_ampere_matrix, solve_ampere, apara_s_evolution, &
       & apara_resplit_and_weight_pullback
  use push_ion_weight_module
  use initial_half_step_for_boris
  use interpolate_module, only: linear_2d_interpolate
  use restart, only: write_data_for_restarting, read_data_for_restarting
  use my_FFTW3
  use spectrum_diagnostic
  use mpi
  implicit none

  integer :: ierr, id_writing_evolution, k, kt, ns
  integer :: file_unit1, file_unit2, file_unit3, file_unit4, file_unit5
  real(p_) :: omega_i_axis, dt_second, rho_i, k_binormal1, k_binormal2, beta_ni
  real(p_) :: minor_r_min, minor_r_max, minor_r_width, va0
  real(p_) :: t1, t2, tarray(2) !store the cpu clock time
  real(p_), allocatable :: xdrift0(:), ydrift0(:), zdrift0(:), zdrift00(:), mirror_force(:)
  real(p_), allocatable :: xdrift1(:), ydrift1(:), zdrift1(:)
  real(p_), allocatable :: phix_ga(:), phiy_ga(:), phiz_ga(:), ax_ga(:), ay_ga(:), az_ga(:)
  real(p_), allocatable :: ahx_ga(:), ahy_ga(:), ahz_ga(:), ah_ga(:) !gyro-averaged perturbations
  real(p_), dimension(:,:), allocatable :: density_left, density_right !gk density
  real(p_), dimension(:,:), allocatable :: jpar_left, jpar_right !gk parallel current
  real(p_), allocatable :: cs(:), vthermal(:)
  character(len=9) :: tmp

  CALL MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  if(myid==0) write(*,*) 'numprocs=', numprocs, 'myid=', myid
  call cpu_time(tarray(1))  !f95 intrinsic subroutine

  call read_parameters()
  nvp = numprocs/ntube
  GCLR = INT(myid/ntube)
  TCLR = MOD(myid,ntube)
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, GCLR, TCLR, GRID_COMM, ierr)
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, TCLR, GCLR, TUBE_COMM, ierr)

  t1 = mpi_wtime() !measure the wall time

  if(mod(numprocs, ntube) .ne. 0) then
     write(*,*) "mod(numprocs,ntube) must be zero, please adjust numprocs or ntube"
     goto 1234 !end the job
  endif
  mpol2 = numprocs/ntube !poloidal grids for perturbed field
  if(mod(mpol-1, mpol2) .ne. 0) then
     write(*,*) 'error: Mod(mpol-1, numprocs/ntube) must be zero, please adjust poloidal gridpoint number in the input namelist'
     goto 1234 !end the job
  endif
  dtheta2=twopi/mpol2 !the poloidal angle spacing of grids for perturbations
  my_right = GCLR+1
  my_right2 = GCLR+2
  if(my_right==mpol2) my_right=0
  if(my_right2==mpol2) my_right2=0
  if(my_right2==mpol2+1) my_right2=1

  my_left=GCLR-1
  my_left2=GCLR-2
  if(my_left==-1) my_left=mpol2-1
  if(my_left2==-1) my_left2=mpol2-1
  if(my_left2==-2) my_left2=mpol2-2
  !Domain decomposition.
  !A mpi process is responsible for the poloidal range [theta_start : theta_start+dtheta2]
  theta_start = -pi+GCLR*dtheta2 

  call read_and_process_equilibrium() !in cylindrical coordinates
  call construct_magnetic_coordinates()

  !dtheta is the poloidal angle spacing of the equilibrium grids,
  multi_eq_cells = NINT(dtheta2/dtheta)
  !equilibrium poloidal index of the present MPI processes
  ipol_eq = 1+nint((theta_start-zgrid(1))/dtheta)
  if(myid==0) write(*,*) 'multi_eq_cells=', multi_eq_cells
  call calculate_dvol(multi_eq_cells, dvol)
  call mapping_table_for_cylindrical_to_magnetic_coordinates()
  !if(myid.eq.0) call field_lines_analyse()
  call prepare_table_in_mc() 

  if(myid==0 .and. diagnosis .eqv. .true.) call visualize_grid()
  minor_r_min = minor_r_radcor(xlow)
  minor_r_max = minor_r_radcor(xupp)
  minor_r_width = minor_r_max-minor_r_min
  if(myid==0) write(*,*) 'minor_r_min, minor_r_max:',  minor_r_min, minor_r_max
  if(myid==0) write(*,*) 'minor_r_width:', minor_r_width
  if(myid==0) write(*,*) 'minor_r center (m):', (minor_r_min+minor_r_max)/two

  call initialize_gk(nsm, baxis, minor_a, dt_omega_i_axis) !set radial profiles etc.
  w_unit = ngk0(1)*vol/nm_gk(1) !gk marker weight unit
  if(myid==0) write(*,*) 'w_unit=', w_unit
  qu = charge_gk(1)
  tu = tgk0(1)
  vu = vn_gk(1)
  nu = ngk0(1)
  dtao_main = dtao_gk(1)
  allocate(vthermal(nsm))
  vthermal(:) = sqrt(2*tgk0(:)*kev/mass_gk(:))
  if(myid==0) write(*,'(A30, 10ES16.4)') 'velocity unit (m/s) =', vu
  if(myid==0) write(*,'(A30, 10ES16.4)') 'vu/vthermal=', vu/vthermal
  if(myid==0) write(*,*) 'mtor, nrad=', mtor, nrad

  call allocate_field_matrix(mtor, nrad)

  if(kstart == 0) then
     potential = 0; phix = 0; phiy = 0; phiz = 0
     apara = 0; ax = 0; ay = 0; az = 0
     apara_s = 0; apara_h = 0; ahx = 0; ahy = 0; ahz = 0
     do ns = 1, nsm
        call load_gk(ns, nm_gk(ns), xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns), &
             & mu_gk(:,ns), v_gk(:,ns), ps_vol_gk(:,ns), w_gk(:,ns))
     enddo
  else if(kstart > 0) then
     call read_data_for_restarting(kstart)
  else
     stop "Error: kstart should >= zero"
  endif

  call set_gyro_phase()

  if(fk_switch==1) then
     call initialize_fk()
     omegan_fk = bn*charge_i/mass_i !cyclotron angular frequency in Hz
     tn_fk = twopi/omegan_fk !time unit used in this program
     vn_fk = Ln/tn_fk !the value of the normalizing velocity in SI unit m/s
     beta_ni = ni0*mass_i*vn_fk**2/two/(bn**2/(two*mu0)) !following the notation used in my notes
     call load_fk()
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,radcor_i,active_i,touch_bdry_i,theta_i,alpha_i)
  endif

  omega_i_axis = abs(baxis)*charge_gk(1)/mass_gk(1)
  rho_i = sqrt(tgk0(1)*kev/mass_gk(1))/omega_i_axis
  dtao_fk = dt_omega_i_axis/omega_i_axis/tn_fk !time step in unit of tn_fk
  dt_second = dt_omega_i_axis/omega_i_axis
  if(myid==0) write(*,*) 'dt (seconds)=', dt_second
  if(myid==0) write(*,'(A30, 10ES16.4)') 'cycltron angular frequency (MHz): ', baxis*charge_gk/mass_gk/10**6
  if(myid==0) write(*,'(A30, 10ES16.4)') 'Parallel CFL condition: dz/(dt*vt)', r_axis*qpsi(1)*(twopi/mpol2)/(dt_second*vthermal)
  if(myid==0) write(*,*) 'first radial sine harmonic: kr*rho_i=',pi/minor_r_width*rho_i
  if(myid==0) write(*,*) 'R0/rho_i=',r_axis/rho_i
  !if(myid==0) write(*,*) radcor_minor_r(0.5d0)
  k_binormal1=nsegment*abs(qfunc(radcor_fixed))/minor_r_radcor(radcor_fixed)*rho_i
  k_binormal2=nsegment*abs(baxis)/(GSpsi_prime*&
       & (grad_psi(1,j_fixed)+grad_psi(mpol/2,j_fixed))/two)*rho_i
  if(myid==0) write(*,*) 'k_binorm*rhoi1=',k_binormal1,'k_binorm*rhoi2=',k_binormal2
  if(myid==0) write(*,*) 'number of radial harmonics that should be included (pi*shear0*k_binormal/kr1)=',&
       & pi*0.84*k_binormal1/(pi/minor_r_width*rho_i)
  !if(myid.eq.0) write(*,*) 'omega_star1/twopi (kHz)=', k_binormal1*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_
  !if(myid.eq.0) write(*,*) 'omega_star2/twopi (kHz)=', k_binormal2*kappa_ti*rho_i*omega_i_axis/twopi/1000._p_
  if(myid==0) write(*,'(A25, 10ES16.4)') 'Species beta at axis=', tgk0*kev*ngk0/(baxis**2/(two*mu0))
  if(myid==0) write(*,'(A20, 10ES16.4)') 'tgk0 (keV)=', tgk0(:)
  if(myid==0) write(*,'(A20, 10ES16.4)') 'ngk0 (m^-3)=', ngk0(:)
  va0 = abs(baxis)/sqrt(mu0*mass_gk(1)*ngk0(1))
  allocate(cs(nsm))
  cs = sqrt(tgk0*kev/mass_gk)
  if(myid==0) write(*,*) 'VA0 (10^6m/s)=',va0/(1.d6)
  if(myid==0) write(*,'(A30, 10ES16.4)') 'sound speed (10^6m/s)=',cs/(1.d6)
  if(myid==0) write(*,'(A30, 10ES16.4)') 'VA0/Cs=',va0/cs
  if(myid==0) write(*,'(A30, 10ES16.4)') 'R0/Cs (seconds) =', r_axis/cs(:)
  if(myid==0) write(*,*) 'compressional Alfven wave freq/omega_i=k_binorma*va0/omega_axis=',k_binormal1/rho_i*va0/omega_i_axis

  if(kstart==0 .and. fk_switch==1) then
     vr_i_integer = vr_i !vr_i_integer is the projection of velocity at t_{n} to the basis vectors at t_{n}
     vz_i_integer = vz_i 
     vphi_i_integer = vphi_i 

     do k = 1, nmarker_i
        !vr_i is initially the the projection of velocity at t_{0} to basis vectors at t_{0}, after the backward pushing, vr_i is the projection of velocity at t_{-1/2} to basis vector at t_{0}
        call backward_half_step_for_boris(sign(1._p_,charge_i),dtao_fk,r_i(k),z_i(k),phi_i(k),vr_i(k),vz_i(k),vphi_i(k)) !push only velocity, to set initial condition for the first step of boris algorithm, using mutiple steps (2nd RK)
!!$        call forward_half_step_for_boris(sign(1._p_,charge_i),dtao_fk,r_i(k),z_i(k),phi_i(k),&
!!$             & vr_i_integer(k),vz_i_integer(k),vphi_i_integer(k),&
!!$             & r_i_mid(k),z_i_mid(k),phi_i_mid(k),vr_i_integer_mid(k),vz_i_integer_mid(k),vphi_i_integer_mid(k)) !push location half-step and then the velocity is projected onto the new local vector basis
        w_i(k) = ps_vol_i(k)*1d-15 !initial weight of markers
     enddo

     call compute_particle_magnetic_coordinates(nmarker_i, r_i, phi_i, z_i, radcor_i, active_i, &
          & touch_bdry_i, theta_i, alpha_i) !so that markers can be sorted and deposited in magnetic coordinates
     call push_ion_orbit_first_step(dtao_fk)
  endif

  call initialize_fft()
  if(kstart > 0 .and. fk_switch==1)  then
     call compute_particle_magnetic_coordinates(nmarker_i,r_i,phi_i,z_i,&
          & radcor_i,active_i,touch_bdry_i,theta_i,alpha_i) !theta of lost markers are not calculated
     call clean_up_lost_markers_fk() !drop lost markers so that we we have smaller particle arrays
  endif

  id_writing_evolution = numprocs/2 !corresponding to theta=0
  if(myid == id_writing_evolution) then
     write(tmp,'(i9.9)') kstart
     open(newunit=file_unit1, file='phi_evolution'//tmp//'.txt')
     open(newunit=file_unit2, file='apara_evolution'//tmp//'.txt')
     open(newunit=file_unit3, file='phi_n_evolution'//tmp//'.txt')
     open(newunit=file_unit4, file='apara_n_evolution'//tmp//'.txt')
  endif

  if(myid == 0) then
     write(tmp,'(i9.9)') kstart
     open(newunit=file_unit5, file='bperp_evolution'//tmp//'.txt')
  endif

  t2=mpi_wtime(); if(myid.eq.0) write(*,*) 'wall_time (seconds) used before calculating polarization matrix=', t2-t1
  call prepare_poisson_matrix()
  call prepare_ampere_matrix()
  t2=mpi_wtime(); if(myid.eq.0) write(*,*) 'wall_time (seconds) used before main time loop=', t2-t1

  allocate(phix_ga(nmmax), phiy_ga(nmmax), phiz_ga(nmmax)) ! _ga stands for "gyro-average"
  allocate(ax_ga(nmmax), ay_ga(nmmax), az_ga(nmmax))
  allocate(ahx_ga(nmmax), ahy_ga(nmmax), ahz_ga(nmmax))
  allocate(ah_ga(nmmax)) 
  allocate(xdrift0(nmmax), ydrift0(nmmax), zdrift0(nmmax))
  allocate(xdrift1(nmmax), ydrift1(nmmax), zdrift1(nmmax))
  allocate(zdrift00(nmmax))
  allocate(mirror_force(nmmax))

  allocate(density_left(mtor,nrad), density_right(mtor,nrad)) 
  allocate(jpar_left(mtor,nrad), jpar_right(mtor,nrad))

  do ns = 1, nsm
     lost_gc(:,ns) = xgc(:,ns)>xupp .or. xgc(:,ns)<xlow
     call gyro_ring(ns, lost_gc(:,ns), xgc(:,ns), ygc(:,ns), zgc(:,ns), &
          & mu_gk(:,ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns)) 
  enddo

  do kt = kstart+1, kend ! time-advancing
     if(myid==0 .and. mod(kt-1, 100)==0) write(*,*) 'time-step No.', kt
     !---- first step of the 2nd order Runge-Kutta-------
     if(fk_switch==1) then
        call fk_push_first_step()
        call deposit_fk(nmarker_i, active_i_mid,radcor_i_mid,theta_i_mid,alpha_i_mid,phi_i_mid,w_i_mid)
     endif
     density_left = zero;  density_right = zero !before the deposition, set the array to zero
     jpar_left = zero; jpar_right = zero
     do ns = 1, nsm ! loop over gk species
        call gyro_average(ns, lost_gc(:,ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns), &
             & phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga)
        call compute_drift(ns, lost_gc(:,ns), xgc(:,ns), zgc(:,ns), mu_gk(:,ns), vpar_gk(:,ns), &
             & phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, &
             & xdrift0, zdrift0, ydrift0, mirror_force, zdrift00, xdrift1, zdrift1, ydrift1)
        if(mod(kt-1, 100)==0)  call compute_particle_and_heat_flux((kt-1)*dt_second, ns, lost_gc(:, ns), &
             & mu_gk(:,ns), vpar_gk(:,ns), w_gk(:,ns), xgc(:,ns), zgc(:,ns), xdrift1(:), xdrift0(:))

        call push_gk_weight(ns, lost_gc(:,ns), one_half*dtao_gk(ns), xgc(:,ns),&
             & vpar_gk(:,ns), v_gk(:,ns), xdrift0, ydrift0, zdrift0, zdrift00, &
             & mirror_force, xdrift1, ydrift1, zdrift1, &
             & phix_ga, phiy_ga, phiz_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga, w_gk(:,ns), w_gk_mid(:,ns))
        call push_gc(ns, lost_gc(:,ns), one_half*dtao_gk(ns), xdrift0, zdrift0, ydrift0, mirror_force, &
             & xdrift1, zdrift1, ydrift1, xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns), &
             & w_gk_mid(:,ns), xgc_mid(:,ns), zgc_mid(:,ns), ygc_mid(:,ns), vpar_gk_mid(:,ns))

        lost_gc(:,ns) = lost_gc(:,ns) .or. xgc_mid(:,ns)>xupp .or. xgc_mid(:,ns)<xlow
        call sort_gk_markers(ns, zgc_mid(:,ns), step=1)

        call gyro_ring(ns, lost_gc(:,ns), xgc_mid(:,ns), ygc_mid(:,ns),  zgc_mid(:,ns), &
             & mu_gk(:,ns),  x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns))
        call deposit_gk(ns, lost_gc(:,ns), vpar_gk_mid(:,ns), w_gk_mid(:,ns), &
             & x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns),  &
             & density_left, density_right, jpar_left, jpar_right)
     enddo

     call solve_poisson(density_left, density_right, potential, phix, phiy, phiz)
     if(adiabatic_electrons .eqv. .false.) then ! EM case, so Ampere's law needs to be solved.
        apara_s_old = apara_s
        call apara_s_evolution(apara_s_old(:,:,1), apara_s(:,:,1), phiz, 0.5_p_*dtao_main)
        call solve_ampere(1, jpar_left, jpar_right, apara_s, apara_h, apara, ax, ay, az, ahx, ahy, ahz)
     endif
     !---------------second step of 2nd order Runge-Kutta------------------------
     if(fk_switch==1) then
        call fk_push_second_step()
        call deposit_fk(nmarker_i, active_i, radcor_i, theta_i, alpha_i, phi_i, w_i_star)
     endif
     density_left = zero; density_right = zero
     jpar_left = zero; jpar_right = zero
     do ns = 1, nsm
        call gyro_average(ns, lost_gc(:,ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns), &
             & phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga)
        call compute_drift(ns, lost_gc(:,ns), xgc_mid(:,ns), zgc_mid(:,ns), mu_gk(:,ns), vpar_gk_mid(:,ns), &
             &  phix_ga, phiy_ga, phiz_ga, ax_ga, ay_ga, az_ga, &
             & xdrift0, zdrift0, ydrift0, mirror_force, zdrift00, xdrift1, zdrift1, ydrift1)

        call push_gk_weight(ns, lost_gc(:,ns), dtao_gk(ns),  xgc_mid(:,ns), &
             & vpar_gk_mid(:,ns), v_gk(:,ns), xdrift0, ydrift0, zdrift0, zdrift00, &
             & mirror_force, xdrift1, ydrift1, zdrift1, &
             & phix_ga, phiy_ga, phiz_ga, ahx_ga, ahy_ga, ahz_ga, ah_ga, w_gk(:,ns), w_gk(:,ns))
        call push_gc(ns, lost_gc(:,ns), dtao_gk(ns), xdrift0, zdrift0, ydrift0, mirror_force, &
             & xdrift1, zdrift1, ydrift1, xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns), &
             & w_gk(:,ns), xgc(:,ns), zgc(:,ns), ygc(:,ns), vpar_gk(:,ns))

        lost_gc(:,ns) = lost_gc(:,ns) .or. xgc(:,ns)>xupp .or. xgc(:,ns)<xlow 
        call sort_gk_markers(ns, zgc(:,ns), step=2)

        call gyro_ring(ns, lost_gc(:,ns), xgc(:,ns), ygc(:,ns), zgc(:,ns), &
             & mu_gk(:,ns), x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns))
        call deposit_gk(ns, lost_gc(:,ns), vpar_gk(:,ns), w_gk(:,ns), &
             & x_ring(:,:,ns), y_ring(:,:,ns), z_ring(:,:,ns), &
             density_left, density_right,  jpar_left, jpar_right)
        !if(mod(kt,1000)==0) call count_lost_markers_gk(ns)
     enddo

     call solve_poisson(density_left, density_right, potential, phix, phiy, phiz)
     if(adiabatic_electrons .eqv. .false.) then
        call apara_s_evolution(apara_s_old(:,:,1), apara_s(:,:,1), phiz, dtao_main)
        call solve_ampere(2, jpar_left, jpar_right, apara_s, apara_h, apara, ax, ay, az, ahx, ahy, ahz)
        call apara_resplit_and_weight_pullback(w_gk, apara_s, apara_h, ahx, ahy, ahz) !at the end of each time-step
     endif

     ! Particle pusher and field solver finish one full time step. The following are diagnosis:
     if(myid==id_writing_evolution .and. mod(kt, 50)==0) then
        call mode_evolution(kt*dt_second, potential(:,:,1), mtor, nrad, file_unit1)
        call mode_evolution(kt*dt_second, apara(:,:,1),     mtor, nrad, file_unit2)
        call mode_evolution6(kt*dt_second, potential(:,:,1), mtor, nrad, file_unit3)
        call mode_evolution6(kt*dt_second, apara(:,:,1),     mtor, nrad, file_unit4) 
     endif

     if(myid==0 .and. mod(kt, 50)==0) call bperp_perturbation(kt*dt_second, file_unit5) 

     if(TCLR==0 .and. mod((kt-1), iplot_mode_structure)==0) then
        call mode_structure_in_xy_plane(kt, GCLR, potential(:,:,1), 'Phi')
        call mode_structure_in_xy_plane(kt, GCLR, apara(:,:,1), 'Apara')
        call mode_structure_in_xz_plane(kt, potential(:,:,1), 'Phi')
        call mode_structure_in_xz_plane(kt, apara(:,:,1), 'Apara')
        call mode_structure_in_yz_plane(kt, potential(:,:,:), 'Phi')
        call mode_structure_in_yz_plane(kt, apara(:,:,:), 'Apara')
        call mode_structure_in_poloidal_plane(kt, potential(:,:, :), 'Phi')
        call mode_structure_in_poloidal_plane(kt, apara(:,:, :), 'Apara')
     endif

  enddo !-----------time advancing loop---------

  if(myid == id_writing_evolution) then
     close(file_unit1)
     close(file_unit2)
     close(file_unit3)
     close(file_unit4)
  endif
  if(myid==0) close(file_unit5)
  if(store_restart_data .eqv. .true.) call write_data_for_restarting(kend)

1234 call cpu_time(tarray(2))  !total CPU time used by all threads of each process

  t2 = mpi_wtime()
  if (myid==0) write(*,*) 'Total CPU time of all threads of MPI process 0 (seconds): ', tarray(2)-tarray(1)
  if (myid==0) write(*,*)  'MPI_wall_time (seconds)=',  t2 - t1, 'Hours =', (t2 - t1)/3600

  call MPI_FINALIZE(ierr)
end program main


subroutine read_parameters()
  use constants, only: p_
  use control_parameters, only: kstart,kend,dt_omega_i_axis, &
       &  poloidal_angle_type,iplot_mode_structure,&
       & filter_radial, store_restart_data, &
       & fk_switch, space_charge_switch, adiabatic_electrons, &
       &   diagnosis, ismooth, nh_min, nh_max
  use magnetic_coordinates, only: nrad,mpol,mtor,pfn_inner, pfn_bdry,nsegment
  use domain_decomposition, only: ntube,myid
  use gk_module, only : nsm
  use mpi
  implicit none
  integer:: u,ierr
  namelist/control_nmlt/kstart,kend,dt_omega_i_axis, &
       & iplot_mode_structure,&
       & filter_radial, store_restart_data, &
       & poloidal_angle_type,nsegment, nrad, mpol,mtor,pfn_inner, pfn_bdry,ntube, &
       & fk_switch, space_charge_switch, adiabatic_electrons, &
       &  nh_min, nh_max, diagnosis, ismooth, nsm

  open(newunit=u,file='input.nmlt')
  read(u,control_nmlt)
  close(u)
  if(myid==0)  write(*,control_nmlt)
end subroutine read_parameters
