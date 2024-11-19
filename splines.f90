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
