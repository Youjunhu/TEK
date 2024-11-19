module rnorm_mod !generate Gaussian random numbers
implicit none
private
public :: dp, rnorm_polar, rnorm_ratio_uni, rnorm_polar_vec, rnorm_ratio_uni_vec, &
          rnorm_box_muller, rnorm_box_muller_vec
integer, parameter :: dp = kind(1.0d0)
real(kind=dp), parameter :: half = 0.5_dp, two_pi = 6.28318530718_dp
contains

function rnorm_polar_vec(n,mu,sigma) result(variates)
integer      , intent(in)           :: n     ! # of normal variates
real(kind=dp), intent(in), optional :: mu    ! target mean 
real(kind=dp), intent(in), optional :: sigma ! target standard deviation
real(kind=dp)                       :: variates(n) ! normal variates
integer                             :: i
do i=1,n
   variates(i) = rnorm_polar()
end do
if (present(sigma)) variates = sigma*variates
if (present(mu)) variates = variates + mu
end function rnorm_polar_vec
!
function rnorm_ratio_uni_vec(n,mu,sigma) result(variates)
integer      , intent(in)           :: n     ! # of normal variates
real(kind=dp), intent(in), optional :: mu    ! target mean 
real(kind=dp), intent(in), optional :: sigma ! target standard deviation
real(kind=dp)                       :: variates(n) ! normal variates
integer                             :: i
do i=1,n
   variates(i) = rnorm_ratio_uni()
end do
if (present(sigma)) variates = sigma*variates
if (present(mu)) variates = variates + mu
end function rnorm_ratio_uni_vec
!
function rnorm_box_muller_vec(n,mu,sigma) result(variates)
integer      , intent(in)           :: n     ! # of normal variates
real(kind=dp), intent(in), optional :: mu    ! target mean 
real(kind=dp), intent(in), optional :: sigma ! target standard deviation
real(kind=dp)                       :: variates(n) ! normal variates
integer                             :: i,j
logical                             :: n_odd
n_odd = mod(n,2) /= 0
do i=1,n/2
   j = 2*i - 1
   variates(j:j+1) = rnorm_box_muller()
end do
if (n_odd) variates(n) = rnorm_box_muller_single_variate()
if (present(sigma)) variates = sigma*variates
if (present(mu)) variates = variates + mu
end function rnorm_box_muller_vec
!
FUNCTION rnorm_ratio_uni() RESULT(fn_val) ! adapted from https://jblevins.org/mirror/amiller/ran_norm.f90
! Adapted from the following Fortran 77 code
!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

!  The function random_normal() returns a normally distributed pseudo-random
!  number with zero mean and unit variance.

!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
!  and J.F. Monahan augmented with quadratic bounding curves.

real(kind=dp) :: fn_val
real(kind=dp), parameter :: s = 0.449871_dp, t = -0.386595_dp, a = 0.19600_dp, b = 0.25472_dp,    &
            r1 = 0.27597_dp, r2 = 0.27846_dp
! local variables
real(kind=dp) :: u, v, x, y, q
!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
DO
  CALL RANDOM_NUMBER(u)
  CALL RANDOM_NUMBER(v)
  v = 1.7156_dp * (v - half)
!     Evaluate the quadratic form
  x = u - s
  y = ABS(v) - t
  q = x**2 + y*(a*y - b*x)
!     Accept P if inside inner ellipse
  IF (q < r1) EXIT
!     Reject P if outside outer ellipse
  IF (q > r2) CYCLE
!     Reject P if outside acceptance region
  IF (v**2 < -4.0_dp*LOG(u)*u**2) EXIT
END DO
!     Return ratio of P's coordinates as the normal deviate
fn_val = v/u
END FUNCTION rnorm_ratio_uni
!
function rnorm_box_muller() result(variates) ! coded formulas from https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
! return two uncorrelated standard normal variates
real(kind=dp) :: variates(2)
real(kind=dp) :: u(2), factor, arg
do
   call random_number(u)
   if (u(1) > 0.0_dp) exit
end do
factor = sqrt(-2 * log(u(1)))
arg = two_pi*u(2)
variates = factor * [cos(arg),sin(arg)]
end function rnorm_box_muller
!
function rnorm_box_muller_single_variate() result(variate)
! return a standard normal variate
real(kind=dp) :: variate
real(kind=dp) :: u(2), factor, arg
call random_number(u)
factor = sqrt(-2 * log(u(1)))
arg = two_pi*u(2)
variate = factor * cos(arg)
end function rnorm_box_muller_single_variate
!
FUNCTION rnorm_polar() RESULT(fn_val) ! adapted from https://jblevins.org/mirror/amiller/rnorm.f90

!   Generate a random normal deviate using the polar method.
!   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
!              normal variables', Siam Rev., vol.6, 260-264, 1964.

IMPLICIT NONE
REAL(kind=dp)  :: fn_val

! Local variables

REAL(kind=dp)   :: u, sum
REAL(kind=dp), SAVE      :: v, sln
LOGICAL, SAVE   :: second = .FALSE.
REAL(kind=dp), PARAMETER :: one = 1.0_dp, vsmall = TINY( one )

IF (second) THEN
! If second, use the second random number generated on last call

  second = .false.
  fn_val = v*sln

ELSE
! First call; generate a pair of random normals

  second = .true.
  DO
    CALL RANDOM_NUMBER( u )
    CALL RANDOM_NUMBER( v )
    u = SCALE( u, 1 ) - one
    v = SCALE( v, 1 ) - one
    sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
    IF(sum < one) EXIT
  END DO
  sln = SQRT(-SCALE(LOG(sum),1)/sum)
  fn_val = u*sln
END IF
END FUNCTION rnorm_polar
end module rnorm_mod
!!$!
!!$program main
!!$use rnorm_mod, only: dp, rnorm_vec, rnorm_box_muller_vec
!!$implicit none
!!$integer          , parameter :: n = 10**8, niter = 5, nmethods = 3
!!$real(kind=dp)    , parameter :: mu = 2.0_dp, sigma = 3.0_dp
!!$character (len=*), parameter :: methods(nmethods) = ["polar     ", &
!!$                                                     "ratio_uni ","box_muller"]
!!$integer                      :: ipow, iter, imethod
!!$real(kind=dp)                :: xmean, old_time, new_time, elapsed_time(nmethods)
!!$real(kind=dp), allocatable   :: x(:)
!!$call random_seed()
!!$elapsed_time = 0.0_dp
!!$allocate (x(n))
!!$print "(*(a10))",    "n", "niter", "mu", "sigma"
!!$print "(2i10,2f10.4)", n ,  niter,   mu  , sigma
!!$do imethod=1,nmethods
!!$   print "(/,'method = ',a)",trim(methods(imethod))
!!$   print "('central moments',/,a10,*(i10))","mean",(ipow,ipow=1,4)
!!$   do iter=1,niter
!!$      call cpu_time(old_time)
!!$      x = rnorm_vec(n,methods(imethod),mu,sigma)
!!$      call cpu_time(new_time)
!!$      elapsed_time(imethod) = elapsed_time(imethod) + new_time - old_time
!!$      xmean = sum(x)/n
!!$      x = x - xmean
!!$      print "(*(f10.4))",xmean,(sum(x**ipow)/n,ipow=1,4)
!!$   end do
!!$   print*,"theoretical"
!!$   print "(*(f10.4))",mu,0.0_dp,sigma**2,0.0_dp,3*sigma**4
!!$end do
!!$print "(/,*(a15))", "method",(trim(methods(imethod)),imethod=1,nmethods)
!!$print "(a15,*(f15.4))", "cpu_time",elapsed_time
!!$print "(/,'check Box-Muller for n odd or even')"
!!$print "(*(f10.4))",rnorm_box_muller_vec(3),rnorm_box_muller_vec(4)
!!$end program main
