!===============================================================================
! ODEPACK Numerical Regression Test Suite
!===============================================================================
!
! PURPOSE:
!   Detect numerical drift across compilers, platforms, and library versions
!   by comparing computed results against validated golden reference values.
!
! RATIONALE:
!   IEEE 754 compliance guarantees representation, not computation order.
!   Different compilers, optimisation levels, BLAS implementations, and
!   hardware (FMA, SSE, AVX, GPU) can produce results that differ at the
!   limits of floating-point precision. These small differences can
!   accumulate in iterative solvers, causing "number drift".
!
! REFERENCES:
!   1. Goldberg, D. (1991). "What every computer scientist should know about
!      floating-point arithmetic." ACM Computing Surveys, 23(1), 5-48.
!      https://doi.org/10.1145/103162.103163
!
!   2. Oberkampf, W.L. & Roy, C.J. (2010). "Verification and Validation in
!      Scientific Computing." Cambridge University Press.
!      ISBN: 978-0521113601
!
!   3. IEEE 754-2019. "IEEE Standard for Floating-Point Arithmetic."
!      https://standards.ieee.org/ieee/754/6210/
!
!   4. Salari, K. & Knupp, P. (2000). "Code verification by the method of
!      manufactured solutions." SAND2000-1444, Sandia National Laboratories.
!      https://www.osti.gov/biblio/759450
!
!   5. Radhakrishnan, K. & Hindmarsh, A.C. (1993). "Description and Use of
!      LSODE." LLNL Report UCRL-ID-113855, pp. 88-91 (Error Control).
!
! GOLDEN REFERENCE METADATA:
!   Generated: 2024-12-14
!   Compiler:  gfortran 13.2.0
!   Flags:     -O2 -fallow-argument-mismatch
!   Platform:  x86_64-w64-mingw32 (Windows)
!   BLAS:      Internal (ODEPACK bundled)
!   Precision: IEEE 754 binary64 (double precision)
!
! TOLERANCE STRATEGY:
!   - Tight (1e-12): Detect compiler/platform changes
!   - Moderate (1e-8): Detect algorithm changes
!   - Loose (1e-4): Detect gross errors only
!
!===============================================================================
program test_regression
   use M_odepack
   implicit none

   integer, parameter :: dp = kind(0.0d0)

   external :: f_decay, f_robertson, jac_robertson, jac_dummy

   logical :: all_passed
   integer :: failures, warnings

   all_passed = .true.
   failures = 0
   warnings = 0

   print '(a)', '==============================================================='
   print '(a)', 'ODEPACK Numerical Regression Test Suite'
   print '(a)', 'Reference: Goldberg (1991), Oberkampf & Roy (2010)'
   print '(a)', '==============================================================='
   print '(a)', ''

   ! Test 1: Simple exponential decay - analytical solution known exactly
   call test_exponential_decay(all_passed, failures, warnings)

   ! Test 2: Robertson chemical kinetics - classic stiff test problem
   call test_robertson(all_passed, failures, warnings)

   ! Test 3: Harmonic oscillator - periodic solution, sensitive to drift
   call test_harmonic_oscillator(all_passed, failures, warnings)

   print '(a)', ''
   print '(a)', '==============================================================='
   if (all_passed) then
      print '(a)', 'ALL REGRESSION TESTS PASSED'
   else
      print '(a,i0,a,i0,a)', 'COMPLETED: ', failures, ' failure(s), ', warnings, ' warning(s)'
   end if
   print '(a)', '==============================================================='

   if (failures > 0) stop 1

contains

!-------------------------------------------------------------------------------
! Test 1: Exponential Decay
! Problem: y' = -y, y(0) = 1
! Exact solution: y(t) = exp(-t)
! Reference: Oberkampf & Roy (2010), Chapter 6 - Exact Solutions
!-------------------------------------------------------------------------------
subroutine test_exponential_decay(all_passed, failures, warnings)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures, warnings

   integer, parameter :: NEQ = 1
   integer, parameter :: LRW = 22 + 9*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, mf

   ! Golden reference values (analytical: y = exp(-t))
   real(kind=dp), parameter :: t_test = 1.0d0
   real(kind=dp), parameter :: y_exact = 0.36787944117144232d0  ! exp(-1)
   real(kind=dp), parameter :: tol_tight = 1.0d-12
   real(kind=dp), parameter :: tol_moderate = 1.0d-8

   real(kind=dp) :: error

   print '(a)', 'Test 1: Exponential Decay (y'' = -y)'
   print '(a)', '  Reference: Analytical solution y = exp(-t)'

   y(1) = 1.0d0
   t = 0.0d0
   tout = t_test

   itol = 2
   rtol = 1.0d-12
   atol(1) = 1.0d-14
   itask = 1
   istate = 1
   iopt = 0
   mf = 22

   rwork = 0.0d0
   iwork = 0

   call dlsode(f_decay, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, mf)

   error = abs(y(1) - y_exact)

   if (error < tol_tight) then
      print '(a,es12.5)', '  PASSED (tight)  - Error: ', error
   else if (error < tol_moderate) then
      print '(a,es12.5)', '  WARNING (drift) - Error: ', error
      print '(a)',        '    Possible compiler/platform difference'
      warnings = warnings + 1
   else
      print '(a,es12.5)', '  FAILED          - Error: ', error
      print '(a,es22.15)','    Expected: ', y_exact
      print '(a,es22.15)','    Got:      ', y(1)
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_exponential_decay

!-------------------------------------------------------------------------------
! Test 2: Robertson Chemical Kinetics
! Classic stiff test problem from ODEPACK documentation
! Reference: UCRL-ID-113855, pp. 3-4; Hindmarsh (1983)
!
! y1' = -0.04*y1 + 1e4*y2*y3
! y2' =  0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
! y3' =  3e7*y2^2
!
! Initial: y1=1, y2=0, y3=0
! This is the canonical ODEPACK test problem.
!-------------------------------------------------------------------------------
subroutine test_robertson(all_passed, failures, warnings)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures, warnings

   integer, parameter :: NEQ = 3
   integer, parameter :: LRW = 22 + 9*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, mf

   ! Golden reference at t=0.4 (from validated ODEPACK runs)
   ! These values are well-established in the literature
   real(kind=dp), parameter :: y1_golden = 0.98517d0
   real(kind=dp), parameter :: y2_golden = 3.3864d-5
   real(kind=dp), parameter :: y3_golden = 1.4794d-2
   real(kind=dp), parameter :: tol_moderate = 1.0d-4  ! Stiff problems have more variation

   real(kind=dp) :: err1, err2, err3, max_err

   print '(a)', 'Test 2: Robertson Chemical Kinetics (stiff)'
   print '(a)', '  Reference: UCRL-ID-113855, pp. 3-4'

   y(1) = 1.0d0
   y(2) = 0.0d0
   y(3) = 0.0d0
   t = 0.0d0
   tout = 0.4d0

   itol = 2
   rtol = 1.0d-4
   atol(1) = 1.0d-6
   atol(2) = 1.0d-10
   atol(3) = 1.0d-6
   itask = 1
   istate = 1
   iopt = 0
   mf = 21

   rwork = 0.0d0
   iwork = 0

   call dlsode(f_robertson, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_robertson, mf)

   err1 = abs(y(1) - y1_golden) / y1_golden
   err2 = abs(y(2) - y2_golden) / y2_golden
   err3 = abs(y(3) - y3_golden) / y3_golden
   max_err = max(err1, err2, err3)

   if (istate < 0) then
      print '(a,i0)', '  FAILED - Solver error, ISTATE = ', istate
      all_passed = .false.
      failures = failures + 1
   else if (max_err < tol_moderate) then
      print '(a,es12.5)', '  PASSED - Max relative error: ', max_err
   else
      print '(a,es12.5)', '  FAILED - Max relative error: ', max_err
      print '(a,3es14.6)','    Expected: ', y1_golden, y2_golden, y3_golden
      print '(a,3es14.6)','    Got:      ', y(1), y(2), y(3)
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_robertson

!-------------------------------------------------------------------------------
! Test 3: Harmonic Oscillator
! Problem: y'' + y = 0  =>  y1' = y2, y2' = -y1
! Exact solution: y1 = cos(t), y2 = -sin(t) for y1(0)=1, y2(0)=0
!
! After t = 2*pi, should return to initial conditions.
! Periodic problems are sensitive to accumulated drift.
! Reference: Goldberg (1991), Section on error accumulation
!-------------------------------------------------------------------------------
subroutine test_harmonic_oscillator(all_passed, failures, warnings)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures, warnings

   integer, parameter :: NEQ = 2
   integer, parameter :: LRW = 22 + 16*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, jt

   ! After one full period, y should return to initial state
   real(kind=dp), parameter :: two_pi = 6.283185307179586d0
   real(kind=dp), parameter :: y1_golden = 1.0d0
   real(kind=dp), parameter :: y2_golden = 0.0d0
   real(kind=dp), parameter :: tol_tight = 1.0d-6
   real(kind=dp), parameter :: tol_moderate = 1.0d-4

   real(kind=dp) :: error

   print '(a)', 'Test 3: Harmonic Oscillator (periodic)'
   print '(a)', '  Reference: Goldberg (1991) - error accumulation'

   y(1) = 1.0d0
   y(2) = 0.0d0
   t = 0.0d0
   tout = two_pi

   itol = 2
   rtol = 1.0d-8
   atol(1) = 1.0d-10
   atol(2) = 1.0d-10
   itask = 1
   istate = 1
   iopt = 0
   jt = 2

   rwork = 0.0d0
   iwork = 0

   call dlsoda(f_harmonic, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, jt)

   error = sqrt((y(1) - y1_golden)**2 + (y(2) - y2_golden)**2)

   if (error < tol_tight) then
      print '(a,es12.5)', '  PASSED (tight)  - Error: ', error
   else if (error < tol_moderate) then
      print '(a,es12.5)', '  WARNING (drift) - Error: ', error
      warnings = warnings + 1
   else
      print '(a,es12.5)', '  FAILED          - Error: ', error
      print '(a,2es14.6)','    Expected: ', y1_golden, y2_golden
      print '(a,2es14.6)','    Got:      ', y(1), y(2)
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_harmonic_oscillator

!-------------------------------------------------------------------------------
! Internal: Harmonic oscillator RHS
!-------------------------------------------------------------------------------
subroutine f_harmonic(neq, t, y, ydot)
   integer, intent(in) :: neq
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: ydot(neq)
   ydot(1) = y(2)
   ydot(2) = -y(1)
end subroutine f_harmonic

end program test_regression

!===============================================================================
! External subroutines (required by ODEPACK interface)
!===============================================================================

! Simple decay: y' = -y
subroutine f_decay(neq, t, y, ydot)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: ydot(neq)
   ydot(1) = -y(1)
end subroutine f_decay

! Robertson chemical kinetics
subroutine f_robertson(neq, t, y, ydot)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: ydot(neq)
   ydot(1) = -0.04d0*y(1) + 1.0d4*y(2)*y(3)
   ydot(3) = 3.0d7*y(2)*y(2)
   ydot(2) = -ydot(1) - ydot(3)
end subroutine f_robertson

! Jacobian for Robertson
subroutine jac_robertson(neq, t, y, ml, mu, pd, nrowpd)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq, ml, mu, nrowpd
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: pd(nrowpd, neq)
   pd(1,1) = -0.04d0
   pd(1,2) = 1.0d4*y(3)
   pd(1,3) = 1.0d4*y(2)
   pd(2,1) = 0.04d0
   pd(2,2) = -1.0d4*y(3) - 6.0d7*y(2)
   pd(2,3) = -1.0d4*y(2)
   pd(3,1) = 0.0d0
   pd(3,2) = 6.0d7*y(2)
   pd(3,3) = 0.0d0
end subroutine jac_robertson

! Dummy Jacobian
subroutine jac_dummy(neq, t, y, ml, mu, pd, nrowpd)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq, ml, mu, nrowpd
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: pd(nrowpd, neq)
   ! Not used
end subroutine jac_dummy
