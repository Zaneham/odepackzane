!===============================================================================
! Comprehensive ODEPACK Method Test Suite
!===============================================================================
! Reference: Radhakrishnan, K. and Hindmarsh, A.C. (1993). "Description and Use
!            of LSODE, the Livermore Solver for Ordinary Differential Equations."
!            LLNL Report UCRL-ID-113855, December 1993.
!            Pages 53-55 (Method Flags), 88-91 (Error Control), 101-105 (MF).
!
! Reference: Hindmarsh, A.C. (2001). "Brief Description of ODEPACK - A Systematized
!            Collection of ODE Solvers." Lawrence Livermore National Laboratory.
!
! This test suite verifies:
!   1. Method flag (MF) behaviour for different problem types
!   2. Automatic stiff/non-stiff switching (DLSODA)
!   3. Tolerance handling
!===============================================================================
program test_methods
   use M_odepack
   implicit none

   integer, parameter :: dp = kind(0.0d0)

   external :: f_harmonic, f_stiff_decay, f_simple_decay
   external :: jac_stiff_decay, jac_dummy

   logical :: all_passed
   integer :: failures

   all_passed = .true.
   failures = 0

   print '(a)', '=============================================='
   print '(a)', 'ODEPACK Comprehensive Method Test Suite'
   print '(a)', '=============================================='
   print '(a)', ''

   ! Test 1: DLSODE with MF=10 (Adams, non-stiff)
   call test_dlsode_mf10(all_passed, failures)

   ! Test 2: DLSODE with MF=21 (BDF, user-supplied Jacobian)
   call test_dlsode_mf21(all_passed, failures)

   ! Test 3: DLSODE with MF=22 (BDF, internally generated Jacobian)
   call test_dlsode_mf22(all_passed, failures)

   ! Test 4: DLSODA automatic switching - stiff problem
   call test_dlsoda_stiff(all_passed, failures)

   ! Test 5: DLSODA automatic switching - non-stiff problem
   call test_dlsoda_nonstiff(all_passed, failures)

   ! Test 6: Tolerance comparison
   call test_tolerances(all_passed, failures)

   print '(a)', ''
   print '(a)', '=============================================='
   if (all_passed) then
      print '(a)', 'ALL TESTS PASSED'
   else
      print '(a,i0,a)', 'TESTS COMPLETED WITH ', failures, ' FAILURE(S)'
   end if
   print '(a)', '=============================================='

   if (.not. all_passed) stop 1

contains

!-------------------------------------------------------------------------------
! Test DLSODE with MF=10 (Adams method, non-stiff)
! Reference: UCRL-ID-113855, pp. 101-105 - Method Flag (MF) values
!-------------------------------------------------------------------------------
subroutine test_dlsode_mf10(all_passed, failures)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures

   integer, parameter :: NEQ = 2
   integer, parameter :: LRW = 20 + 16*NEQ
   integer, parameter :: LIW = 20
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, mf
   real(kind=dp) :: y_expected, error

   print '(a)', 'Test 1: DLSODE MF=10 (Adams, non-stiff)'

   ! Simple harmonic oscillator: y'' + y = 0
   ! Written as system: y1' = y2, y2' = -y1
   ! Solution: y1 = cos(t), y2 = -sin(t) for y1(0)=1, y2(0)=0

   y(1) = 1.0d0
   y(2) = 0.0d0
   t = 0.0d0
   tout = 6.283185307179586d0  ! 2*pi

   itol = 2
   rtol = 1.0d-6
   atol(1) = 1.0d-8
   atol(2) = 1.0d-8
   itask = 1
   istate = 1
   iopt = 0
   mf = 10  ! Adams method, no Jacobian needed

   rwork = 0.0d0
   iwork = 0

   call dlsode(f_harmonic, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, mf)

   ! After 2*pi, should return to initial conditions
   y_expected = 1.0d0
   error = abs(y(1) - y_expected)

   if (istate > 0 .and. error < 1.0d-4) then
      print '(a,es10.3)', '  PASSED - Final y(1) error: ', error
   else
      print '(a,i0,a,es10.3)', '  FAILED - ISTATE=', istate, ', error=', error
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_dlsode_mf10

!-------------------------------------------------------------------------------
! Test DLSODE with MF=21 (BDF, user-supplied full Jacobian)
! Reference: UCRL-ID-113855, pp. 102-103 - MF = 21 for stiff problems
!-------------------------------------------------------------------------------
subroutine test_dlsode_mf21(all_passed, failures)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures

   integer, parameter :: NEQ = 1
   integer, parameter :: LRW = 22 + 9*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, mf
   real(kind=dp) :: y_expected, error

   print '(a)', 'Test 2: DLSODE MF=21 (BDF, user Jacobian)'

   ! Stiff decay: y' = -1000*y, y(0) = 1
   ! Solution: y = exp(-1000*t)

   y(1) = 1.0d0
   t = 0.0d0
   tout = 0.01d0

   itol = 2
   rtol = 1.0d-6
   atol(1) = 1.0d-10
   itask = 1
   istate = 1
   iopt = 0
   mf = 21  ! BDF with user-supplied Jacobian

   rwork = 0.0d0
   iwork = 0

   call dlsode(f_stiff_decay, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_stiff_decay, mf)

   y_expected = exp(-1000.0d0 * tout)
   error = abs(y(1) - y_expected) / y_expected

   if (istate > 0 .and. error < 1.0d-4) then
      print '(a,es10.3)', '  PASSED - Relative error: ', error
   else
      print '(a,i0,a,es10.3)', '  FAILED - ISTATE=', istate, ', rel error=', error
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_dlsode_mf21

!-------------------------------------------------------------------------------
! Test DLSODE with MF=22 (BDF, internally generated Jacobian)
! Reference: UCRL-ID-113855, pp. 102-103 - MF = 22 (internal difference quotient)
!-------------------------------------------------------------------------------
subroutine test_dlsode_mf22(all_passed, failures)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures

   integer, parameter :: NEQ = 1
   integer, parameter :: LRW = 22 + 9*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, mf
   real(kind=dp) :: y_expected, error

   print '(a)', 'Test 3: DLSODE MF=22 (BDF, internal Jacobian)'

   ! Same stiff decay problem, but with internal Jacobian

   y(1) = 1.0d0
   t = 0.0d0
   tout = 0.01d0

   itol = 2
   rtol = 1.0d-6
   atol(1) = 1.0d-10
   itask = 1
   istate = 1
   iopt = 0
   mf = 22  ! BDF with internal Jacobian

   rwork = 0.0d0
   iwork = 0

   call dlsode(f_stiff_decay, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, mf)

   y_expected = exp(-1000.0d0 * tout)
   error = abs(y(1) - y_expected) / y_expected

   if (istate > 0 .and. error < 1.0d-4) then
      print '(a,es10.3)', '  PASSED - Relative error: ', error
   else
      print '(a,i0,a,es10.3)', '  FAILED - ISTATE=', istate, ', rel error=', error
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_dlsode_mf22

!-------------------------------------------------------------------------------
! Test DLSODA automatic method switching - stiff problem
! Reference: UCRL-ID-113855, pp. 53, 119 - LSODA automatic switching
!-------------------------------------------------------------------------------
subroutine test_dlsoda_stiff(all_passed, failures)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures

   integer, parameter :: NEQ = 1
   integer, parameter :: LRW = 22 + 16*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, jt
   real(kind=dp) :: y_expected, error

   print '(a)', 'Test 4: DLSODA auto-switch (stiff problem)'

   ! Stiff decay - DLSODA should switch to BDF

   y(1) = 1.0d0
   t = 0.0d0
   tout = 0.01d0

   itol = 2
   rtol = 1.0d-6
   atol(1) = 1.0d-10
   itask = 1
   istate = 1
   iopt = 0
   jt = 2  ! Internal Jacobian

   rwork = 0.0d0
   iwork = 0

   call dlsoda(f_stiff_decay, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, jt)

   y_expected = exp(-1000.0d0 * tout)
   error = abs(y(1) - y_expected) / y_expected

   if (istate > 0 .and. error < 1.0d-4) then
      print '(a,es10.3)', '  PASSED - Relative error: ', error
   else
      print '(a,i0,a,es10.3)', '  FAILED - ISTATE=', istate, ', rel error=', error
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_dlsoda_stiff

!-------------------------------------------------------------------------------
! Test DLSODA automatic method switching - non-stiff problem
! Reference: UCRL-ID-113855, pp. 53-54 - LSODA uses Adams for non-stiff
!-------------------------------------------------------------------------------
subroutine test_dlsoda_nonstiff(all_passed, failures)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures

   integer, parameter :: NEQ = 2
   integer, parameter :: LRW = 22 + 16*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, jt
   real(kind=dp) :: y_expected, error

   print '(a)', 'Test 5: DLSODA auto-switch (non-stiff problem)'

   ! Simple harmonic oscillator - should use Adams method

   y(1) = 1.0d0
   y(2) = 0.0d0
   t = 0.0d0
   tout = 6.283185307179586d0  ! 2*pi

   itol = 2
   rtol = 1.0d-6
   atol(1) = 1.0d-8
   atol(2) = 1.0d-8
   itask = 1
   istate = 1
   iopt = 0
   jt = 2

   rwork = 0.0d0
   iwork = 0

   call dlsoda(f_harmonic, [NEQ], y, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, jt)

   y_expected = 1.0d0
   error = abs(y(1) - y_expected)

   if (istate > 0 .and. error < 1.0d-4) then
      print '(a,es10.3)', '  PASSED - Final y(1) error: ', error
   else
      print '(a,i0,a,es10.3)', '  FAILED - ISTATE=', istate, ', error=', error
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_dlsoda_nonstiff

!-------------------------------------------------------------------------------
! Test tolerance handling
! Reference: UCRL-ID-113855, pp. 88-91 - Error control parameters
!-------------------------------------------------------------------------------
subroutine test_tolerances(all_passed, failures)
   logical, intent(inout) :: all_passed
   integer, intent(inout) :: failures

   integer, parameter :: NEQ = 1
   integer, parameter :: LRW = 22 + 9*NEQ + NEQ**2
   integer, parameter :: LIW = 20 + NEQ
   real(kind=dp) :: y_tight(NEQ), y_loose(NEQ), rwork(LRW), atol(NEQ)
   integer :: iwork(LIW)
   real(kind=dp) :: t, tout, rtol
   integer :: itol, itask, istate, iopt, mf
   real(kind=dp) :: y_exact, err_tight, err_loose

   print '(a)', 'Test 6: Tolerance comparison (tight vs loose)'

   ! Solve y' = -y with y(0) = 1, exact solution y = exp(-t)
   tout = 1.0d0
   y_exact = exp(-tout)

   ! Tight tolerance
   y_tight(1) = 1.0d0
   t = 0.0d0
   itol = 2
   rtol = 1.0d-10
   atol(1) = 1.0d-12
   itask = 1
   istate = 1
   iopt = 0
   mf = 22
   rwork = 0.0d0
   iwork = 0

   call dlsode(f_simple_decay, [NEQ], y_tight, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, mf)
   err_tight = abs(y_tight(1) - y_exact)

   ! Loose tolerance
   y_loose(1) = 1.0d0
   t = 0.0d0
   rtol = 1.0d-4
   atol(1) = 1.0d-6
   istate = 1
   rwork = 0.0d0
   iwork = 0

   call dlsode(f_simple_decay, [NEQ], y_loose, t, tout, itol, [rtol], atol, itask, &
               istate, iopt, rwork, LRW, iwork, LIW, jac_dummy, mf)
   err_loose = abs(y_loose(1) - y_exact)

   ! Tight tolerance should give smaller error
   if (err_tight < err_loose) then
      print '(a,es10.3,a,es10.3)', '  PASSED - Tight err: ', err_tight, ', Loose err: ', err_loose
   else
      print '(a,es10.3,a,es10.3)', '  FAILED - Tight err: ', err_tight, ', Loose err: ', err_loose
      all_passed = .false.
      failures = failures + 1
   end if

end subroutine test_tolerances

end program test_methods

!===============================================================================
! External test problem right-hand sides
!===============================================================================

! Simple harmonic oscillator: y1' = y2, y2' = -y1
subroutine f_harmonic(neq, t, y, ydot)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: ydot(neq)
   ydot(1) = y(2)
   ydot(2) = -y(1)
end subroutine f_harmonic

! Stiff decay: y' = -1000*y
subroutine f_stiff_decay(neq, t, y, ydot)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: ydot(neq)
   ydot(1) = -1000.0d0 * y(1)
end subroutine f_stiff_decay

! Jacobian for stiff decay
subroutine jac_stiff_decay(neq, t, y, ml, mu, pd, nrowpd)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq, ml, mu, nrowpd
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: pd(nrowpd, neq)
   pd(1,1) = -1000.0d0
end subroutine jac_stiff_decay

! Simple decay: y' = -y
subroutine f_simple_decay(neq, t, y, ydot)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: ydot(neq)
   ydot(1) = -y(1)
end subroutine f_simple_decay

! Dummy Jacobian (not used when MF has internal Jacobian)
subroutine jac_dummy(neq, t, y, ml, mu, pd, nrowpd)
   implicit none
   integer, parameter :: dp = kind(0.0d0)
   integer, intent(in) :: neq, ml, mu, nrowpd
   real(kind=dp), intent(in) :: t
   real(kind=dp), intent(in) :: y(neq)
   real(kind=dp), intent(out) :: pd(nrowpd, neq)
   ! Not used
end subroutine jac_dummy
