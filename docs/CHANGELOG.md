## odepack Changelog

The intent of this log is to keep everyone in the loop about what is new
in the odepack project. It is a curated, chronologically ordered list
of notifications of notable events such as bug fixes, new features,
and usage changes.

"Do unto others as you would have them do unto you", as they say. When I
find OS (Open Source) resources, I am hoping a lot of these boxes can be
checked ...
   - [x] git repository on WWW (github)
   - [x] annotated source files with an open licence
   - [x] unit test
   - [x] make(1) build
   - [x] fpm(1) build
   - [x] user manual (on-line)
   - [x] man-page
   - [x] app program
   - [x] demo program for public procedures
   - [x] developer documents (ford(1))
   - [x] CI/CD (Continuous Integration/Development) verification (github actions)
   - [ ] registered in fpm(1) repository

---
**2024-12-14**  Modernisation Update

### :green_circle: ADD:
- Replaced 8 `dnrm2` BLAS calls with Fortran 2008 `norm2()` intrinsic
  - Files: datv.inc, dorthog.inc, dspiom.inc, dspigmr.inc
  - Reference: Hindmarsh (2001), opkd-sum.txt lines 229-237
- Added EQUALITY_AUDIT.md documenting 46 real equality tests (no changes needed)
- Added TYPE_MISMATCH_ANALYSIS.md documenting 5 excluded files and rationale
- Updated M_odepack.f90 developer notes with [DONE] status markers
- Added comprehensive test suite (test_methods.f90) verifying:
  - DLSODE MF=10 (Adams non-stiff), MF=21 (BDF user Jacobian), MF=22 (BDF internal)
  - DLSODA automatic stiff/non-stiff method switching
  - Tolerance handling (tight vs loose comparison)
  - Reference: UCRL-ID-113855, pp. 53-55, 88-91, 101-105

### :orange_circle: DIFF:
- Developer notes reorganised into "DONE" and "REMAINING WORK" sections

### :memo: NOTES:
- dcopy replacement was already complete (found commented with !X! markers)
- JROOT initialisation was already fixed in lsodkr.f90 example
- All 8 tests pass after changes (7 original + test_methods)
- Type mismatch files remain outside module (documented, intentional design)

### :book: REFERENCES:
1. Hindmarsh, A.C. (2001). "Brief Description of ODEPACK - A Systematized 
   Collection of ODE Solvers." Lawrence Livermore National Laboratory.
   20 June 2001. (opkd-sum.txt)

2. Radhakrishnan, K. and Hindmarsh, A.C. (1993). "Description and Use of 
   LSODE, the Livermore Solver for Ordinary Differential Equations."
   LLNL Report UCRL-ID-113855, December 1993. Pages 53-55 (Method Flags), 88-91 (Error Control), 101-105 (MF values).

3. NIST (1985). "FIPS PUB 69-1: Federal Information Processing Standards
   Publication - FORTRAN." National Institute of Standards and Technology.

---
**2021-02-10**  John S. Urban  <https://github.com/urbanjost>

### :green_circle: ADD:
     initial release on github
---

<!--
### :orange_circle: DIFF:
       + renamed ADVICE(3f) to ALERT(3f)
### :green_circle: ADD:
       + advice(3f) was added to provide a standardised message format simply.
### :red_circle: FIX:
       + </bo> did not work on several terminal types, changed it to a more
         universally accepted value.
-->
