# Type Mismatch Analysis

**Reference**: Hindmarsh, A.C. (2001). opkd-sum.txt, lines 248-264, 272-285.

## Problem Statement

Five files are excluded from the M_odepack module interface due to type mismatches:

- `src/M_da1/dprep.inc`
- `src/M_da1/dainvgs.inc`
- `src/M_da1/dprepi.inc`
- `src/M_da1/dstodi.inc`
- `src/M_da1/dstode.inc`

## Root Cause

From opkd-sum.txt (Hindmarsh, 2001):
> "In various places in the LSODES and LSODIS solvers, a call to a subroutine
> has a subscripted real array as an argument where the subroutine called
> expects an integer array... This is done in order to use work space in an
> efficient manner, as the same space is sometimes used for real work space
> and sometimes for integer work space."

### Technical Details

1. **RWORK array** is declared as `REAL(8)` (double precision)
2. **Subroutines** like DPREP expect both:
   - `Wk(*)` - real work array
   - `Iwk(*)` - integer work array
3. **Callers** pass `RWORK(offset)` to BOTH parameters

Example from dprep.inc documentation:
```
!! WK : a real work array of length LENWK
!! IWK: integer work array, assumed to occupy the same space as WK.
```

### LENRAT Constant

The code uses LENRAT (real-to-integer wordlength ratio):
- Usually 2 for double precision (8-byte real, 4-byte integer)
- Controls how integer space is carved from real array
- Set in DLSODES, DLSODIS, and CDRV subroutines

## Possible Solutions

### Option 1: Separate Arrays (Breaking Change)
- Split RWORK into RWORK + IWORK_EXTRA
- Requires API changes to all callers
- Most correct but incompatible with existing code

### Option 2: TRANSFER() Intrinsic (Fortran 2003)
- Use TRANSFER() to reinterpret memory
- Maintains same memory layout
- Still non-standard but more explicit

### Option 3: C Interoperability (Fortran 2003)
- Use ISO_C_BINDING with c_ptr
- Type-safe memory reinterpretation
- More complex implementation

### Option 4: Keep Current Approach
- Use `-fallow-argument-mismatch` flag (gfortran)
- Document the non-standard behaviour
- Maintains backward compatibility

## Recommendation

**Option 4 (Keep Current)** for this PR:
- The type punning is intentional and documented
- Changing it would be a major refactoring effort
- Current approach works with compiler flags
- Backward compatibility is maintained

For future work, Option 2 or 3 could be explored in a dedicated PR
after comprehensive testing with real ODEPACK users.

## Files Currently Outside Module

These are included via `include` at module end, after `end module`:
```fortran
!- TYPE MISMATCH
include "M_da1/dprep.inc"
include "M_da1/dainvgs.inc"
include "M_da1/dprepi.inc"
include "M_da1/dstodi.inc"
include "M_da1/dstode.inc"
```

This allows them to compile without interface checking.

## References

1. Hindmarsh, A.C. (2001). "Brief Description of ODEPACK - A Systematized 
   Collection of ODE Solvers." Lawrence Livermore National Laboratory.
   20 June 2001.

2. Radhakrishnan, K. and Hindmarsh, A.C. (1993). "Description and Use of 
   LSODE, the Livermore Solver for Ordinary Differential Equations."
   LLNL Report UCRL-ID-113855, December 1993.
