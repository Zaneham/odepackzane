# Real Equality Tests Audit

**Reference**: M_odepack.f90 developer note - "many equality tests for real values; could replace"

## Analysis Summary

After reviewing 46 occurrences of `==0.0D0` comparisons across 25 files, these fall into categories:

### 1. User Parameter Checks (SAFE - Intentional)
- `if (h0==0.0D0)` - Checking if user specified initial step size = 0 (let solver choose)
- `if (dlpk%delt==0.0D0)` - Checking default preconditioner parameter

**Reason to keep**: User explicitly sets 0.0 to request default behaviour.

### 2. Matrix Singularity Detection (SAFE - Intentional)
- `if (A(l,nm1)==0.0D0)` - Pivot checks in dhefa.inc, dheqr.inc
- `if (A(N,N)==0.0D0)` - Diagonal singularity checks
- `if (t2==0.0D0)` - Householder transformation checks

**Reason to keep**: Exact zero indicates singular matrix. Epsilon comparison would miss actual singularities.

### 3. Convergence and Early Exit Conditions (SAFE - Intentional)
- `if (sumdsq==0.0D0)` - Gram-Schmidt orthogonalisation in dorthog.inc
- `if (ztr0==0.0D0)` - PCG iteration check in dpcg.inc, dpcgs.inc
- `if (ptw==0.0D0)` - Preconditioner check

**Reason to keep**: Mathematical conditions where exact zero has meaning.

### 4. Sparse Matrix Entry Checks (SAFE - Intentional)
- `if (Wk(ljfo)==0.0D0)` - Checking for structural zeros
- `if (Savr(i)==0.0D0)` - Checking for zero entries

**Reason to keep**: Sparsity detection requires exact zero check.

## Recommendation

**No changes recommended.** All equality tests serve valid purposes:
- User-specified defaults
- Matrix singularity detection
- Mathematical convergence conditions
- Sparsity pattern detection

Changing these to epsilon comparisons could introduce subtle bugs.

## Reference

Hindmarsh, A.C. (2001). "Brief Description of ODEPACK - A Systematized Collection of ODE Solvers."
Lawrence Livermore National Laboratory, UCRL-CODE-113855. 20 June 2001.
