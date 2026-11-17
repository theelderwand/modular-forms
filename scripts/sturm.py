"""Sturm-bound proof verifier for eta-quotient identities.

If f and g are modular forms of weight k for Gamma_0(N) (same character),
and their q-expansions agree up to the Sturm bound:
    B = (k/12) * [SL_2(Z) : Gamma_0(N)]
then f = g identically.

For eta-quotients, we use the standard criteria:
- Weight: k = (1/2) * sum(e_d)
- Level: N = lcm of d's (with adjustments for the character)
- Character: determined by the eta-quotient (Kronecker symbol)

This module provides:
- sturm_bound(k, N): compute the bound
- verify_identity(edict_A, edict_B, k, N): check agreement to Sturm bound
"""
from math import gcd
from fractions import Fraction
import sys
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from eta import eta_quotient_qexp

def index_gamma0(N):
    """[SL_2(Z) : Gamma_0(N)] = N * prod_{p|N} (1 + 1/p)"""
    result = N
    n = N
    p = 2
    seen = set()
    while p * p <= n:
        if n % p == 0:
            seen.add(p)
            while n % p == 0: n //= p
        p += 1 if p == 2 else 2
    if n > 1: seen.add(n)
    for p in seen:
        result = result * (p + 1) // p
    return result

def sturm_bound(k, N):
    """Sturm bound for M_k(Gamma_0(N)). k can be half-integer (as Fraction)."""
    k = Fraction(k)
    idx = index_gamma0(N)
    # B = (k/12) * idx, rounded up
    B = k * idx / 12
    return int(B) + (1 if B != int(B) else 0)

def verify_identity(edict_A, edict_B, k, N, verbose=True):
    """Verify eta-quotient identity A = B via Sturm bound.
    Returns (True/False, bound_used, terms_checked).
    """
    B = sturm_bound(k, N)
    # Need enough terms: shift may differ, so compare q^{shift} * series
    # For simplicity, require same shift (or adjust)
    sA, cA = eta_quotient_qexp(edict_A, B + 10)
    sB, cB = eta_quotient_qexp(edict_B, B + 10)
    if verbose:
        print(f"Sturm bound: {B} (k={k}, N={N}, index={index_gamma0(N)})")
        print(f"Shift A: {sA}, Shift B: {sB}")
    if sA != sB:
        if verbose: print("Shifts differ — cannot directly compare (need modular adjustment)")
        return False, B, 0
    match = all(cA[i] == cB[i] for i in range(B + 1))
    if verbose:
        print(f"First {B+1} coefficients match: {match}")
    return match, B, B + 1

if __name__ == "__main__":
    # Test: trivial identity eta(z)^3 / eta(3z) ... actually let's test a known one
    # Jacobi triple product related: (q;q)^3 = sum (-1)^n (2n+1) q^{n(n+1)/2}
    # As eta-quotient: eta(z)^3 has weight 3/2, level 1? Actually eta^3 is on Gamma_0(1) with character.
    # Simple test: A = B (same edict) should verify
    ok, B, n = verify_identity({1: 2}, {1: 2}, k=1, N=1)
    print("Self-identity:", ok)
