"""Eta-quotient q-expansions with exact integer arithmetic.
eta(z) = q^{1/24} prod_{n>=1} (1-q^n)
We work with q-series directly: prod_d (q^d; q^d)_inf^{e_d}, tracking the
q^{sum d*e_d / 24} prefactor via overall shift.
"""
from fractions import Fraction
from math import gcd

def qpoch_powers(terms, N):
    """terms: dict {d: e} for prod_d (1-q^d)^e... actually prod_d (q^d;q^d)_inf^{e_d}.
    Returns list coeff[0..N] of the product (as integers, may be negative).
    (q^d;q^d)_inf = prod_{k>=1} (1 - q^{d*k}).
    """
    # log approach with exact rationals is slow; use direct recurrence via
    # generalized pentagonal for each (q^d;q^d) then exponentiate by e_d via series power.
    from collections import defaultdict
    # base series for (q^d;q^d)_inf via pentagonal numbers
    def base_series(d, N):
        c = [0]*(N+1); c[0]=1
        k=1
        while True:
            g1 = k*(3*k-1)//2; g2 = k*(3*k+1)//2
            done=True
            for g in (g1,g2):
                idx=g*d
                if idx<=N:
                    c[idx]+= -1 if k%2==1 else 1
                    done=False
            if done and g1*d>N: break
            k+=1
            if k>10000: break
        return c
    # multiply series
    def mul(a,b,N):
        c=[0]*(N+1)
        for i in range(N+1):
            if a[i]:
                for j in range(N+1-i):
                    if b[j]: c[i+j]+=a[i]*b[j]
        return c
    def series_pow(s,e,N):
        # s[0]=1
        if e==0: return [1]+[0]*N
        if e<0:
            # invert then positive power
            inv=[0]*(N+1); inv[0]=1
            for n in range(1,N+1):
                inv[n] = -sum(s[k]*inv[n-k] for k in range(1,n+1) if s[k])
            s=inv; e=-e
        r=[1]+[0]*N
        base=s[:]
        while e:
            if e&1: r=mul(r,base,N)
            base=mul(base,base,N); e>>=1
        return r
    result=[1]+[0]*N
    for d,e in terms.items():
        if e==0: continue
        b=base_series(d,N)
        p=series_pow(b,e,N)
        result=mul(result,p,N)
    return result

def eta_quotient_qexp(edict, N):
    """edict: {d: exponent}. Returns (shift, coeffs) where the true q-expansion is
    q^{shift} * sum coeffs[n] q^n, shift = sum(d*e_d)/24 as Fraction."""
    shift = sum(Fraction(d*e,24) for d,e in edict.items())
    coeffs = qpoch_powers(edict, N)
    return shift, coeffs

if __name__ == "__main__":
    # test: eta(z)^1 -> q^{1/24}(1 - q - q^2 + q^5 + ...)
    s,c = eta_quotient_qexp({1:1}, 10)
    print("shift", s, "coeffs", c)
    # test: (q;q)_inf^3 = sum (-1)^n (2n+1) q^{n(n+1)/2}
    s,c = eta_quotient_qexp({1:3}, 10)
    print("shift", s, "coeffs", c)
