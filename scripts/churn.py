"""Churn search engine: eta-quotient enumeration (Newman + Ligozat),
character via squarefree kernel, exact nullspace relations, theta-series tools.
"""
import sys, itertools
from fractions import Fraction
from math import gcd, isqrt
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from eta import eta_quotient_qexp
import sympy as sp

def divisors(n):
    ds = []
    for d in range(1, n + 1):
        if n % d == 0:
            ds.append(d)
    return ds

def newman_ok(edict, N):
    """edict {d: r_d}; both congruences mod 24."""
    s1 = sum(r * d for d, r in edict.items())
    s2 = sum(r * (N // d) for d, r in edict.items())
    return s1 % 24 == 0 and s2 % 24 == 0

def ligozat_orders(edict, N):
    """Order of vanishing at each cusp 1/c, c|N. Returns dict or None if negative.
    Integer arithmetic: ord = (sum_d r_d*gcd(d,c)^2*(N/d)) / (24*gcd(c^2,N))."""
    orders = {}
    for c in divisors(N):
        g = gcd(c * c, N)
        num = 0
        for d, r in edict.items():
            num += r * gcd(d, c) ** 2 * (N // d)
        if num < 0:
            return None
        orders[c] = Fraction(num, 24 * g)
    return orders

def _sqfree_int(n):
    """Squarefree kernel of a nonzero integer (sign preserved)."""
    if n == 0:
        return 0
    sign = 1 if n > 0 else -1
    n = abs(n)
    m = 1
    p = 2
    while p * p <= n:
        cnt = 0
        while n % p == 0:
            n //= p
            cnt += 1
        if cnt % 2 == 1:
            m *= p
        p += 1 if p == 2 else 2
    if n > 1:
        m *= n
    return sign * m

def squarefree_kernel(n):
    """Squarefree integer m with n = (rational square) * m. Accepts int/Fraction."""
    n = Fraction(n)
    if n == 0:
        return 0
    return _sqfree_int(n.numerator) * _sqfree_int(n.denominator)

def character_label(edict, weight):
    """Nebentypus discriminant label: D from s = (-1)^k prod d^{r_d} (exact)."""
    s = Fraction((-1) ** weight)
    for d, r in edict.items():
        s *= Fraction(d) ** r
    m = squarefree_kernel(s)
    if m == 0:
        return "trivial"
    # fundamental discriminant
    if m % 4 == 1:
        D = m
    else:
        D = 4 * m
    if D == 1:
        return "trivial"
    return f"chi_{D}"

def enumerate_edicts(N, weight, max_abs=3, trivial_only=True):
    """Enumerate holomorphic eta-quotients of given weight/level.
    Recursive over divisors with sum pruning + suffix residue pruning
    (Newman mod-24). Returns list of (edict, char)."""
    divs = divisors(N)
    target = 2 * weight
    out = []
    edict = {}
    # suffix achievable (sum, s1 mod 24, s2 mod 24) sets
    n = len(divs)
    suff = [None] * (n + 1)
    suff[n] = {(0, 0, 0)}
    for i in range(n - 1, -1, -1):
        d = divs[i]
        s = set()
        for (sm, a1, a2) in suff[i + 1]:
            for r in range(-max_abs, max_abs + 1):
                s.add((sm + r, (a1 + r * d) % 24, (a2 + r * (N // d)) % 24))
        suff[i] = s

    def rec(i, remaining, a1, a2):
        if i == len(divs):
            if remaining == 0 and edict:
                if not newman_ok(edict, N):
                    return
                if ligozat_orders(edict, N) is None:
                    return
                ch = character_label(edict, weight)
                if trivial_only and ch != "trivial":
                    return
                out.append((dict(edict), ch))
            return
        d = divs[i]
        nleft = len(divs) - i - 1
        for r in range(-max_abs, max_abs + 1):
            rem2 = remaining - r
            if abs(rem2) > max_abs * nleft:
                continue
            # suffix must be able to complete sum and Newman residues (O(1) lookup)
            if (rem2, (-a1 - r * d) % 24, (-a2 - r * (N // d)) % 24) not in suff[i + 1]:
                continue
            if r != 0:
                edict[d] = r
            rec(i + 1, rem2, (a1 + r * d) % 24, (a2 + r * (N // d)) % 24)
            if r != 0:
                del edict[d]

    rec(0, target, 0, 0)
    return out

def find_relations(edicts, nterms=60, max_support=4):
    """Given list of edicts (all same weight/level/char), compute q-expansions
    aligned on a common shift, find exact rational nullspace relations.
    Returns list of (coeffs, edict_indices) with small support, sorted by support."""
    # align shifts: all should have same shift for a relation to make sense
    # (different shifts -> multiply by q-power; only group equal shifts)
    from collections import defaultdict
    groups = defaultdict(list)
    for idx, e in enumerate(edicts):
        s, c = eta_quotient_qexp(e, nterms)
        groups[s].append((idx, c))
    relations = []
    for shift, items in groups.items():
        if len(items) < 2:
            continue
        M = sp.Matrix([c for _, c in items])
        ns = M.T.nullspace()
        for v in ns:
            nz = [(int(v[i]), items[i][0]) for i in range(len(items)) if v[i] != 0]
            if 2 <= len(nz) <= max_support:
                relations.append((nz, shift))
    # dedupe by canonical form: divide by gcd, flip so the coefficient of
    # the lexicographically smallest edict is positive
    def canon(nz):
        g = 0
        for c, _ in nz:
            g = gcd(g, abs(c))
        nn = tuple(sorted((c // g, i) for c, i in nz))
        smallest = min(nn, key=lambda x: x[1])
        if smallest[0] < 0:
            nn = tuple((-c, i) for c, i in nn)
        return nn
    seen = set()
    uniq = []
    for nz, shift in relations:
        cn = canon(nz)
        key = (shift, cn)
        if key not in seen:
            seen.add(key)
            uniq.append((cn, shift))
    uniq.sort(key=lambda t: (len(t[0]), t[0]))
    return uniq

if __name__ == "__main__":
    # sanity: MF-001 terms at N=16, weight 2, trivial char
    t1 = {2: 4, 4: -2, 8: -2, 16: 4}
    t2 = {1: 2, 2: -1, 4: -2, 8: 3, 16: 2}
    t3 = {1: 4, 2: -2, 8: -2, 16: 4}
    for t in (t1, t2, t3):
        print(t, "newman:", newman_ok(t, 16), "ligozat:", ligozat_orders(t, 16) is not None,
              "char:", character_label(t, 2))
    # enumerate N=16 weight 2 trivial
    es = enumerate_edicts(16, 2, max_abs=3, trivial_only=True)
    print("N=16 w=2 trivial edicts:", len(es))
