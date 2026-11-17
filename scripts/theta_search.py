"""Theta-function identity search.
Work in units of q^(1/4): series index i corresponds to exponent i/4.
theta3(q^d) = sum_n q^{d n^2}            -> index 4*d*n^2, coeff 1 (n=0) else 2
theta4(q^d) = sum_n (-1)^n q^{d n^2}     -> index 4*d*n^2, coeff (-1)^n*2 (n>0)
theta2(q^d) = sum_n q^{d (n+1/2)^2}     -> index d*(2n+1)^2, coeff 2
Products give weight-1 forms; 4th powers give weight-2 forms.
Find exact rational nullspace relations; report all, then classify vs classical.
"""
import numpy as np
import sympy as sp

NIDX = 400  # up to q^100

def theta3(d):
    s = np.zeros(NIDX + 1, dtype=object)
    n = 0
    while True:
        i = 4 * d * n * n
        if i > NIDX:
            break
        s[i] += 1 if n == 0 else 2
        n += 1
    return s

def theta4(d):
    s = np.zeros(NIDX + 1, dtype=object)
    n = 0
    while True:
        i = 4 * d * n * n
        if i > NIDX:
            break
        s[i] += (1 if n == 0 else 2) * (1 if n % 2 == 0 else -1)
        n += 1
    return s

def theta2(d):
    s = np.zeros(NIDX + 1, dtype=object)
    n = 0
    while True:
        i = d * (2 * n + 1) ** 2
        if i > NIDX:
            break
        s[i] += 2
        n += 1
    return s

def prod(a, b):
    # convolution, truncated
    out = np.zeros(NIDX + 1, dtype=object)
    ia = np.nonzero(a)[0]
    ib = np.nonzero(b)[0]
    for i in ia:
        ai = a[i]
        for j in ib:
            if i + j > NIDX:
                break
            out[i + j] += ai * b[j]
    return out

def nullspace_relations(named, max_support=4):
    """named: list of (name, series). Returns small-support exact relations
    with proper integer coefficients (scaled by lcm of denominators)."""
    from math import gcd
    from fractions import Fraction
    names = [n for n, _ in named]
    M = sp.Matrix([[int(x) for x in s] for _, s in named])
    ns = M.T.nullspace()
    out = []
    for v in ns:
        # scale to integers
        dens = [sp.denom(x) for x in v]
        L = 1
        for d in dens:
            L = L * int(d) // gcd(L, int(d))
        iv = [int(x * L) for x in v]
        g = 0
        for x in iv:
            g = gcd(g, abs(x))
        if g > 1:
            iv = [x // g for x in iv]
        if iv[0] < 0:
            iv = [-x for x in iv]
        nz = [(iv[i], names[i]) for i in range(len(names)) if iv[i] != 0]
        if 2 <= len(nz) <= max_support:
            out.append(nz)
    return out

if __name__ == "__main__":
    ds = [1, 2, 4, 8, 16]
    A = []
    for d in ds:
        A.append((f"th3(q^{d})", theta3(d)))
        A.append((f"th4(q^{d})", theta4(d)))
        A.append((f"th2(q^{d})", theta2(d)))
    print("=== weight 1/2: relations among 15 theta series ===")
    for r in nullspace_relations(A):
        print(r)
    print("=== weight 1: pairwise products ===")
    B = []
    for i in range(len(A)):
        for j in range(i, len(A)):
            B.append((f"{A[i][0]}*{A[j][0]}", prod(A[i][1], A[j][1])))
    print(f"{len(B)} products")
    for r in nullspace_relations(B, max_support=4):
        print(r)
    print("=== weight 2: fourth powers and mixed ===")
    C = []
    for d in ds:
        C.append((f"th3(q^{d})^4", prod(prod(theta3(d), theta3(d)), prod(theta3(d), theta3(d)))))
        C.append((f"th2(q^{d})^4", prod(prod(theta2(d), theta2(d)), prod(theta2(d), theta2(d)))))
        C.append((f"th4(q^{d})^4", prod(prod(theta4(d), theta4(d)), prod(theta4(d), theta4(d)))))
    for r in nullspace_relations(C, max_support=4):
        print(r)
