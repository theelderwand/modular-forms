"""Echelon + sparsify search for large edict sets (N=96 w=2).
Over GF(p): process edicts in order, maintain RREF basis; each dependency
gives a support set; exact QQ nullspace on the support finds sparse relations.
"""
import sys, time
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from fractions import Fraction
from math import gcd
from eta import eta_quotient_qexp
import sympy as sp

P = 1000003  # prime for echelon

def qexp_mod(edict, nterms):
    s, c = eta_quotient_qexp(edict, nterms)
    return s, [x % P for x in c]

def sparsify(idxs, series, max_support=4):
    """Exact QQ nullspace on the given edict indices; return sparse relations.
    series: list of coefficient lists aligned with idxs."""
    M = sp.Matrix([[int(x) for x in s] for s in series])
    ns = M.T.nullspace()
    out = []
    for v in ns:
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
        nz = [(iv[i], idxs[i]) for i in range(len(idxs)) if iv[i] != 0]
        if 2 <= len(nz) <= max_support:
            out.append(nz)
    return out

def echelon_search(edicts, nterms=40):
    """edicts: list of dicts, all same shift. Returns sparse relations (coeff, edict)."""
    # precompute integer series once; reduce mod P for echelon
    all_series = []
    for e in edicts:
        _, c = eta_quotient_qexp(e, nterms)
        all_series.append(c)
    # basis: pivot_col -> (row vector mod P, support dict edict_idx->coeff)
    basis = {}
    found = []
    seen_canon = set()
    for idx, e in enumerate(edicts):
        v = [x % P for x in all_series[idx]]
        support = {idx: 1}
        # reduce
        for col in sorted(basis):
            if v[col] != 0:
                prow, psup = basis[col]
                f = v[col] * pow(prow[col], -1, P) % P
                for j in range(col, nterms):
                    v[j] = (v[j] - f * prow[j]) % P
                for k, cv in psup.items():
                    support[k] = (support.get(k, 0) - f * cv) % P
        # find first nonzero
        piv = next((j for j in range(nterms) if v[j] != 0), None)
        if piv is None:
            # dependency: sparsify support
            sup = sorted(support)
            for nz in sparsify(sup, [all_series[i] for i in sup]):
                # canonical dedup: coeff of smallest index positive
                g = 0
                for c_, _ in nz:
                    g = gcd(g, abs(c_))
                nn = tuple(sorted((c_ // g, i) for c_, i in nz))
                if min(nn, key=lambda x: x[1])[0] < 0:
                    nn = tuple((-c_, i) for c_, i in nn)
                if nn not in seen_canon:
                    seen_canon.add(nn)
                    found.append(nn)
        else:
            inv = pow(v[piv], -1, P)
            v = [(x * inv) % P for x in v]
            support = {k: (cv * inv) % P for k, cv in support.items()}
            basis[piv] = (v, support)
    return found

if __name__ == "__main__":
    import json
    from collections import defaultdict
    from run_direction import is_substitution_of_known, is_primitive, verify_terms, canon_terms
    def parse(line):
        return {int(p.split(':')[0]): int(p.split(':')[1]) for p in line.strip().split(',') if p}
    def lcm_sup(ed):
        l = 1
        for d in ed:
            l = l * d // gcd(l, d)
        return l
    eds = [parse(l) for l in open('/tmp/e96w2.txt') if l.strip()]
    new96 = [e for e in eds if lcm_sup(e) == 96]
    groups = defaultdict(list)
    for e in new96:
        s, _ = eta_quotient_qexp(e, 0)
        groups[s].append(e)
    print(f"{len(new96)} new-at-96 edicts in {len(groups)} shift groups", flush=True)
    all_found = []
    for s0, g0 in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        t0 = time.time()
        rels = echelon_search(g0)
        print(f"shift {s0}: {len(g0)} edicts -> {len(rels)} sparse in {time.time()-t0:.1f}s", flush=True)
        for nz in rels:
            terms = [{"coeff": c, "edict": {str(k): v for k, v in g0[i].items()}} for c, i in nz]
            all_found.append(terms)
    # dedupe + screen + verify
    seen = set()
    out = []
    for terms in all_found:
        k = canon_terms(terms)
        if k in seen:
            continue
        seen.add(k)
        if not is_primitive(terms):
            continue
        sub, kid, m = is_substitution_of_known(terms)
        if sub:
            print(f"skip substitution of {kid} m={m}", flush=True)
            continue
        ok, sh = verify_terms(terms, 200)
        print(f"verified200={ok}: {terms}", flush=True)
        if ok:
            out.append({"id": None, "direction": "higher_level", "terms": terms,
                        "weight": 2, "level": 96, "character": "trivial",
                        "terms_verified": 200, "sparsity": len(terms),
                        "novelty_notes": "Primitive; lcm(support)=96; not a q->q^m substitution of MF-001..015 / W-001..030."})
    for i, r in enumerate(out, 1):
        r["id"] = f"L96-{i:03d}"
    json.dump(out, open("/home/hatch/workspace/modular-forms/identities/candidates_churn_level96w2.json", "w"), indent=1)
    print(f"WROTE {len(out)} candidates")
