"""Runner for churn directions 2 and 3. Usage:
  python run_direction.py level64_96
  python run_direction.py weight5
Writes JSON to ~/workspace/modular-forms/identities/candidates_churn_<tag>.json
"""
import sys, json, time
from math import gcd
from fractions import Fraction
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from churn import enumerate_edicts, find_relations, character_label
from eta import eta_quotient_qexp

KNOWN = []
for fn in ("candidates.json", "candidates_wide.json"):
    try:
        KNOWN.extend(json.load(open(f"/home/hatch/workspace/modular-forms/identities/{fn}")))
    except FileNotFoundError:
        pass

def known_signatures():
    sigs = []
    for x in KNOWN:
        terms = x["terms"]
        sigs.append(frozenset((t["coeff"], tuple(sorted((int(k), v) for k, v in t["edict"].items()))) for t in terms))
    return sigs

def is_primitive(terms):
    g = 0
    for t in terms:
        for d in t["edict"]:
            g = gcd(g, int(d))
    return g == 1

def canon_terms(terms):
    """Sign-canonical normalized form of a relation's term list.
    Divide coeffs by gcd; flip so the coefficient of the lexicographically
    smallest edict is positive. Returns frozenset of (coeff, edict-tuple)."""
    g = 0
    for t in terms:
        g = gcd(g, abs(t["coeff"]))
    nn = []
    for t in terms:
        nn.append((t["coeff"] // g, tuple(sorted((int(k), v) for k, v in t["edict"].items()))))
    # canonical sign: coefficient of smallest edict must be positive
    smallest = min(nn, key=lambda x: x[1])
    if smallest[0] < 0:
        nn = [(-c, e) for c, e in nn]
    return frozenset(nn)

def is_substitution_of_known(terms):
    """Check whether terms == scale of a known relation by q->q^m for m>=1."""
    new_n = canon_terms(terms)
    for x in KNOWN:
        known_terms = x["terms"]
        # try scalings m=1..8 (m=1 catches exact duplicates / oldforms at higher level)
        for m in range(1, 9):
            scaled = []
            for t in known_terms:
                ed = {str(int(k) * m): v for k, v in t["edict"].items()}
                scaled.append({"coeff": t["coeff"], "edict": ed})
            if canon_terms(scaled) == new_n:
                return True, x.get("id", "?"), m
    return False, None, None

def verify_terms(terms, nterms=200):
    """Verify sum coeff*eta == 0 to nterms. Returns (ok, shift)."""
    from collections import defaultdict
    total = None
    shift0 = None
    for t in terms:
        ed = {int(k): v for k, v in t["edict"].items()}
        s, c = eta_quotient_qexp(ed, nterms)
        if shift0 is None:
            shift0 = s
            total = [t["coeff"] * x for x in c]
        else:
            if s != shift0:
                return False, None
            total = [a + t["coeff"] * b for a, b in zip(total, c)]
    return all(v == 0 for v in total), shift0

def run_case(N, weight, trivial_only, tag, max_abs=4):
    t0 = time.time()
    es = enumerate_edicts(N, weight, max_abs=max_abs, trivial_only=trivial_only)
    print(f"[{tag}] N={N} w={weight} trivial={trivial_only}: {len(es)} edicts ({time.time()-t0:.1f}s)", flush=True)
    # group by character
    from collections import defaultdict
    bychar = defaultdict(list)
    for e, ch in es:
        bychar[ch].append(e)
    out = []
    for ch, eds in bychar.items():
        rels = find_relations(eds, nterms=40, max_support=4)
        print(f"[{tag}] char {ch}: {len(eds)} edicts, {len(rels)} raw relations", flush=True)
        for nz, shift in rels:
            terms = [{"coeff": c, "edict": {str(k): v for k, v in eds[i].items()}} for c, i in nz]
            if not is_primitive(terms):
                continue
            sub, kid, m = is_substitution_of_known(terms)
            if sub:
                print(f"[{tag}] skip: substitution of {kid} (m={m})", flush=True)
                continue
            ok, sh = verify_terms(terms, 200)
            print(f"[{tag}] candidate w={weight} N={N} char={ch} verified200={ok}: {terms}", flush=True)
            if ok:
                out.append({
                    "id": None, "direction": tag,
                    "terms": terms, "weight": weight, "level": N,
                    "character": ch, "terms_verified": 200,
                    "sparsity": len(terms),
                    "novelty_notes": "Primitive; not a q->q^m substitution of MF-001..015 / W-001..030.",
                })
    return out

if __name__ == "__main__":
    mode = sys.argv[1]
    results = []
    if mode == "level64_96":
        for N in (64, 96):
            for w in (2, 4):
                results.extend(run_case(N, w, True, "higher_level"))
    elif mode == "weight5":
        for N in (12, 24):
            results.extend(run_case(N, 5, False, "weight5"))
    for i, r in enumerate(results, 1):
        r["id"] = f"C-{i:03d}"
    out = f"/home/hatch/workspace/modular-forms/identities/candidates_churn_{mode}.json"
    json.dump(results, open(out, "w"), indent=1)
    print(f"WROTE {len(results)} candidates to {out}")
