"""Direction 2 pipeline: load C-enumerated edicts, find sparse relations.
Usage: python dir2.py <edict_file> <N> <weight> <tag>
Writes candidates to identities/candidates_churn_<tag>.json
"""
import sys, json
sys.path.insert(0, "/home/hatch/workspace/modular-forms/scripts")
from churn import find_relations, character_label
from run_direction import is_substitution_of_known, is_primitive, verify_terms

def parse(line):
    return {int(p.split(':')[0]): int(p.split(':')[1]) for p in line.strip().split(',') if p}

def main():
    fn, N, w, tag = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    eds = [parse(l) for l in open(fn) if l.strip()]
    print(f"[{tag}] loaded {len(eds)} edicts", flush=True)
    # group by character (all should be trivial; check)
    from collections import defaultdict
    bychar = defaultdict(list)
    for e in eds:
        bychar[character_label(e, w)].append(e)
    print(f"[{tag}] chars: { {k: len(v) for k, v in bychar.items()} }", flush=True)
    out = []
    for ch, group in bychar.items():
        if ch != "trivial":
            print(f"[{tag}] skipping non-trivial char {ch}", flush=True)
            continue
        rels = find_relations(group, nterms=40, max_support=4)
        print(f"[{tag}] {len(rels)} raw relations", flush=True)
        for nz, shift in rels:
            terms = [{"coeff": c, "edict": {str(k): v for k, v in group[i].items()}} for c, i in nz]
            if not is_primitive(terms):
                continue
            sub, kid, m = is_substitution_of_known(terms)
            if sub:
                print(f"[{tag}] skip substitution of {kid} m={m}", flush=True)
                continue
            ok, sh = verify_terms(terms, 200)
            print(f"[{tag}] verified200={ok}: {terms}", flush=True)
            if ok:
                out.append({"id": None, "direction": tag, "terms": terms,
                            "weight": w, "level": N, "character": ch,
                            "terms_verified": 200, "sparsity": len(terms),
                            "novelty_notes": "Primitive; not a q->q^m substitution of MF-001..015 / W-001..030."})
    for i, r in enumerate(out, 1):
        r["id"] = f"L-{i:03d}"
    json.dump(out, open(f"/home/hatch/workspace/modular-forms/identities/candidates_churn_{tag}.json", "w"), indent=1)
    print(f"[{tag}] WROTE {len(out)} candidates")

if __name__ == "__main__":
    main()
