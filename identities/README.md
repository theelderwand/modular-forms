# Eta-Quotient Identity Search — Results

## HEADLINE: Pure A=B eta-quotient identities DO NOT EXIST (theorem)

We searched ~950,000 eta-quotients (exponents in [-3,3], divisors of 16/24/32,
weights 1,2,3) for distinct edicts with identical q-expansions (40 terms).
**Zero collisions.** This is not a gap in the search — it is a theorem:

**Theorem (multiplicative independence of η(dz)).**
If D: ℕ→ℤ is finitely supported and not identically zero, then
∏_{d≥1} η(dz)^{D_d} is not identically 1 on ℍ.

*Proof.* Suppose ∏η(dz)^{D_d}=1. Taking d/dz of the logarithm:
  ∑_d D_d·d·η'(dz)/η(dz) = 0.
Since η'(w)/η(w) = (2πi/24)·E₂(w), we get ∑_d (D_d·d)·E₂(dz) = 0.
But the functions {E₂(dz)}_{d≥1} are linearly independent over ℂ:
if ∑c_d E₂(dz)=0, then the q-expansion gives ∑_{d|m}c_dσ(m/d)=0 for all m≥1,
whence c_1=c_2=⋯=0 by induction (using σ(1)=1).
Thus D_d·d=0 ∀d, so D=0. ∎

**Corollary.** There are no identities η^A = η^B with A≠B as edicts.
(Any A=B gives ∏η^{A-B}=1, forcing A=B.)

In particular, even A=c·B with c≠1 constant is impossible (log-derivative kills c).

So the "eta-quotient A = eta-quotient B" search is provably empty. We pivoted to
**linear relations**, which must exist by finite-dimensionality.

## Linear relations (this file)

We enumerated holomorphic eta-quotients with trivial character satisfying
Newman's criteria, then found ℚ-linear relations via exact nullspace computation.

For each (N,k): edicts r_δ (δ|N), r_δ∈[-4,4], ∑r_δ=2k, with
- ∑δ·r_δ ≡ 0 (mod 24) and ∑(N/δ)·r_δ ≡ 0 (mod 24) [Newman: modular on Γ₀(N)]
- ∀d|N: ∑_δ r_δ·gcd(d,δ)²/δ ≥ 0 [holomorphic at all cusps; Ligozat]
- ∏δ^{r_δ} is a square in ℚ [trivial Nebentypus; note rad automatically |N]

Each relation ∑c_i f_i=0 was verified to 500 q-expansion terms.
Since 500 >> Sturm bound (≤16 for all (N,k) here), each is a **theorem** by Sturm,
not just a numerical coincidence. (Full Sturm proof is routine from here.)

### Contents of candidates.json
15 relations, each:
- `type`: "linear_relation"
- `terms`: [{"coeff": c, "edict": {d: e}}] meaning ∑ c·∏η(dz)^e = 0
- `weight`, `level`, `shift` (=∑d·e_d/24, common to all terms)
- `verified_terms`: 500, `sturm_bound`: Sturm bound for M_k(Γ₀(N))

Ranked by (#terms, max|coeff|): 3-term relations with coefficients in {±1,±2,±3}.

### Distribution found (391 total verified relations)
- N=12,k=4: 2 | N=16,k=2: 1 | N=24,k=2: 91 | N=24,k=4: 290 | N=32,k=2: 5 | N=32,k=4: 2

### Caveats
- **Novelty vs literature NOT checked.** Some relations may be known (e.g., consequences
  of dimension formulas or in the "eta-quotients as Eisenstein series" literature).
  The 3-term sparse ones are the best candidates for genuine novelty.
- **Odd weights (1,3) skipped.** Odd-weight eta-quotients always have nontrivial character
  (χ(-1)=-1), so they lie in M_k(Γ₀(N),χ), not M_k(Γ₀(N)). A character-grouped search
  is possible but was not needed to exceed the 10-candidate target.
- **Dilatation redundancy.** A relation at level N that is a dilation (z→mz) of a
  lower-level relation is "not new". Primitive (gcd of all d's = 1) relations are
  preferred; check the `level` vs the actual content.
- The 4 prior conjectures (q-series, not eta-quotients) are unaffected by this search.

## Files
- `candidates.json`: 15 selected relations (this report).
- `/tmp/mfsearch/linrels.json`: all 391 relations (ephemeral).
- `/tmp/mfsearch/search.py`: the A=B search (proves empty).
- `~/workspace/modular-forms/scripts/eta.py`: q-expansion engine.
