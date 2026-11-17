# Modular Forms: Original Identities

A collection of **new, proven identities** in the theory of modular forms —
eta-quotients, theta functions, and related q-series.

Each identity in this repository is:
- **New**: not found in the standard literature (Slater, Martin–Ono, LMFDB, etc.)
- **Proven**: with a complete, lean proof (Sturm bound, modularity criteria, or elementary methods)
- **General**: parametric families where possible, not just isolated numerical curiosities

## Structure

- `identities/` — Statement of each identity (LaTeX + q-expansion data)
- `proofs/` — Lean proofs (mathematical argument + Sturm-bound verification code)
- `pages/` — Rendered HTML pages for each identity
- `scripts/` — Computational tools (eta-quotient engine, Sturm-bound verifier)

## Method

New identities are discovered via systematic computational search over
eta-quotients, then **proven** using the Sturm bound: if two modular forms
of the same weight, level, and character agree on sufficiently many
Fourier coefficients, they are identical.

## Identities

| # | Name | Type | Status |
|---|------|------|--------|
| # | Name | Type | Weight | Level | Status |
|---|------|------|--------|-------|--------|
| MF-001 | Eta-quotient linear relation | 3-term | 2 | 16 | Proven (Sturm) |
| MF-002 | Eta-quotient linear relation | 3-term | 2 | 24 | Proven (Sturm) |
| MF-003 | Eta-quotient linear relation | 3-term | 2 | 24 | Proven (Sturm) |
| MF-004 | Eta-quotient linear relation | 3-term | 2 | 24 | Proven (Sturm) |
| MF-005 | Eta-quotient linear relation | 3-term | 4 | 24 | Proven (Sturm) |
| MF-006 | Eta-quotient linear relation | 3-term | 2 | 32 | Proven (Sturm) |
| MF-007 | Eta-quotient linear relation | 3-term | 2 | 32 | Proven (Sturm) |
| MF-008 | Eta-quotient linear relation | 3-term | 2 | 32 | Proven (Sturm) |
| MF-009 | Eta-quotient linear relation | 3-term | 2 | 32 | Proven (Sturm) |
| MF-010 | Eta-quotient linear relation | 3-term | 4 | 32 | Proven (Sturm) |
| MF-011 | Eta-quotient linear relation | 3-term | 4 | 32 | Proven (Sturm) |
| MF-012 | Eta-quotient linear relation | 3-term | 2 | 24 | Proven (Sturm) |
| MF-013 | Eta-quotient linear relation | 3-term | 2 | 24 | Proven (Sturm) |
| MF-014 | Eta-quotient linear relation | 3-term | 2 | 24 | Proven (Sturm) |
| MF-015 | Eta-quotient linear relation | 3-term | 4 | 24 | Proven (Sturm) | |

## References

- Martin–Ono databases of eta-quotients
- LMFDB (L-functions and Modular Forms Database)
- Slater's list of Rogers–Ramanujan-type identities
