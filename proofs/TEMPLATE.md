# Proof Template: [Identity Name]

## Statement

[LaTeX statement of the identity]

**Parameters:** [parameter ranges]
**Type:** [eta-quotient / theta / etc.]
**Weight:** [k] | **Level:** [N] | **Character:** [chi]

## Proof

### Step 1: Modularity
Show that both sides are modular forms of weight `k` for `Γ₀(N)` with character `χ`.

[Details: eta-quotient criteria, or theta transformation laws]

### Step 2: Sturm Bound
The Sturm bound for `M_k(Γ₀(N))` is:
```
B = (k/12) · [SL₂(Z) : Γ₀(N)] = [computed value]
```

### Step 3: Verification
The q-expansions of both sides agree for `0 ≤ n ≤ B`:
```
[table or statement of coefficients]
```

Since both sides are modular forms of the same weight, level, and character,
and their Fourier coefficients agree up to the Sturm bound, they are identical.
∎

## Computational Verification

```python
# See scripts/sturm.py
from sturm import verify_identity
ok, B, n = verify_identity(A_edict, B_edict, k=k, N=N)
assert ok
```

## Novelty

- [ ] Checked against Martin–Ono eta-quotient database
- [ ] Checked against LMFDB
- [ ] Checked against Slater's list (if RR-type)
- [ ] Web search for the q-expansion (OEIS, etc.)

**Status:** [New / Known — with citation]

## References

[Relevant papers, database entries]
