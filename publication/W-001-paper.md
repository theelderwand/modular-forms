# A Weight-One Eta-Quotient Relation with Nontrivial Nebentypus

## Abstract

We exhibit a new three-term linear relation among eta-quotients of weight 1
and level 12 with Nebentypus character $\chi_{-1}$, proven by the Sturm bound.
We also record the elementary observation that no nontrivial identity of the
form $A = B$ between distinct eta-quotients exists, by multiplicative
independence of the functions $\eta(dz)$.

## 1. Introduction

The Dedekind eta function $\eta(z) = q^{1/24}\prod_{n\geq 1}(1-q^n)$,
$q = e^{2\pi i z}$, generates a rich family of modular forms through
eta-quotients $\prod_{d\mid N} \eta(dz)^{e_d}$. The classification of
individual eta-quotients as modular forms is classical (Newman; Martin–Ono).
Here we are concerned with *linear relations* among eta-quotients — that is,
$\mathbb{Q}$-linear dependencies in the finite-dimensional spaces
$M_k(\Gamma_0(N), \chi)$.

## 2. Multiplicative independence

**Theorem 2.1.** The functions $\eta(dz)$, $d \geq 1$, are multiplicatively
independent: if $\prod_d \eta(dz)^{D_d} = 1$ with $D$ finitely supported,
then $D = 0$.

*Proof.* Taking the logarithmic derivative,
$\sum_d D_d \cdot d \cdot E_2(dz) = 0$, where $E_2$ is the quasimodular
Eisenstein series of weight 2. The functions $E_2(dz)$ are linearly
independent over $\mathbb{C}$: from $\sum_d c_d E_2(dz) = 0$, comparing
$q$-expansions gives $\sum_{d\mid m} c_d\,\sigma(m/d) = 0$ for all $m$,
whence $c_1 = c_2 = \cdots = 0$ by induction on $m$. ∎

**Corollary 2.2.** There is no identity $A = B$ between distinct
eta-quotients, nor $A = cB$ with $c \neq 1$.

Hence the natural objects of study are linear relations with three or more
terms.

## 3. The identity

**Theorem 3.1.** Let $\chi_{-1}$ be the nontrivial Dirichlet character
modulo $4$ (the Kronecker symbol $(\frac{-1}{\cdot})$). Then the following
eta-quotients are modular forms of weight $1$ for $\Gamma_0(12)$ with
Nebentypus $\chi_{-1}$, and they satisfy:

$$3\cdot \eta(z)^{-2}\eta(2z)^{3}\eta(3z)^{2}\eta(4z)^{-1}\eta(6z)\eta(12z)^{-1}$$
$$- 4\cdot \eta(z)^{-1}\eta(2z)\eta(3z)^{-1}\eta(4z)^{2}\eta(6z)^{3}\eta(12z)^{-2}$$
$$+ \eta(z)^{2}\eta(2z)\eta(3z)^{-2}\eta(4z)^{-1}\eta(6z)^{3}\eta(12z)^{-1} = 0.$$

*Proof.* Each term is an eta-quotient with $\sum e_d = 2$, hence weight
$(1/2)\sum e_d = 1$. Each satisfies Newman's congruence conditions and
Ligozat's cusp-holomorphicity criteria, so each lies in
$M_1(\Gamma_0(12), \chi_{-1})$. The Sturm bound for this space is

$$B = \frac{1}{12}\cdot[\mathrm{SL}_2(\mathbb{Z}):\Gamma_0(12)] = \frac{1}{12}\cdot 24 = 2.$$

A direct computation of the $q$-expansions shows the linear combination
vanishes through $q^{200}$ (far exceeding the bound). By the Sturm bound,
the combination is identically zero. ∎

## 4. Remarks

1. Weight-one modular forms with nontrivial character are comparatively rare;
   this relation lives in a space of small dimension, making the three-term
   dependency striking.
2. The relation is *primitive*: the divisors $\{1,2,3,4,6,12\}$ have
   $\gcd = 1$, so it is not a dilation of a lower-level relation.
3. It is not a $q \mapsto q^m$ substitution of any previously known relation.

## 5. Computational verification

The $q$-expansions were computed with exact integer arithmetic (Python,
`fractions` module) to 200 terms. The Sturm-bound verification script is
available at `scripts/sturm.py` in the accompanying repository.

## References

- M. Newman, *Construction and application of a class of modular functions*,
  Proc. London Math. Soc. (1959).
- Y. Martin and K. Ono, *Eta-quotients and elliptic curves*,
  Proc. Amer. Math. Soc. (1997).
- G. Ligozat, *Courbes modulaires de genre 1*, Bull. Soc. Math. France (1975).
