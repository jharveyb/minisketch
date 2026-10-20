# Fast root finding for large sketches: an introduction

This note explains the root-finding optimizations of the 2026-10 "coset
Montgomery / power-sum" series (see [optimizations.md](optimizations.md)
for costs and measurements) for a reader with a discrete-math background.
The same material is covered from a coding-theory angle in
[root-finding-coding-theory.md](root-finding-coding-theory.md), and in the
style of [math.md](math.md) in [math-root-finding.md](math-root-finding.md).

## The setup

A sketch of a set M = {m₁, …, mₙ} is the list of **power sums**
s_k = m₁ᵏ + m₂ᵏ + … + mₙᵏ, computed in the finite field GF(2⁶⁴). Elements
of that field are 64-bit strings. Addition is XOR, and multiplication is
polynomial multiplication modulo a fixed degree-64 polynomial.

Decoding has two steps:

1. **Berlekamp–Massey (BM)** turns the power sums into a polynomial
   f(x) = (x − m₁)(x − m₂)…(x − mₙ) of degree d = n.
2. **Root finding** recovers the mᵢ from f.

You can't just try all 2⁶⁴ values of x, so step 2 is the interesting one.

## How root finding works (unchanged idea, new arithmetic)

The **trace** is the map Tr(y) = y + y² + y⁴ + … + y^(2⁶³). It has three
properties we use:

- It always outputs 0 or 1.
- It is linear: Tr(a + b) = Tr(a) + Tr(b).
- It sends exactly half the field to 0.

Pick a constant β. Then the polynomial t(x) = Tr(βx) mod f is 0 at roughly
half of f's roots and 1 at the rest. So gcd(f, t) is the product of
(x − mᵢ) over the roots where t is 0. That splits f into two polynomials of
about half the size. Recurse with a new β at each level of the split tree
until every piece has degree 1.

The expensive part is computing t(x). It needs x, x², x⁴, …, x^(2⁶³) mod f,
which is 63 repeated squarings modulo a degree-d polynomial. This chain is
called the **Frobenius orbit**.

Before this series, the orbit was already computed once, at the root, and
β was chosen so that traces for several tree levels could be read off that
one orbit.

## Change 1: squaring mod f using point values ("Montgomery on a coset")

**Coefficients vs. point values.** A polynomial of degree < N can be stored
either as its coefficients or as its values at N points. In point-value
form, multiplying two polynomials is just multiplying the values pointwise.
An FFT converts between the two forms in O(N log N).

The difficulty is reducing mod f, which is natural in coefficient form. We
borrow an idea from integer arithmetic: **Montgomery multiplication**. For
integers, Montgomery reduces mod p by dividing by R = 2ᵏ, because dividing
by a power of two is easy. Here we use a polynomial R that is easy in
point-value form:

- **The choice of R.** Let C be the set of N evaluation points. Let
  R = Z_C(x) = ∏_{w∈C} (x − w), so R is zero at every point of C.
- **Step 1: find the correction q.** For a product a = r², we want a
  polynomial q such that a + q·f is divisible by R. That means
  q ≡ a/f mod R. "Mod R" is the same as "the values on C", so we get q by
  pointwise division: q(w) = a(w)/f(w) at each point w.
- **Step 2: divide by R.** The result is r_new = (a + q·f)/R. Its values on
  C would be 0/0, since R vanishes there.

**The trick to resolve 0/0.** Write h = a + q·f = r_new·R and
differentiate both sides:

  h′ = r_new′·R + r_new·R′.

At a point w of C, R(w) = 0, so h′(w) = r_new(w)·R′(w). Now two facts about
characteristic 2 help:

- **R′ = 1.** C is chosen as an affine subspace (bit strings forming a
  shifted "subspace" under XOR). Then R has the form x^(2^m) + … + x + c,
  with only power-of-two exponents. The derivative of x^(2^j) is
  2^j·x^(2^j − 1), and 2^j = 0 in characteristic 2 for j ≥ 1. Only the x
  term survives, so R′ = 1.
- **(r²)′ = 2rr′ = 0.** So a′ = 0, and h′ = q′·f + q·f′.

Putting these together:

  r_new(w) = q′(w)·f(w) + q(w)·f′(w).

Everything there is pointwise except q′. Getting q′ costs one inverse FFT
(to get q's coefficients), a cheap differentiation, and one FFT back. So one
modular squaring costs **2 FFTs of size N ≈ d**. The old method needed about
8 FFTs of size 2d.

The transforms got cheaper too:

- **Cheaper derivative.** The FFT works internally in a special "novel"
  polynomial basis, and in that basis the derivative is just an XOR pass.
- **No basis conversion.** We never convert back to ordinary coefficients
  in the middle of the squaring chain.

## Change 2: pass a half-computed trace down the tree

Each child polynomial g in the split tree divides f. So anything reduced
mod f, reduced further mod g, is correct mod g.

We exploit a tower of subfields: GF(16) sits inside GF(2⁶⁴). Instead of the
full trace, we compute T(x) = Tr₆₄→₁₆(βx), a "partial trace" whose values
lie in GF(16). One T gives the split polynomial for 4 tree levels, at 3
squarings each, instead of 63 squarings for a fresh trace. A few such T
(for different β) cover all the levels a tree of degree d actually reaches.

The tree hands each T (in point-value form) from parent to child. To shrink
the point set as the polynomials shrink:

- Split C into halves C₀ and C₁. Then Z_C = Z_C₀ · Z_C₁, and Z_C₀ = 1 at
  every point of C₁.
- Converting to the child's smaller representation is therefore one more
  Montgomery step, costing 2 half-size FFTs.

## Change 3: division by pointwise division

After a split g = g₁·g₂, we know g₁ and want g₂. Instead of long division,
compute g₂(w) = g(w)/g₁(w) at each point and interpolate back to
coefficients.

If the answer comes out monic (leading coefficient 1) with the right
degree, it really divides g. That check is also what makes Change 4 safe.

## Change 4: replace each gcd with a smaller decode

Let A be the set of roots where t = 1. Its power sums are

  P_k = Σ_{α∈A} αᵏ = Σ_α t(α)·αᵏ.

Write t(x) = Σⱼ tⱼxʲ. Expanding gives

  P_k = Σⱼ tⱼ·Q_{j+k},

where Q is the list of power sums of all of g's roots. **So the sketch of
the subset A can be computed directly from the sketch of the whole set.**
At the root, Q is literally the syndromes we were handed.

There's a bonus. Since t(α) ∈ {0, 1}, we have t(α)² = t(α), and so
P₂ₖ = Pₖ². That is the same "even terms are squares of earlier terms"
structure that real sketches have. So P is a genuine sketch of A, and we
decode it with the same fast BM to get A's polynomial directly.

Why this is faster:

- A has about n/2 roots, so its BM costs about a quarter of a full-size BM.
- That's about half of what Euclid's gcd costs on the degree-n polynomial.

Safety: if BM had too few terms, or produced something bogus, the
divisibility check from Change 3 catches it. We then retry on the
complement (power sums Q − P), and finally fall back to the gcd. Children
inherit P and Q − P, so nobody recomputes power sums.
