# How the faster decoder works

An explanation of the 2026-10 decode optimizations for a reader with
undergraduate discrete math. Optimizations 1–6 speed up decoding in general;
7–10 are the root-finding techniques for large, fairly full sketches (root
degree in the hundreds or more). [decode_coding_theory.md](decode_coding_theory.md)
covers the same ideas from the coding-theory side, and
[math_decode.md](math_decode.md) in the style of [math.md](math.md).

## The setup

Each set element is a 64-bit number, treated as an element of the finite field
GF(2⁶⁴). Addition in that field is XOR, and multiplication is polynomial
multiplication of the bit patterns, modulo a fixed degree-64 polynomial.

A sketch of a set {m₁,…,mₙ} with capacity c stores the power sums
sₖ = m₁ᵏ + m₂ᵏ + … + mₙᵏ for k = 1…2c. XOR-ing two sketches cancels the elements
they have in common, which is why sketches subtract sets.

Decoding has two steps:

1. **Berlekamp–Massey (BM)** finds the polynomial
   L(x) = (1 − m₁x)(1 − m₂x)⋯(1 − mₙx) from the power sums. The key fact is
   that the sequence s₁, s₂, … obeys a linear recurrence whose coefficients are
   exactly L's coefficients, like Fibonacci with n terms.
2. **Root finding** recovers the mᵢ from L.

The old decoder spent its time in two places: BM walked all 2c sketch entries
even when n ≪ c, and root finding redid the same expensive computation over and
over.

## Optimization 1: half of BM's steps are provably useless

The key identity is (a + b)² = a² + b², because the cross term 2ab is 0 when
1 + 1 = 0. So s₂ₖ = Σ mᵢ²ᵏ = (Σ mᵢᵏ)² = sₖ². Every even entry is determined by
an earlier one.

BM processes the sequence one term at a time and asks: does my current
recurrence predict this term? If not, it fixes the recurrence. The prediction
error is called the *discrepancy*. A classical theorem (Berlekamp, 1968) says
that for sequences with the s₂ₖ = sₖ² structure, the discrepancy at every even
position is exactly 0. A term that is always predicted correctly teaches
nothing, so those steps are skipped. That is half the work gone, with an
identical answer. (The unit tests compare against the unsimplified algorithm,
exhaustively on small fields.)

## Optimization 2: stop early, then check the rest in bulk

If the sketch holds n ≪ c elements, BM has learned L after about 2n terms.
After that, every remaining term just agrees with the recurrence. We can't just
stop, though. An adversary could build a sketch whose first 2n terms look like
n elements but whose later terms don't, and the old decoder would correctly
reject it. So the whole tail still has to be verified. But verifying differs
from learning:

- **Learning is sequential**: each step depends on the previous fix.
- **Verifying is embarrassingly parallel**: "does term j match Σ lᵢ·sⱼ₋ᵢ?" is
  independent for each j. It's a convolution of the short filter L with the
  long sequence, so fast multiplication applies.

There's a nice extra trick. The *odd* entries alone form a recurrence:
uⱼ = s₂ⱼ₊₁ = Σ mᵢ · (mᵢ²)ʲ, a sum of geometric sequences with ratios mᵢ². Its
recurrence polynomial is L with every coefficient squared. So the check only
needs the raw stored sketch values, not the reconstructed even ones.

The identity L·C = C² + L²·S_odd shows that this odd-only check is exactly
equivalent to the full check. Here S is the generating function of the sketch
and C = L·S. (The derivation is in [optimizations.md](optimizations.md).)

The bulk check uses a Karatsuba trick that turns 4 half-size multiplications
into 3. Applied recursively, that is about 33 instead of 65 multiplications per
position at n = 64.

**Result:** cost went from (capacity × n) sequential steps to about n²
sequential steps, plus a fast parallel check.

## Optimization 3: one hashing computation for the whole root-finding tree

Root finding works by recursive splitting, like randomized quicksort on the
roots. It needs a "random bit" for each root, but without knowing the roots.
The **trace** gives one: Tr(y) = y + y² + y⁴ + ⋯ + y^(2⁶³). By the same squaring
identity it's linear (Tr(a+b) = Tr(a) + Tr(b)), and it always outputs 0 or 1.
For a random β, the bit Tr(β·m) acts like which side of a random hyperplane m
lands on.

As a polynomial, T(x) = Tr(βx) mod L(x) takes value Tr(β·mᵢ) at each root mᵢ.
So gcd(L, T) collects exactly the roots with bit 0, and splits L in two.
Recursing on both halves gives all the roots.

Computing T needs x^(2ⁱ) mod L for i < 64: 63 "square, then reduce mod L" steps,
each costing about deg(L)². The old code redid this at every node of the
recursion tree.

**The observation:** every node's polynomial g divides the root polynomial f,
and reducing mod f then mod g is the same as reducing mod g. So the 63
squarings are computed once at the root, and each node just takes one
remainder. (For large polynomials, Optimizations 7–10 below replace these
remainders and most of the gcds.) As a bonus, x^(2⁶⁴) ≡ x mod f holds iff f has all its roots in the
field and none repeated. That's the field's version of Fermat's little
theorem, so one extra squaring tells immediately whether decoding can succeed.

## Optimization 4: choosing β cleverly

β^(2ⁱ) has 64 different values in general. But if β = γ·z with z inside the
small subfield GF(2⁸), then z^(2⁸) = z, so the powers repeat every 8 steps.
Combining all the precomputed squarings for 8 different β's then needs 128
multiplications per coefficient instead of 512. The β's still form a valid
"independent hyperplanes" basis, because the products γ·xˡ·zʲ form a tower
basis.

## Optimization 5: reduce once, not every time (lazy reduction)

Multiplying two field elements happens in two stages. First a carry-less
product (one CPU instruction, giving a 127-bit result). Then a reduction mod
the field polynomial, which is about two more instructions.

Reduction is linear: (a·b mod p) + (c·d mod p) = (a·b + c·d) mod p. So in any
sum of products, which describes most of the hot loops, the 127-bit products
are added first and reduced once at the end. That is 1 instruction per term
instead of 3. It's the same trick as computing Σ aᵢbᵢ mod q with integers:
reduce at the end, not after each term.

## Optimization 6: smaller algebra wins

- **Fraction-free Euclid.** In the gcd, instead of dividing by the leading
  coefficient at each step, scale the other polynomial so the leading terms
  cancel. That's like fraction-free Gaussian elimination. It costs one extra
  multiplication per coefficient, but saves a field inversion, which costs far
  more than a multiplication.
- **Dot-product long division.** Each quotient and remainder coefficient
  becomes one accumulated sum, which pairs with lazy reduction.

## Root finding for large sketches

When the sketch holds hundreds or thousands of elements, the root polynomial
f(x) = (x − m₁)(x − m₂)⋯(x − mₙ) (L with its coefficients reversed, so its
roots are the elements themselves) has degree d = n in the hundreds or
thousands. Then the 63 squarings mod f, the remainders at every node and the
gcds dominate decoding. The next four optimizations attack those.

## Optimization 7: squaring mod f using point values ("Montgomery on a coset")

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
modular squaring costs **2 FFTs of size N ≈ d**. The previous method needed
about 8 FFTs of size 2d.

The transforms got cheaper too:

- **Cheaper derivative.** The FFT works internally in a special "novel"
  polynomial basis, and in that basis the derivative is just an XOR pass.
- **No basis conversion.** We never convert back to ordinary coefficients
  in the middle of the squaring chain.

## Optimization 8: pass a half-computed trace down the tree

Each child polynomial g in the split tree divides f. So anything reduced
mod f, reduced further mod g, is correct mod g.

We exploit a tower of subfields, as in Optimization 4, but with the smaller
subfield GF(16) inside GF(2⁶⁴). Instead of the full trace, we compute
T(x) = Tr₆₄→₁₆(βx), a "partial trace" whose values lie in GF(16). One T
gives the split polynomial for 4 tree levels, at 3 squarings each, instead
of 63 squarings for a fresh trace. A few such T (for different β) cover all
the levels a tree of degree d actually reaches.

The tree hands each T (in point-value form) from parent to child. To shrink
the point set as the polynomials shrink:

- Split C into halves C₀ and C₁. Then Z_C = Z_C₀ · Z_C₁, and Z_C₀ = 1 at
  every point of C₁.
- Converting to the child's smaller representation is therefore one more
  Montgomery step, costing 2 half-size FFTs.

## Optimization 9: division by pointwise division

After a split g = g₁·g₂, we know g₁ and want g₂. Instead of long division,
compute g₂(w) = g(w)/g₁(w) at each point and interpolate back to
coefficients.

If the answer comes out monic (leading coefficient 1) with the right
degree, it really divides g. That check is also what makes Optimization 10
safe.

## Optimization 10: replace each gcd with a smaller decode

Let A be the set of roots where t = 1, for the node's split polynomial t.
Its power sums are

  P_k = Σ_{α∈A} αᵏ = Σ_α t(α)·αᵏ.

Write t(x) = Σⱼ tⱼxʲ. Expanding gives

  P_k = Σⱼ tⱼ·Q_{j+k},

where Q is the list of power sums of all of the node's roots. **So the
sketch of the subset A can be computed directly from the sketch of the
whole set.** At the root, Q is literally the sketch we were handed.

There's a bonus. Since t(α) ∈ {0, 1}, we have t(α)² = t(α), and so
P₂ₖ = Pₖ². That is the same "even terms are squares of earlier terms"
structure that real sketches have (Optimization 1). So P is a genuine
sketch of A, and we decode it with the same fast BM to get A's polynomial
directly.

Why this is faster:

- A has about n/2 roots, so its BM costs about a quarter of a full-size BM.
- That's about half of what Euclid's gcd costs on the degree-n polynomial.

Safety: if BM had too few terms, or produced something bogus, the
divisibility check from Optimization 9 catches it. We then retry on the
complement (power sums Q − P), and finally fall back to the gcd. Children
inherit P and Q − P, so nobody recomputes power sums.

## Everything is exact

Every change returns the same output as the old decoder for every input,
including garbage and adversarial sketches. Only the amount of work changed.
(Optimizations 7–10 use different split parameters, so they return the same
set of elements, possibly in a different order.)

See [optimizations.md](optimizations.md) for the derivations, the measurements
on an i7-4770 and a Raspberry Pi 5, and the ideas that were tried and dropped,
and [math.md](math.md) for the underlying construction.
