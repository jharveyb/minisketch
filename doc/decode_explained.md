# How the faster decoder works

An explanation of the 2026-10 decode optimizations for a reader with
undergraduate discrete math. [decode_coding_theory.md](decode_coding_theory.md)
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
remainder. As a bonus, x^(2⁶⁴) ≡ x mod f holds iff f has all its roots in the
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

## Everything is exact

Every change returns the same output as the old decoder for every input,
including garbage and adversarial sketches. Only the amount of work changed.

See [optimizations.md](optimizations.md) for the derivations, the measurements
on an i7-4770 and a Raspberry Pi 5, and the ideas that were tried and dropped,
and [math.md](math.md) for the underlying construction.
