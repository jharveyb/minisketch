# More optimizations to Minisketch decoding

This continues the "Optimizations" section of [math.md](math.md), using the same notation, and describes how decoding was made faster in 2026-10. [decode_explained.md](decode_explained.md) covers the same ideas for a reader with undergraduate discrete math, and [decode_coding_theory.md](decode_coding_theory.md) from the coding-theory side.

**Most of the Berlekamp-Massey steps can be skipped**

BM consumes the coefficients *s<sub>1</sub>, s<sub>2</sub>, s<sub>3</sub>, ...* one at a time. At each step it checks whether the current *L* correctly predicts the next coefficient, and corrects *L* if it doesn't. Thanks to the Frobenius endomorphism, *s<sub>2k</sub> = s<sub>k</sub><sup>2</sup>*. It turns out that this makes every prediction of an even coefficient correct, for any list of coefficients with this structure, not just the sketches of actual sets. So BM only has to look at the odd ones.

**Stopping early, but checking everything**

When the number of elements *n* is much smaller than the capacity *c*, *L* no longer changes once about *2n* coefficients have been processed. All the other coefficients merely have to agree with it, and we must still check that they do. If we didn't, a sketch whose first coefficients look like those of a small set would decode to the wrong answer instead of being rejected.

Checking is cheaper than solving, because each check is independent: coefficient *j* must satisfy *s<sub>j</sub> + s<sub>j-1</sub>l<sub>1</sub> + ... + s<sub>j-n</sub>l<sub>n</sub> = 0*. These are exactly the coefficients of degree *n* and higher of the product *S(M)L*, which can be computed with fast polynomial multiplication instead of one equation at a time.

Better still, only the odd coefficients are needed. Their sequence is *s<sub>1</sub> = m<sub>1</sub> + m<sub>2</sub> + ...*, *s<sub>3</sub> = m<sub>1</sub>m<sub>1</sub><sup>2</sup> + m<sub>2</sub>m<sub>2</sub><sup>2</sup> + ...*, *s<sub>5</sub> = m<sub>1</sub>(m<sub>1</sub><sup>2</sup>)<sup>2</sup> + ...*: for each element, a geometric series with ratio *m<sub>i</sub><sup>2</sup>*. Those series are cancelled by *(1 - m<sub>1</sub><sup>2</sup>y)(1 - m<sub>2</sub><sup>2</sup>y)...*, which, by the Frobenius endomorphism again, is *L* with each coefficient squared. One can show that once *L* has been found, checking the odd coefficients against this squared polynomial is equivalent to checking all of them. So the even coefficients are only ever computed for the short initial part.

**One trace for all the roots**

To find the roots of *L*, the polynomial is repeatedly split into two factors. A split uses the trace function *Tr(y) = y + y<sup>2</sup> + y<sup>4</sup> + ... + y<sup>2<sup>b-1</sup></sup>*. This function only ever produces *0* or *1*, and *Tr(a + b) = Tr(a) + Tr(b)*. For a random *β*, the polynomial *T = Tr(βx)* reduced modulo *L* is *0* at the roots *m* with *Tr(βm) = 0* and *1* at the others. So *gcd(L, T)* contains exactly the former roots, and *L* divided by it the latter. Computing *T* requires *x<sup>2</sup>, x<sup>4</sup>, x<sup>8</sup>, ...* modulo *L*, i.e. *b-1* squarings.

Every polynomial produced by splitting divides *L*. And for any divisor *G* of *L*, reducing a polynomial modulo *L* and then modulo *G* gives the same result as reducing it modulo *G* directly. So the squarings are computed once, for *L* itself, and every factor then only needs a single remainder. The traces for *b*-bit elements are linear combinations of *x<sup>2<sup>i</sup></sup>* modulo *L* with coefficients *β<sup>2<sup>i</sup></sup>*. If *β* is chosen as a fixed random value times an element of the 256-element subfield, those coefficients repeat every 8 steps, and the combination becomes much cheaper.

**A free test for undecodable sketches**

One more squaring gives *x<sup>2<sup>b</sup></sup>* modulo *L*. Every element *a* of *GF(2<sup>b</sup>)* satisfies *a<sup>2<sup>b</sup></sup> = a*. Therefore *x<sup>2<sup>b</sup></sup> - x* is the product of *(x - a)* over all field elements, and *L* divides it exactly when its roots are distinct field elements. So *x<sup>2<sup>b</sup></sup> = x mod L* holds if and only if the sketch is decodable, which we now know before any splitting starts.

**Delaying the modular reduction**

Multiplying two field elements consists of multiplying them as polynomials over *GF(2)* (one CLMUL instruction) and reducing the result modulo the field's polynomial. That reduction is linear, so a sum *a<sub>1</sub>b<sub>1</sub> + a<sub>2</sub>b<sub>2</sub> + ...*, which is what almost every loop in the decoder computes, can be reduced once at the end rather than after every product.

See [optimizations.md](optimizations.md) for the derivations, the measurements and the ideas that did not pay off.
