# More optimizations to Minisketch decoding

This continues the "Optimizations" section of [math.md](math.md), using the same notation, and describes how decoding was made faster in 2026-10, including the techniques used to find the roots when a sketch holds many elements. [decode_explained.md](decode_explained.md) covers the same ideas for a reader with undergraduate discrete math, and [decode_coding_theory.md](decode_coding_theory.md) from the coding-theory side.

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

One more squaring gives *x<sup>2<sup>b</sup></sup>* modulo *L*. Every element *a* of *GF(2<sup>b</sup>)* satisfies *a<sup>2<sup>b</sup></sup> = a*. Therefore *x<sup>2<sup>b</sup></sup> - x* is the product of *(x - a)* over all field elements, and *L* divides it exactly when its roots are distinct field elements. So *x<sup>2<sup>b</sup></sup> = x mod L* holds if and only if the roots of *L* are distinct field elements. Given that Berlekamp-Massey succeeded within the capacity, that is exactly when the sketch decodes (and 0 is then never a root), which we now know before any splitting starts.

**Delaying the modular reduction**

Multiplying two field elements consists of multiplying them as polynomials over *GF(2)* (one CLMUL instruction) and reducing the result modulo the field's polynomial. That reduction is linear, so a sum *a<sub>1</sub>b<sub>1</sub> + a<sub>2</sub>b<sub>2</sub> + ...*, which is what almost every loop in the decoder computes, can be reduced once at the end rather than after every product.

**Root finding for large sketches**

Write *f* for *L* with its coefficients reversed: its roots are the set elements themselves (see "Avoiding inversions" in [math.md](math.md)), and it has degree *n*. When *n* is in the hundreds or more, the *b - 1* squarings modulo *f*, the remainders at every factor and the gcds dominate decoding. The following three ideas replace most of that work.

**Squaring by evaluation**

A polynomial of degree less than *N* is determined by its values at *N* points, and products are pointwise in that form. The FFT used here evaluates on an affine subspace *C* of the field (*N = 2<sup>k</sup>* elements closed under XOR, shifted by a constant). The polynomial *Z(x)*, the product of *(x - w)* over all *w &isin; C*, has a special shape: *x<sup>N</sup> + ... + x + c*, with only power-of-two exponents.

We keep every residue *r* in the scaled form *rZ mod f*, as is done for integers in Montgomery multiplication. To reduce a product *a*:

* Find *q* with *a + qf* divisible by *Z*. Since *Z* vanishes on *C*, that just means *q(w) = a(w)/f(w)* at every point *w* of *C*.
* The reduced value is *(a + qf)/Z*. That is 0/0 on *C*, but differentiating *a + qf = rZ* gives *(a + qf)' = r'Z + rZ'*, so at the points of *C*, *r(w) = (a + qf)'(w) / Z'(w)*.

In fields where *x + x = 0*, the derivative of every *x<sup>2<sup>j</sup></sup>* with *j &ge; 1* is zero. So *Z' = 1*, and the derivative of a square is zero too. For squaring this leaves:

* *r(w) = q'(w)f(w) + q(w)f'(w)*

Only *q'* needs a transform: an inverse FFT to get *q*, a differentiation, and an FFT back. A modular squaring costs two transforms of the size of *f*.

**Traces for every factor from one partial trace**

The factors all divide *f*, so a residue mod *f* is also one mod each factor. Write the trace in two steps through the subfield *GF(16)*:

* *T = Tr<sub>GF(2<sup>b</sup>)/GF(16)</sub>(&beta;x)*, whose values at the roots lie in *GF(16)*.
* *Tr(&beta;z<sup>k</sup>x) = Tr<sub>GF(16)/GF(2)</sub>(z<sup>k</sup>T) = z<sup>k</sup>T + (z<sup>k</sup>T)<sup>2</sup> + (z<sup>k</sup>T)<sup>4</sup> + (z<sup>k</sup>T)<sup>8</sup>*, for *z* generating *GF(16)*.

The splits of four consecutive levels therefore need only three squarings each of *T* reduced mod the factor. *T* is passed from each factor to its own factors, still in evaluated form. When the factors get small enough, the evaluation set is halved: *C* splits into halves *C<sub>0</sub>* and *C<sub>1</sub>*, *Z<sub>C</sub> = Z<sub>C<sub>0</sub></sub>Z<sub>C<sub>1</sub></sub>*, and *Z<sub>C<sub>0</sub></sub> = 1* on *C<sub>1</sub>*, so the change of representation is one more reduction step. Once a factor *g* is split into *g<sub>1</sub>g<sub>2</sub>*, *g<sub>2</sub>* follows from *g<sub>2</sub>(w) = g(w)/g<sub>1</sub>(w)* without any long division.

**Sketches of halves**

Let *t = Tr(&beta;x) mod g*, and let *A* be the roots of *g* where *t = 1*. The power sums of *A* are:

* *p<sub>k</sub> = &sum;<sub>m&isin;A</sub> m<sup>k</sup> = &sum;<sub>m</sub> t(m)m<sup>k</sup> = &sum;<sub>j</sub> t<sub>j</sub>q<sub>j+k</sub>*

where *q<sub>k</sub>* are the power sums of all roots of *g*. For *f* itself, those are the sketch's own coefficients. Because *t(m)* is 0 or 1, *t(m)<sup>2</sup> = t(m)*, and so *p<sub>2k</sub> = p<sub>k</sub><sup>2</sup>*, just like a real sketch. The list *p* is therefore a sketch of *A*, computed from the sketch of the whole set, and Berlekamp-Massey decodes it directly into the polynomial of *A*. That's a problem of about half the size, cheaper than the gcd it replaces.

The other half follows by division, and its power sums are *q<sub>k</sub> - p<sub>k</sub>*. If the division doesn't come out exact (too few terms to determine *A*), the same is tried with the other half, and finally with a gcd.

See [optimizations.md](optimizations.md) for the derivations, the measurements and the ideas that did not pay off.
