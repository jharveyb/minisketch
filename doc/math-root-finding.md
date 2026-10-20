# The mathematics of root finding for large sketches

This continues [math.md](math.md), in the same spirit, for the root-finding
techniques used when a sketch holds many elements. Costs and measurements
are in [optimizations.md](optimizations.md).

## Root finding for large sketches

**The cost of a trace**

The roots of *L* are found with the Berlekamp trace algorithm. The trace *Tr(y) = y + y<sup>2</sup> + y<sup>4</sup> + ... + y<sup>2<sup>b-1</sup></sup>* only takes the values 0 and 1, and is linear. So for a constant *&beta;*, *gcd(f, Tr(&beta;x) mod f)* collects exactly the roots *m* with *Tr(&beta;m) = 0*, about half of them. Repeating this with different *&beta;* splits *f* into degree-1 pieces.

Computing *Tr(&beta;x) mod f* needs *x, x<sup>2</sup>, x<sup>4</sup>, ..., x<sup>2<sup>b-1</sup></sup> mod f*: *b - 1* squarings modulo a polynomial of degree *n*. For large *n* these squarings, the gcds, and the divisions between pieces dominate decoding.

**Squaring by evaluation**

A polynomial of degree less than *N* is determined by its values at *N* points, and products are pointwise in that form. The FFT used here evaluates on an affine subspace *C* of the field (*N = 2<sup>k</sup>* elements closed under XOR, shifted by a constant). The polynomial *Z(x)*, the product of *(x - w)* over all *w &isin; C*, has a special shape: *x<sup>N</sup> + ... + x + c*, with only power-of-two exponents.

We keep every residue *r* in the scaled form *rZ mod f*, as is done for integers in Montgomery multiplication. To reduce a product *a*:

* Find *q* with *a + qf* divisible by *Z*. Since *Z* vanishes on *C*, that just means *q(w) = a(w)/f(w)* at every point *w* of *C*.
* The reduced value is *(a + qf)/Z*. That is 0/0 on *C*, but differentiating *a + qf = rZ* gives *(a + qf)' = r'Z + rZ'*, so at the points of *C*, *r(w) = (a + qf)'(w) / Z'(w)*.

In fields where *x + x = 0*, the derivative of every *x<sup>2<sup>j</sup></sup>* with *j &ge; 1* is zero. So *Z' = 1*, and the derivative of a square is zero too. For squaring this leaves:

* *r(w) = q'(w)f(w) + q(w)f'(w)*

Only *q'* needs a transform: an inverse FFT to get *q*, a differentiation, and an FFT back. A modular squaring costs two transforms of the size of *f*.

**Traces for every node from one partial trace**

The split pieces all divide *f*, so a residue mod *f* is also one mod each piece. Write the trace in two steps through the subfield *GF(16)*:

* *T = Tr<sub>GF(2<sup>b</sup>)/GF(16)</sub>(&beta;x)*, whose values at the roots lie in *GF(16)*.
* *Tr(&beta;z<sup>k</sup>x) = Tr<sub>GF(16)/GF(2)</sub>(z<sup>k</sup>T) = z<sup>k</sup>T + (z<sup>k</sup>T)<sup>2</sup> + (z<sup>k</sup>T)<sup>4</sup> + (z<sup>k</sup>T)<sup>8</sup>*, for *z* generating *GF(16)*.

The splits of four consecutive levels therefore need only three squarings each of *T* reduced mod the piece. *T* is passed from each piece to its factors, still in evaluated form. When the pieces get small enough, the evaluation set is halved: *C* splits into halves *C<sub>0</sub>* and *C<sub>1</sub>*, *Z<sub>C</sub> = Z<sub>C<sub>0</sub></sub>Z<sub>C<sub>1</sub></sub>*, and *Z<sub>C<sub>0</sub></sub> = 1* on *C<sub>1</sub>*, so the change of representation is one more reduction step. Once a piece *g* is split into *g<sub>1</sub>g<sub>2</sub>*, *g<sub>2</sub>* follows from *g<sub>2</sub>(w) = g(w)/g<sub>1</sub>(w)* without any long division.

**Sketches of halves**

Let *t = Tr(&beta;x) mod g*, and let *A* be the roots of *g* where *t = 1*. The power sums of *A* are:

* *p<sub>k</sub> = &sum;<sub>m&isin;A</sub> m<sup>k</sup> = &sum;<sub>m</sub> t(m)m<sup>k</sup> = &sum;<sub>j</sub> t<sub>j</sub>q<sub>j+k</sub>*

where *q<sub>k</sub>* are the power sums of all roots of *g*. For *f* itself, those are the sketch's own coefficients. Because *t(m)* is 0 or 1, *t(m)<sup>2</sup> = t(m)*, and so *p<sub>2k</sub> = p<sub>k</sub><sup>2</sup>*, just like a real sketch. The list *p* is therefore a sketch of *A*, computed from the sketch of the whole set, and Berlekamp-Massey decodes it directly into the polynomial of *A*. That's a problem of about half the size, cheaper than the gcd it replaces.

The other half follows by division, and its power sums are *q<sub>k</sub> - p<sub>k</sub>*. If the division doesn't come out exact (too few terms to determine *A*), the same is tried with the other half, and finally with a gcd.
