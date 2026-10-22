# 5. Root finding

`FindRoots(f, basis, field, odd_power_sums)` (`src/sketch_impl.h:2355`)
returns the roots of a monic f if f is square-free and splits into linear
factors over K, and the empty vector (failure) otherwise. It is Berlekamp's
trace algorithm with a shared Frobenius orbit, tower split parameters and,
for large f on b ∈ {16, 32, 64} CLMUL fields, a coset Montgomery remainder
tree with power-sum splits. Every variant returns the same root set as the
textbook algorithm (`poly_find_roots` in `tests/pyminisketch.py`); only the
order and the work differ.

The oracle `find_roots_fast` in `tests/pyminisketch.py` implements the same
ideas in about 80 lines and is the recommended companion to this chapter.

## The textbook algorithm

For f square-free and split, and any β ∈ K, T = Tr(βx) mod f takes the
value Tr(βα) ∈ {0, 1} at each root α, so gcd(f, T) is the product of the
x − α with Tr(βα) = 0 (C20). Recursing on gcd(f, T) and f/gcd(f, T) with
GF(2)-independent β_0, β_1, … separates all roots (C25). Computing T needs
the Frobenius orbit X_i = x^(2^i) mod f, i < b:
Tr(βx) mod f = Σ_{i<b} β^(2^i) X_i. One more squaring gives X_b, and
X_b = x iff f is square-free and splits (C22).

## Base cases

Degree 1, f = x + a: the root is a.

Degree 2, f = x^2 + a_1 x + a_0 (`src/sketch_impl.h:2097`):

```
if a_1 == 0: return FAIL            # x^2 + a_0 has a double root
c = a_0 · (1/a_1)^2                 # substitute x = a_1·y: y^2 + y = c
y = Qrt(c)
if y^2 + y != c: return FAIL        # no root in K (irreducible quadratic)
s = y · a_1
roots += [s, s + a_1]
```

## Split parameters

The split at depth k (the k-th trace attempt along a path, counting
attempts whose split was trivial) uses β_k. Three schemes, all GF(2)-bases
of K (C25, C27):

1. **Plain** (root degree < 8): β_k = basis·x^k, where x is the field
   generator (`Mul2`).
2. **GF(2^8) tower** (root degree ≥ 8 and 8 | b): with z = the generator of
   GF(2^8) below, β_(8l+j) = basis·x^l·z^j for l < b/8, j < 8 (here, as in
   all split parameters, x^l is the field element x = 2 raised to l, i.e.
   `Mul2` applied l times, not the polynomial variable of Tr(βx) mod f). Otherwise
   (root degree ≥ 8, 8 ∤ b) the plain scheme.
3. **GF(16) tower** (the Montgomery path: root degree ≥ 224, b ∈ {16, 32,
   64}, CLMUL): with z_4 = the generator of GF(16), β_(4l+j) =
   basis·x^l·z_4^j for l < b/4, j < 4. `MONT_TOWER_BITS` = 4 is the
   subfield degree s.

`SubfieldGenerator(sub)` (`src/sketch_impl.h:1296`), for sub a power of two
dividing b:

```
if sub == 1: return 1
y = 1
loop:
    y = Mul2(y)                                   # y = x, x^2, x^3, ...
    w = Σ_{l < b/sub} y^(2^(sub·l))               # relative trace into GF(2^sub)
    if w^(2^(sub/2)) != w: return w               # not in GF(2^(sub/2)): generates GF(2^sub)
```

Why the towers help: for β = basis·z^k with z in GF(2^s), (z^k)^(2^i)
depends only on i mod s, so the traces of the first s depths are
Tr(β_k x) = Σ_{r<s} (z^k)^(2^r)·Z_r with s bucket sums
Z_r = Σ_{i ≡ r mod s} basis^(2^i)·X_i (C27, C28).

## FindRoots

```
FindRoots(f, basis, Q_odd):              # f monic, d = deg f
    if f empty or basis == 0: return FAIL
    if d == 0: return []
    if d >= max(ORBIT_MIN_DEGREE = 8, 3):
        depths = min(b, 2 + floor(log2 b))     # ORBIT_EXTRA_DEPTHS = 2
        ok, traces, betas, tree, carry = ComputeTraceOrbit(f, basis, depths)
        if not ok: return FAIL                 # X_b != x: not square-free or not split (C22, C36)
        psum = none
        if Q_odd given and carry.tree (Montgomery path taken) and d >= 256:   # PSUM_SPLIT_MIN
            len = min(2·|Q_odd|, 2d + 4·ceil(sqrt d) + 64)
            psum[k-1] = Q_k for k = 1..len: Q_odd[k/2] for odd k, Q_(k/2)^2 for even k
        if not RecFindRoots(f, ff = true, depth = 0, randv = basis,
                            orbit = traces, orbit_degree = d, betas, carry, psum):
            return FAIL
    else:                                      # 1 <= d <= 7
        if not RecFindRoots(f, ff = false, depth = 0, randv = basis): return FAIL
    if number of roots != d: return FAIL
    return roots
```

`ff` ("fully factorizable") records that the node is known to be
square-free and split. Decode passes the sketch's odd syndromes as Q_odd:
they are the odd power sums of the roots f is expected to have.

## ComputeTraceOrbit (`src/sketch_impl.h:1782`)

Computes, from one pass over the orbit of x mod f, the split traces
traces[k] = Tr(β_k x) mod f for k < depths, the parameters β_0..β_(b−1),
and the full-split test.

```
ComputeTraceOrbit(f, basis, depths):        # d = deg f >= 3
    tower = (b mod 8 == 0) and depths <= 8  # always true when 8 | b
    betas = GF(2^8)-tower betas if tower else plain betas
    if tower and b mod 4 == 0 and d >= 224 and FFT-capable and cheap multiplier:
        res = MontgomeryOrbit(f, basis)     # below
        if res == NOT_SPLIT: return false
        if res == OK:
            betas = GF(16)-tower betas; tree.z = z_4
            return true, [trace0], betas, tree, carry   # deeper traces come from the carry
        # res == NOT_APPLICABLE: fall through, keep the GF(2^8) betas
    X = x                                    # X_0
    acc = zero accumulators: 8 buckets if tower, else one per depth (length d each)
    bcoef = basis; coef[k] = betas[k]
    for i in 0 .. b-1:
        if tower: acc[i mod 8] += bcoef · X;  bcoef = bcoef^2
        else:     acc[k] += coef[k] · X for k < depths;  coef[k] = coef[k]^2
        X = X^2 mod f                        # TraceMod.SquareAndReduce
    traces = CombineTowerTraces(acc) if tower else acc     # stripped
    return (X == x), traces, betas
```

`CombineTowerTraces`: traces[k] = Σ_{r<8} (z^k)^(2^r)·acc[r] for k < depths
(depths ≤ 8, so all within l = 0). CLMUL fields accumulate unreduced
(one `Wide` per coefficient per accumulator); cost (b + 8·depths)·d
products instead of b·depths·d.

Any node g of the tree divides f, so (h mod f) mod g = h mod g (C26): a node
at depth k < depths takes traces[k] mod g, one division, if
`ORBIT_NODE_FACTOR`·d_f ≤ deg g · b (that is 2·d_f ≤ deg g · b); smaller
nodes compute their own trace.

## Computing a trace at a node

### TraceMod (`src/sketch_impl.h:1030`)

A reducer for repeated squaring modulo a fixed monic g of degree d ≥ 3.
`SquareAndReduce(v)` maps v (reduced) to v^2 mod g; `trace(param)` returns
Tr(param·x) mod g:

```
trace(param):
    t = param·x
    repeat b - 1 times:
        t = t^2 mod g
        t += param·x          # coefficient 1 (the square's x coefficient is 0)
    return t                  # Σ_{i<b} (param x)^(2^i) mod g
```

Two modes, chosen by d against `TraceModTableCutoff` (704 for FFT-capable
cheap fields, 512 otherwise):

- **Square table** (d < cutoff). Precompute R_e = x^e mod g for even e
  from d (or d + 1 if d is odd) to 2d − 2, by repeated multiplication by
  x. Then v^2 mod g = Σ_{2i<d} v_i^2 x^(2i) + Σ_{2i≥d} v_i^2·R_(2i). CLMUL
  fields store the rows column-block-major (for each block of 8 output
  columns, all rows contiguous; the d mod 8 tail columns row-major) and
  accumulate 8 columns unreduced in registers, streaming the table once per
  squaring (about 2.1 cycles per product on Haswell, the PCLMULQDQ
  throughput limit). Generic fields add Multiplier(v_i^2)(R_(2i)) row by
  row.
- **Reciprocal (Barrett)** (d ≥ cutoff). inv = inv_series(rev(g), d) once.
  For v of length n ≥ d + 1 (after squaring), with k = n − d:
  q_rev = rev(v)[0:k]·inv mod x^k; q = rev(q_rev); v = (v + q·g mod x^d)
  truncated to d coefficients, stripped. If FFT-capable and 2d − 1 ≥ 128
  and t = ⌈log2 2d⌉ ≤ b, the 2^t-point transforms of inv and g are prepared
  once (4 transforms per reduction instead of 6) and used when the product
  would otherwise be an FFT of that size (la + lb − 1 > 2^(t−1)).

### SmallTrace (`src/sketch_impl.h:1913`)

For 3 ≤ d ≤ 8 (`SMALL_TRACE_MAX`) and only when ff: the same trace
iteration on fixed-size local arrays, with rows x^e mod g for every
e ∈ [d, 2d − 2] (lazy accumulation on CLMUL fields). Avoids the TraceMod
setup cost at tiny degrees.

### Trace from a carried relative trace

Below (MontTowerTrace, TowerTraceSmall).

## Coset Montgomery arithmetic

Used for roots of degree ≥ 224 (`MONT_ORBIT_MIN_DEGREE`) on FFT-capable
cheap fields, i.e. implementation 1 at b ∈ {16, 32, 64}. Notation from
chapter 3: C = C(m, B) a coset of 2^m = N points, Z_C = S^m + η_B its
vanishing polynomial, Z_C' = 1 (C29).

### MontCoset (`src/sketch_impl.h:1370`)

For a monic g with deg g ≤ N and no root on C, Z_C is invertible mod g, and
a residue r is held in **Montgomery form** r·Z_C mod g, as a polynomial of
degree < N represented by its N values on C (index j ↔ η_(B·2^m+j)).

```
Eval(g):                       # values of g on C, deg g <= N
    v = g[0:N] zero-padded to N
    if deg g == N:             # reduce mod Z_C: subtract the monic Z_C
        v[2^j] ^= 1 for j < m with (m AND j) == j
        v[0]   ^= η_B          # Z_C(0) = η_B = ⊕ β_k over bits k of B
    return NovelFFT(ToNovel(v), m, B)

SetModulus(g, values = Eval(g)):
    ga = values; inv_ga = BatchInvert(ga) or FAIL (g has a root on C)
    dga = values of g' on C    # g'_(i-1) = g_i for odd i, else 0

Derivative(v):                 # values of p' where p (deg < N) has values v
    return NovelFFT(NovelDerivative(NovelIFFT(v)))

Reduce(v, square):             # v := v / Z_C mod g, values in and out
    if not square: da = Derivative(v)
    q  = v · inv_ga            # pointwise: q = v / g mod Z_C
    dq = Derivative(q)
    v  = dq·ga + q·dga         # pointwise, one lazy reduction per point
    if not square: v += da

Square(v): v = v^2 pointwise; Reduce(v, square = true)    # Montgomery squaring

ToNormal(v, d): Reduce(v, false); return FromNovel(NovelIFFT(v))[0:d], stripped
ToCoefficients(v): return FromNovel(NovelIFFT(v)), stripped   # the Montgomery-form polynomial itself
```

Why `Reduce` is right (C30): with q ≡ a/g on C (deg q < N) and h = a + q·g,
Z_C divides h, r = h/Z_C ≡ a·Z_C^(−1) mod g, and differentiating
h = r·Z_C gives r(w) = h'(w) = a'(w) + q'(w)g(w) + q(w)g'(w) on C. For a
square a, a' = 0. Degrees: deg r < N if deg a < 2N; deg r < deg g if
deg a < N. A Montgomery squaring is therefore two size-N transforms plus
pointwise work, against about eight size-2N transforms for the Barrett
reducer. `mont_square` in `tests/pyminisketch.py` is the same computation.

### MontgomeryOrbit (`src/sketch_impl.h:1634`)

```
MontgomeryOrbit(f, basis):                  # d = deg f >= 224
    m = least m >= 1 with 2^m >= d;  N = 2^m
    if m + 2 > b: return NOT_APPLICABLE
    extend the Cantor basis to m + 2 elements
    for B in (1, 2, 3, 0):                  # the four cosets of W_m in W_(m+2)
        if SetModulus on C(m, B) succeeds: break
    else: return NOT_APPLICABLE             # f has a root on all four
    # Montgomery form of x: x·Z_C mod f
    c = x · Z_C = Σ_{j <= m, (m AND j) == j} x^(2^j + 1) + η_B·x
    r = values on C of (c mod f)            # c mod f has degree < d <= N
    r0 = r
    n_traces = min(4, ceil((m + CARRY_DEPTH_MARGIN) / 4))      # margin 2
    A = 0 (N values); T_l = 0 for l < n_traces       # unreduced accumulators
    bcoef = basis; tcoef[l] = basis·x^l for l < 4
    for i in 0 .. b-1:
        A += bcoef · r
        if i mod 4 == 0: T_l += tcoef[l] · r  for l < n_traces
        bcoef = bcoef^2; tcoef[l] = tcoef[l]^2 for all l < 4
        r = Square(r)                       # r = Montgomery form of X_(i+1)
    # full-split test (C36): X_b == x iff the reduced difference vanishes
    v = Reduce(r + r0, square = false)
    if any value of v != 0: return NOT_SPLIT
    carry = {coset (m, B), normal = false, t[l] = T_l (reduced) for l < n_traces,
             others empty, ga = values of f on C}
    trace0 = ToNormal(A, d)                 # Tr(basis·x) mod f, normal form
    return OK
```

(The code reduces c mod f with `FastDivMod(f, c, cutoff = 32)`.) The
accumulators hold the Montgomery forms of Tr(basis·x) mod f and of the
relative traces T_l = Tr_(K/GF(16))(basis·x^l·x) mod f
= Σ_{i ≡ 0 mod 4} (basis·x^l)^(2^i)·X_i. The test of `v` is exact: v has
degree < d ≤ N and vanishes on N points iff it is 0.

### Carrying relative traces (C28, C31)

`MontCarry` holds, for a node g, the relative traces T_l mod g for
l < `CARRY_TRACES` = 16/4 = 4 (empty = not carried), either as Montgomery
values on the node's coset (with g's values `ga` if known) or, for small
nodes, in normal form (`normal = true`; coefficient vectors, with [0]
standing for the zero residue so that empty keeps meaning "absent"). The
tree's `z` is z_4.

At depth k = 4l + j < `CARRY_DEPTHS` = 16 the split trace is
Tr(β_k x) = Tr_(GF(16)/GF(2))(z_4^j·T_l) = Σ_{r<4} (z_4^j T_l)^(2^r) mod g:
three modular squarings (C28).

```
MontTowerTrace(node coset, v = T_l values, deg g, zk, force_normal):
    y = acc = zk · v
    repeat 3 times: y = Square(y); acc += y
    if force_normal or deg g >= 512: return ToNormal(acc, deg g)       # MONT_NORMALIZE_MIN
    else: return ToCoefficients(acc)     # trace·Z_C mod g: a unit multiple (C32)

TowerTraceSmall(g, T, zk):               # normal form, schoolbook
    y = zk · T; out = y
    repeat 3 times: y = y^2 mod g; out += y
    return out (stripped)
```

The unnormalized form is valid because the trace is only used in a gcd with
g, and gcd(g, T·Z_C mod g) = gcd(g, T) when gcd(Z_C, g) = 1 (C32).

`MontCarryToFactor(parent, g, l_min, gv)` (`src/sketch_impl.h:1519`), for a
factor g of the node with 3 ≤ deg g, deg g > 16: copy the parent's coset and
traces l ≥ l_min (clear l < l_min), then halve the coset while m > 1 and
deg g ≤ 2^(m−1):

```
halve:                               # C = C0 ∪ C1, C0 = C(m-1, 2B), C1 = C(m-1, 2B+1)
    h = 2^(m-1)
    g0 = values of g on C0;  ga = values of g on C1   # from gv on the first halving
    inv_g0 = BatchInvert(g0)         # cannot fail: no root of f on C
    for each carried l with values v (length 2h; v[0:h] on C0, v[h:2h] on C1):
        q = NovelFFT_C1(NovelIFFT_C0(v[0:h] · inv_g0))    # q = v/g mod Z_C0, values on C1
        v = v[h:2h] + q · ga                              # pointwise
    (m, B) = (m-1, 2B+1)
```

Correct because Z_C = Z_C0·Z_C1 and Z_C0 ≡ 1 on C1 (C31): the child's form
(v + q·g)/Z_C0 has values v + q·g on C1. If no halving happens and gv is
known, ga = gv. (`mont_carry_to_factor` in `tests/pyminisketch.py`.)

### CosetDivide (`src/sketch_impl.h:1974`) (C33)

```
CosetDivide(coset, gv = values of g, deg g, h monic):
    if deg h > deg g: return FAIL
    hv = Eval(h); inv = BatchInvert(hv) or return FAIL
    qv = gv · inv                    # pointwise
    q = ToCoefficients(qv)
    if |q| == deg g - deg h + 1 and q's leading coefficient == 1: return q, hv, qv
    return FAIL
```

If g and h are monic, deg g ≤ N, and the interpolated q is monic of degree
deg g − deg h, then q·h = g (both sides agree on N points and have degree
≤ N with equal leading terms). So success certifies divisibility.

### PowerSumSplit (`src/sketch_impl.h:2008`) (C35)

For a Montgomery-form node g of degree n with normal-form split trace t
(t(α) ∈ {0, 1} at every root α) and power sums q[k−1] = Q_k = Σ α^k of its
roots, the roots A = {t = 1} have power sums
P_k = Σ_α t(α)α^k = Σ_{j<n} t_j Q_(j+k), and P_(2k) = P_k^2 since
t(α)^2 = t(α). So P is a syndrome sequence of A, and BM on its odd terms
returns Λ_A when there are enough terms.

```
PowerSumSplit(g, gv, t, q):              # n = deg g
    if |q| + 1 < 2n: return GCD_FALLBACK
    sq = ceil(sqrt n)
    two_m = min(|q| + 1 - n, n + 4·sq + 32) rounded down to even
    if two_m < n: return GCD_FALLBACK
    M = two_m / 2
    P_odd[i] = P_(2i+1) = Σ_{j<n} t_j Q_(j+2i+1)   for i < M
    for side in (1, 0):                  # t = 1 side first, then its complement
        osyn = P_odd if side == 1 else P_odd + Q_odd[0:M]   # complement: Q - P
        conn = BerlekampMasseyOdd(osyn, two_m, max_degree = M)
        if FAIL: continue
        if |conn| == 1: return TRIVIAL       # this side has no roots
        if |conn| == n + 1: return TRIVIAL   # all roots on this side
        gs = rev(conn)
        go = CosetDivide(node coset, gv, n, gs) or continue
        q1 = power sums of the t = 1 side: q1[2i] = P_odd[i], q1[k-1] = q1[k/2-1]^2 for even k <= two_m
        q0 = q[0:two_m] + q1
        return SPLIT: g1 = side-1 factor, g0 = side-0 factor, their coset values, q0, q1
    return GCD_FALLBACK
```

P_odd is computed as two half-length middle products split by the parity
of j: P_(2i+1) = Σ_a t_(2a) Q_(2(a+i)+1) + Σ_a t_(2a+1) Q_(2(a+i)+2), in
blocks of h outputs (h = ⌈n/2⌉ rounded up to even) by `LazyMiddleProduct`,
the outputs past the last full block as direct dot products.

Exactness (C35): with 2M ≥ n, one side has at most n/2 ≤ M roots; BM on a
side with |A| ≤ M returns Λ_A (C3, C5). BM's minimality (|B| ≤ |A| for its
result B) plus the divisibility certified by CosetDivide forces B = A, so a
SPLIT result is the true split and the children inherit correct power
sums. A BM result of length 0 means A = ∅ (C35c), and the complement has
power sums Q_k − P_k (C35d). The 4·sq + 32 slack gives the children
(which need twice their degree) enough sums unless the split is
unbalanced by more than about four standard deviations; FindRoots' root
slack is 4·sq + 64.

## RecFindRoots (`src/sketch_impl.h:2084`)

```
RecFindRoots(g, ff, depth, randv, orbit, orbit_degree, betas, carry, psum):
    d = deg g                                       # g monic, d >= 1
    if d <= 2: base case (above)
    use_psum = psum and carry.tree and not carry.normal and d >= 256 and |psum| + 1 >= 2d
    loop (iter = 0, 1, ...):
        if betas:
            if depth >= |betas|: return FAIL
            randv = betas[depth]
        at_root = orbit and depth < |orbit| and d == orbit_degree
        # 1. the split trace t ≡ Tr(randv·x) mod g (or a unit multiple)
        if carry.tree and not at_root and depth < 16 and carry.t[depth/4] present:
            zk = z_4^(depth mod 4)
            t = TowerTraceSmall(g, carry.t[depth/4], zk)             if carry.normal
                MontTowerTrace(node coset, carry.t[depth/4], d, zk,
                               force_normal = use_psum)              otherwise
        elif orbit and depth < |orbit| and 2·orbit_degree <= d·b:
            t = orbit[depth] mod g                                    # FastDivMod
        elif ff and d <= 8:
            t = SmallTrace(g, randv)
        else:
            t = TraceMod(g).trace(randv)          # reducer built once per node
        # 2. fail-fast test, only on the small-root path (C21)
        if iter >= 1 and not ff:
            if (t^2 mod g) + t != 0: return FAIL
            ff = true
        # 3. depth bound (C24)
        if ff and b - depth < 64 and d - 1 >= 2^(b - depth): return FAIL
        depth += 1; randv = Mul2(randv)
        # 4. split
        if use_psum:
            r = PowerSumSplit(g, values of g on the node coset, t, psum)
            if r == SPLIT: (g0, g1) = its factors; break
            if r == TRIVIAL: continue
        h = FastGCD(t, g)
        if 1 <= deg h < d: g0 = h; break
    # 5. factors: g0 holds the roots with trace 0, g1 the rest
    if not psum split:
        g0 = monic(g0)
        if carry.tree and not carry.normal and d >= 64:          # COSET_DIV_MIN
            g1 = CosetDivide(node coset, values of g, d, g0)     # always succeeds here
        else:
            g1 = g / g0                                          # FastDivMod
    # 6. carry the relative traces still needed (depth is the children's depth)
    if carry.tree and depth < 16:
        l_min = depth / 4
        for (factor, its coset values if known) in (g0, g1):
            if deg factor < 3: no carry
            elif not carry.normal and deg factor > 16:            # SMALL_CARRY_MAX
                child carry = MontCarryToFactor(carry, factor, l_min, values)
            else:
                child carry = normal form: for l >= l_min with carry.t[l] present,
                    (carry.t[l] if carry.normal else ToNormal(carry.t[l], d), once per node)
                    mod factor
    release this node's coset state
    # 7. recurse: the t = 1 side first
    if not RecFindRoots(g1, ff, depth, randv, ..., carry of g1, q1 if psum split): return FAIL
    return RecFindRoots(g0, true, depth, randv, ..., carry of g0, q0 if psum split)
```

Notes:

- On the Montgomery path the root's orbit has one entry, so the root uses
  trace0 at depth 0 and the carry from depth 1 on; a trivial split at the
  root moves to the carried trace of depth 1.
- On the GF(2^8)-tower path nodes use the orbit traces at depths
  < depths when 2·d_f ≤ d·b, and compute their own trace otherwise.
- A Montgomery-form split trace of a node below degree 512 is the
  unnormalized ToCoefficients form, whose degree can reach N − 1 ≥ d; the
  gcd handles that.
- When use_psum is false at a node with carried psum (degree below 256,
  or too few sums), the node uses the gcd; its children get no power sums.

### The second recursion never fails (C23)

g0 = gcd(g, t) collects the roots with t = 0. At a root α ∉ K of an
irreducible factor of g, t + t^2 = randv·(α + α^(2^b)) ≠ 0 (C21), so
t(α) ∉ {0, 1}: g0 has only roots in K, all simple (t has derivative randv at
each root). So g0 is square-free and split, and its recursion is called with
ff = true and cannot fail. Irreducible factors always end up in g1, which
inherits ff.

### Depth bound (C24)

After k attempts with independent β_0..β_(k−1) the roots of a node share
Tr(β_i α) for i < k, so they lie in an affine subspace of dimension b − k
and number at most 2^(b−k). A split node with d > 2^(b−depth) is therefore
impossible; the check `(d − 1) >> (b − depth) == 0` (skipped when
b − depth ≥ 64) guards against it and also bounds the loop. A trivial split
at degree d has probability about 2^(1−d) for random parameters (C41).

## Fallback edges and failure returns

| where | condition | effect |
|---|---|---|
| FindRoots | f empty or basis = 0 | FAIL (unreachable from decode) |
| FindRoots | d < 8 | small path: plain betas, no orbit, ff = false |
| ComputeTraceOrbit | not (8 ∣ b, d ≥ 224, b ∈ {16,32,64}, CLMUL) | standard orbit with TraceMod |
| MontgomeryOrbit | m + 2 > b (only b = 16, d > 2^14) | NOT_APPLICABLE → standard orbit, GF(2^8) betas |
| MontgomeryOrbit | f has a root on all four cosets | NOT_APPLICABLE → standard orbit, GF(2^8) betas (common at b = 16: each coset of 2^m ≥ d points expects d·2^m/2^16 roots, so at d = 300 all four are hit about two times in three; rare at b = 32, 64) |
| ComputeTraceOrbit / MontgomeryOrbit | X_b ≠ x | FAIL (not square-free or not split) |
| FindRoots | no Q_odd, or d < 256, or no Montgomery carry | no power-sum splits |
| RecFindRoots | degree 2 with a_1 = 0 | FAIL |
| RecFindRoots | degree 2, y^2 + y ≠ c | FAIL |
| RecFindRoots | depth ≥ b with betas | FAIL |
| RecFindRoots | fail-fast test nonzero (small path, iter ≥ 1) | FAIL |
| RecFindRoots | depth bound violated | FAIL |
| PowerSumSplit | too few sums, both sides fail BM or CosetDivide | gcd |
| PowerSumSplit | one side empty or full | next depth (no gcd) |
| RecFindRoots | gcd trivial (degree 0 or d) | next depth |
| RecFindRoots | first (g1) recursion fails | FAIL |
| FindRoots | root count ≠ d | FAIL (defensive) |

On the orbit path (d ≥ 8) the root is already known to be square-free and
split, so after the orbit test no FAIL edge is reachable for a correct
implementation (the remaining ones are `CHECK_RETURN` guards); all real
failures on that path come from the orbit test. On the small path
(3 ≤ d ≤ 7) failures come from the degree-2 base case, the fail-fast test
and the first recursion.

## Output order

Roots are appended depth-first: at each node the t = 1 factor's roots
before the t = 0 factor's, and at a quadratic leaf s before s + a_1. The
order is a deterministic function of the basis (seed, through `FromSeed`,
whose key depends on the implementation), the split-parameter scheme
(plain, GF(2^8) or GF(16): root degree, b, implementation and whether the
Montgomery path ran), trivial splits moving on to the next parameter, the
quadratic-leaf rule, and the field's internal representation (b = 58,
implementation 1). It does not depend on the Cantor basis, the coset, the
carried traces, power-sum splits or CosetDivide: those change how a split
is computed, not which split it is. (A clean-room port without them
reproduced libminisketch's order exactly.)
