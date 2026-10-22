#!/usr/bin/env python3
# Copyright (c) 2026 The minisketch developers
# Distributed under the MIT software license, see the accompanying
# file LICENSE or http://www.opensource.org/licenses/mit-license.php.

"""Conformance test: libminisketch's internal decode stages vs. the oracle.

Runs seeded random cases through src/oracle_dump.cpp (the C++ stages:
BerlekampMasseyOdd, the textbook BerlekampMasseyRef, the orbit's full-split
test and traces, FindRoots with and without power sums, and the full decode)
and compares every stage against the optimized-decode oracle in
pyminisketch.py (bm_odd, bm_reference, frobenius_orbit / trace_from_orbit,
find_roots_fast, decode_fast) and the textbook Minisketch.decode.

Unlike differential.py (public API, capacity <= 16), the case classes here
reach the high-degree tiers: the Karatsuba (L >= 64) and FFT (L >= 512, or
>= 40 on table fields) BM tail checks, the Montgomery orbit (degree >= 224,
with roots planted off the first candidate coset to force it), power-sum
splits (>= 256) and CosetDivide (>= 64). oracle_dump reports MINISKETCH_STAT
counters per stage, and the run fails if a tier a class is meant to reach
never ran.

Usage: conformance.py --dump path/to/oracle-dump [--iters N] [--big] [-v]

Reproduce a failure by setting MINISKETCH_TEST_SEED to the printed seed.
"""

import argparse
import collections
import json
import os
import random
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyminisketch import (  # noqa: E402
    Minisketch, bm_odd, bm_reference, cantor_basis, coset_points, decode_fast,
    elem_trace, find_roots_fast, frobenius_orbit, gf2_rank, make_gf, poly_eval,
    poly_from_roots, poly_mul, poly_trim, power_sums, trace_from_orbit,
)

GENERIC, CLMUL, CLMUL_TRI = 0, 1, 2
IMPL_NAMES = {GENERIC: 'generic', CLMUL: 'clmul', CLMUL_TRI: 'clmultri'}
SMALL_FIELDS = list(range(2, 17)) + [32]


class Mismatch(Exception):
    pass


class Dump:
    """Batches case lines through one oracle-dump process."""

    def __init__(self, path):
        self.path = path

    def run(self, lines):
        if not lines:
            return []
        proc = subprocess.run([self.path], input='\n'.join(lines) + '\n',
                              capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            raise Mismatch(f"oracle-dump exited with {proc.returncode}: {proc.stderr.strip()}")
        out = [json.loads(line) for line in proc.stdout.splitlines()]
        assert len(out) == len(lines)
        return out


def line(op, bits, impl, param, a, b=None):
    parts = [op, bits, impl, param, len(a)] + list(a)
    if b is not None:
        parts += [len(b)] + list(b)
    return ' '.join(str(x) for x in parts)


def random_set(rng, bits, n, avoid=()):
    n = min(n, (1 << bits) - 1 - len(avoid))
    out = set()
    while len(out) < n:
        v = rng.randrange(1, 1 << bits)
        if v not in avoid:
            out.add(v)
    return sorted(out)


def odd_power_sums(roots, count, gf):
    """[Q_1, Q_3, ..., Q_(2 count - 1)]: the odd syndromes of a sketch."""
    out = [0] * count
    for r in roots:
        sq, p = gf.sqr(r), r
        for i in range(count):
            out[i] ^= p
            p = gf.mul(p, sq)
    return out


def all_power_sums(odd, gf):
    """[0, Q_1, ..., Q_(2 len(odd))] from the odd ones (Q_2k = Q_k^2)."""
    q = [0] * (2 * len(odd) + 1)
    for i, v in enumerate(odd):
        q[2 * i + 1] = v
    for k in range(2, 2 * len(odd) + 1, 2):
        q[k] = gf.sqr(q[k // 2])
    return q


def irreducible_quadratic(rng, gf):
    """x^2 + x + c with Tr(c) = 1 (no roots in the field)."""
    while True:
        c = rng.randrange(1, 1 << gf.field_size)
        if elem_trace(c, gf) == 1:
            return [c, 1, 1]


def expect_eq(what, got, want, ctx):
    if got != want:
        raise Mismatch(f"{what}: C++ {short(got)} != {short(want)}\n  case: {short(ctx)}")


def short(v, limit=400):
    s = json.dumps(v) if not isinstance(v, str) else v
    return s if len(s) <= limit else s[:limit] + f'... ({len(s)} chars)'


class Conformance:
    def __init__(self, dump, rng, verbose=False, big=False):
        self.dump = dump
        self.rng = rng
        self.verbose = verbose
        self.big = big
        self.gfs = {}
        self.stats = collections.defaultdict(collections.Counter)  # class -> stat -> count
        self.cases = collections.Counter()
        self.impls = self.available_impls()
        self.cxx_cantor = {}

    def gf(self, bits):
        if bits not in self.gfs:
            self.gfs[bits] = make_gf(bits)
        return self.gfs[bits]

    def available_impls(self):
        probes = [(b, i) for b in SMALL_FIELDS for i in IMPL_NAMES]
        res = self.dump.run([line('cantor', b, i, 0, []) for b, i in probes])
        impls = collections.defaultdict(list)
        for (b, i), r in zip(probes, res):
            if 'error' not in r:
                impls[b].append(i)
        return impls

    def record(self, cls, result):
        self.cases[cls] += 1
        for key, val in result.items():
            if key.startswith('stats_'):
                for name, n in val.items():
                    self.stats[cls][f"{key[6:]}:{name}"] += n

    # --- Cantor basis -------------------------------------------------------

    def check_cantor(self):
        """The additive FFT's basis; also needed to plant roots off the first
        Montgomery coset. pyminisketch's choice of Artin-Schreier root may
        legitimately differ (C++ uses field.Qrt), so only the defining
        recurrence is required; equality is reported."""
        for bits in (16, 32):
            gf = self.gf(bits)
            py = cantor_basis(gf)
            for impl in self.impls[bits]:
                basis = self.dump.run([line('cantor', bits, impl, 0, [])])[0]['basis']
                if basis[0] != 1 or gf2_rank(basis) != bits or any(
                        gf.sqr(basis[j]) ^ basis[j] != basis[j - 1] for j in range(1, bits)):
                    raise Mismatch(f"cantor basis {bits}/{IMPL_NAMES[impl]} is not a Cantor basis: {basis}")
                self.cxx_cantor[bits, impl] = basis
                note = "equal to" if basis == py else "differs from"
                if self.verbose or basis != py:
                    print(f"  cantor {bits}/{IMPL_NAMES[impl]}: {note} pyminisketch.cantor_basis")

    def avoid_first_coset(self, bits, impl, d):
        """Points of the first coset MontgomeryOrbit tries (block 1) for a
        degree-d polynomial: roots planted off it force the Montgomery orbit."""
        basis = self.cxx_cantor[bits, impl]
        m = max(1, (d - 1).bit_length())
        return set(coset_points(basis, basis[m], m))

    # --- decode cases -------------------------------------------------------

    def decode_case(self, bits, impl, odd, max_count, planted=None, textbook=True, mont=True):
        """One decode line, and the checks on its result."""
        seed = self.rng.randrange(1 << 64) if self.rng.randrange(4) else (1 << 64) - 1
        return (line('decode', bits, impl, max_count, odd, [seed]),
                dict(kind='decode', bits=bits, impl=impl, odd=odd, max_count=max_count,
                     planted=planted, textbook=textbook, mont=mont))

    def check_decode(self, case, res):
        gf = self.gf(case['bits'])
        odd, max_count = case['odd'], case['max_count']
        ctx = dict(bits=case['bits'], impl=IMPL_NAMES[case['impl']], max_count=max_count,
                   capacity=len(odd), odd=odd if len(odd) <= 40 else f"<{len(odd)} syndromes>")
        want_bm = bm_odd(odd, max_count, gf)
        expect_eq("BerlekampMasseyOdd vs bm_odd", res['bm'], want_bm, ctx)
        expect_eq("BerlekampMasseyRef vs bm_reference", res['bm_ref'], bm_reference(odd, max_count, gf), ctx)
        expect_eq("bm_odd vs bm_reference (oracle self-check)", want_bm, bm_reference(odd, max_count, gf), ctx)
        want = decode_fast(odd, max_count, gf, self.rng.randrange(1, 1 << gf.field_size), case['mont'])
        expect_eq("decode vs decode_fast", res['decode'], want, ctx)
        if case['textbook']:
            sketch = Minisketch(case['bits'], len(odd))
            sketch.odd_syndromes = list(odd)
            expect_eq("decode vs Minisketch.decode", res['decode'], sketch.decode(max_count), ctx)
        planted = case['planted']
        if planted is not None and len(planted) <= max_count:
            expect_eq("decode vs planted set", res['decode'], sorted(planted), ctx)

    # --- roots cases --------------------------------------------------------

    def roots_case(self, bits, impl, roots, extra=None, psum_count=None, perturb=False,
                   oracle=True, mont=True):
        """FindRoots on prod (x - r) (times extra, if given: then not split)."""
        gf = self.gf(bits)
        f = poly_from_roots(roots, gf)
        if extra is not None:
            f = poly_mul(f, extra, gf)
        odd = None
        if psum_count:
            odd = odd_power_sums(roots, psum_count, gf)
            if perturb:
                odd[self.rng.randrange(len(odd))] ^= self.rng.randrange(1, 1 << bits)
        basis = self.rng.randrange(1, 1 << bits)
        return (line('roots', bits, impl, basis, f, odd),
                dict(kind='roots', bits=bits, impl=impl, roots=roots, f=f, split=extra is None,
                     odd=odd, oracle=oracle, mont=mont, basis=basis))

    def check_roots(self, case, res):
        bits = case['bits']
        gf = self.gf(bits)
        f, roots = case['f'], case['roots']
        d = len(f) - 1
        ctx = dict(bits=bits, impl=IMPL_NAMES[case['impl']], degree=d, basis=case['basis'],
                   f=f if d <= 16 else f"<degree {d}>")
        want = sorted(roots) if case['split'] else None
        if d >= 3:
            # The full-split test (C22): X_b = x mod f.
            if d <= 64:
                expect_eq("full-split flag vs frobenius_orbit", res['split'],
                          frobenius_orbit(f, gf)[bits] == [0, 1], ctx)
            expect_eq("full-split flag vs construction", res['split'], case['split'], ctx)
            if res['split']:
                betas = res['betas']
                if gf2_rank(betas) != bits:
                    raise Mismatch(f"split parameters are not a basis: {betas}\n  case: {short(ctx)}")
                orbit = frobenius_orbit(f, gf) if d <= 128 else None
                for k, t in enumerate(res['traces']):
                    t = poly_trim(t)
                    if orbit is not None:
                        expect_eq(f"trace {k} vs trace_from_orbit", t,
                                  trace_from_orbit(orbit, betas[k], f, gf), ctx)
                    else:
                        # deg t < d and the d roots are distinct: its values
                        # there, Tr(beta_k r), determine it.
                        for r in roots:
                            if poly_eval(t, r, gf) != elem_trace(gf.mul(betas[k], r), gf):
                                raise Mismatch(f"trace {k} wrong at root {r}\n  case: {short(ctx)}")
        expect_eq("FindRoots vs construction", res['roots'], want, ctx)
        if case['odd'] is not None:
            expect_eq("FindRoots with power sums vs construction", res['roots_psum'], want, ctx)
        if case['oracle']:
            q = all_power_sums(case['odd'], gf) if case['odd'] is not None else None
            got = find_roots_fast(f, gf, case['basis'], q, case['mont'])
            expect_eq("FindRoots vs find_roots_fast", res['roots'], got, ctx)

    # --- running ------------------------------------------------------------

    def run_class(self, cls, cases, require=()):
        """Run (line, case) pairs; require: stats ('stage:name') that must
        have fired at least once in this class."""
        t0 = time.time()
        results = self.dump.run([c[0] for c in cases])
        t1 = time.time()
        for (_, case), res in zip(cases, results):
            if 'error' in res:
                raise Mismatch(f"oracle-dump error {res['error']} on {short(case)}")
            if case['kind'] == 'decode':
                self.check_decode(case, res)
            else:
                self.check_roots(case, res)
            self.record(cls, res)
        for stat in require:
            if not self.stats[cls][stat]:
                raise Mismatch(f"class {cls}: tier {stat} never ran (stats: {dict(self.stats[cls])})")
        print(f"  {cls}: {len(cases)} cases OK (C++ {t1 - t0:.1f}s, oracle {time.time() - t1:.1f}s)", flush=True)
        if self.verbose:
            for name, n in sorted(self.stats[cls].items()):
                print(f"      {name}: {n}")

    def small_decode(self, iters):
        """All fields and implementations, capacity < 40: sparse sets (early
        termination), perturbed syndromes, garbage, over-capacity and exact
        sets, random max_count."""
        cases = []
        for _ in range(iters):
            for bits in SMALL_FIELDS:
                gf = self.gf(bits)
                for impl in self.impls[bits]:
                    cap = self.rng.randrange(1, 40 if bits <= 16 else 12)
                    kind = self.rng.randrange(5)
                    planted = None
                    if kind in (0, 1):  # sparse, maybe perturbed
                        planted = random_set(self.rng, bits, self.rng.randrange(cap // 3 + 1))
                        odd = odd_power_sums(planted, cap, gf)
                        if kind == 1:
                            odd[self.rng.randrange(cap)] ^= self.rng.randrange(1, 1 << bits)
                            planted = None
                    elif kind == 2:  # garbage
                        odd = [self.rng.randrange(1 << bits) for _ in range(cap)]
                    else:  # exact or over capacity
                        n = self.rng.randrange(cap + 1) if kind == 3 else self.rng.randrange(cap + 1, 2 * cap + 3)
                        planted = random_set(self.rng, bits, n)
                        odd = odd_power_sums(planted, cap, gf)
                    max_count = cap if self.rng.randrange(2) else self.rng.randrange(cap + 1)
                    cases.append(self.decode_case(bits, impl, odd, max_count, planted))
        self.run_class('decode-small', cases, require=['bm:bm-early-stop'])

    def small_roots(self, iters):
        """FindRoots on planted root sets of degree < 40, split and not split
        (an irreducible quadratic factor, or a double root at degree >= 8,
        where the orbit test applies), with true and perturbed power sums."""
        cases = []
        for _ in range(iters):
            for bits in SMALL_FIELDS:
                gf = self.gf(bits)
                for impl in self.impls[bits]:
                    n = self.rng.randrange(min(40 if bits <= 16 else 24, (1 << bits) - 1) + 1)
                    roots = random_set(self.rng, bits, n)
                    kind = self.rng.randrange(4)
                    extra = None
                    if kind == 1:
                        extra = irreducible_quadratic(self.rng, gf)
                    elif kind == 2 and len(roots) >= 8:
                        extra = [roots[0], 1]
                    psum = len(roots) + self.rng.randrange(1, 20) if kind != 1 else None
                    cases.append(self.roots_case(bits, impl, roots, extra, psum, perturb=kind == 3))
        self.run_class('roots-small', cases)

    def large_16(self, iters):
        """GF(2^16), degrees ~64, 224-300 and ~1100."""
        bits, gf = 16, self.gf(16)
        impls = [i for i in (GENERIC, CLMUL) if i in self.impls[bits]]
        for _ in range(iters):
            # BM tail-check tiers: Karatsuba (lazy, 64 <= L < 512), FFT
            # (table fields from L >= 40, lazy from 512), with true sketches,
            # one perturbed tail syndrome (the check rejects, BM resumes) and
            # over-capacity sets.
            cases = []
            for impl in impls:
                for n in (self.rng.randrange(64, 100), self.rng.randrange(224, 300)):
                    cap = n + self.rng.randrange(40, 80)
                    planted = random_set(self.rng, bits, n)
                    odd = odd_power_sums(planted, cap, gf)
                    cases.append(self.decode_case(bits, impl, odd, cap, planted, textbook=n < 100))
                    bad = list(odd)
                    bad[self.rng.randrange(n + 2, cap)] ^= self.rng.randrange(1, 1 << bits)
                    cases.append(self.decode_case(bits, impl, bad, cap, None, textbook=False))
                planted = random_set(self.rng, bits, self.rng.randrange(70, 90))
                cases.append(self.decode_case(bits, impl, odd_power_sums(planted, 64, gf), 64, planted, textbook=False))
            req = ['bm:bm-verify-kara', 'bm:bm-verify-fft', 'bm:bm-verify-rejected'] if CLMUL in impls else []
            self.run_class('decode-16-tiers', cases, require=req)

            # Root finding at degree 224-300: random roots (the Montgomery
            # orbit usually finds a root on every candidate coset and falls
            # back), roots off the first coset (Montgomery orbit, power-sum
            # splits, CosetDivide), and not split (the Montgomery orbit's
            # full-split test fails).
            cases = []
            for impl in impls:
                d = self.rng.randrange(224, 300)
                psum = d + 60
                cases.append(self.roots_case(bits, impl, random_set(self.rng, bits, d), psum_count=psum))
                n = self.rng.randrange(256, 300)
                avoid = self.avoid_first_coset(bits, impl, n) if impl == CLMUL else ()
                roots = random_set(self.rng, bits, n, avoid)
                cases.append(self.roots_case(bits, impl, roots, psum_count=len(roots) + 60))
                cases.append(self.roots_case(bits, impl, roots, psum_count=len(roots) + 60, perturb=True))
                roots = random_set(self.rng, bits, d, self.avoid_first_coset(bits, impl, d + 2) if impl == CLMUL else ())
                cases.append(self.roots_case(bits, impl, roots, irreducible_quadratic(self.rng, gf)))
                cases.append(self.roots_case(bits, impl, random_set(self.rng, bits, self.rng.randrange(64, 100)),
                                             psum_count=120))
            req = ['orbit:mont-orbit', 'orbit:mont-orbit-not-split', 'roots:mont-orbit',
                   'roots_psum:psum-split', 'roots:coset-divide', 'roots_psum:coset-divide'] if CLMUL in impls else []
            self.run_class('roots-16-montgomery', cases, require=req)

            # Degree ~1100 (clmul: FFT tail check, Montgomery orbit at
            # m = 11, power-sum splits down the tree): one decode and one
            # FindRoots with power sums, roots off the first coset.
            cases = []
            if CLMUL in impls:
                n = self.rng.randrange(1050, 1150)
                avoid = self.avoid_first_coset(bits, CLMUL, n)
                planted = random_set(self.rng, bits, n, avoid)
                cap = n + 40
                cases.append(self.decode_case(bits, CLMUL, odd_power_sums(planted, cap, gf), cap, planted,
                                              textbook=False))
                self.run_class('decode-16-1100', cases,
                               require=['bm:bm-verify-fft', 'bm:bm-early-stop', 'decode:mont-orbit',
                                        'decode:psum-split', 'decode:coset-divide'])

    def large_other(self, iters):
        """The remaining BM tail-check tier: lazily-reducing fields whose size
        is not a power of two have no FFT, so L >= 512 uses the generic
        polynomial product (GF(2^12), clmul)."""
        bits, gf = 12, self.gf(12)
        if CLMUL not in self.impls[bits]:
            return
        cases = []
        for _ in range(iters):
            n = self.rng.randrange(512, 560)
            planted = random_set(self.rng, bits, n)
            cap = n + 40
            cases.append(self.decode_case(bits, CLMUL, odd_power_sums(planted, cap, gf), cap, planted,
                                          textbook=False))
        self.run_class('decode-12-polymul', cases, require=['bm:bm-verify-polymul'])

    def large_32(self, iters):
        """GF(2^32): BM tiers at L ~ 64-100 (Karatsuba on clmul, FFT on the
        table field), and root finding at degree ~230 / ~260 (Montgomery
        orbit; power-sum splits from 256). The pure-Python root finder is
        slow here (GF2Ops), so the high-degree cases are checked against the
        planted roots, and against find_roots_fast only with --big."""
        bits, gf = 32, self.gf(32)
        impls = [i for i in (GENERIC, CLMUL) if i in self.impls[bits]]
        for _ in range(iters):
            cases = []
            for impl in impls:
                n = self.rng.randrange(64, 80)
                cap = n + 40
                planted = random_set(self.rng, bits, n)
                odd = odd_power_sums(planted, cap, gf)
                cases.append(self.decode_case(bits, impl, odd, cap, planted, textbook=False, mont=False))
                bad = list(odd)
                bad[self.rng.randrange(n + 2, cap)] ^= self.rng.randrange(1, 1 << bits)
                cases.append(self.decode_case(bits, impl, bad, cap, None, textbook=False, mont=False))
            req = ['bm:bm-verify-kara'] if CLMUL in impls else []
            self.run_class('decode-32-tiers', cases, require=req + ['bm:bm-verify-fft'])

            cases = []
            for impl in impls:
                for d in (self.rng.randrange(224, 240), self.rng.randrange(256, 270)):
                    roots = random_set(self.rng, bits, d)
                    cases.append(self.roots_case(bits, impl, roots, psum_count=d + 50,
                                                 oracle=self.big, mont=False))
            req = ['roots:mont-orbit', 'roots_psum:psum-split'] if CLMUL in impls else []
            self.run_class('roots-32-montgomery', cases, require=req)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dump', required=True, help="path to the oracle-dump executable")
    parser.add_argument('--iters', type=int, default=1, help="rounds of every case class (default 1)")
    parser.add_argument('--big', action='store_true', help="also run the slow pure-Python root finder on the 32-bit high-degree cases")
    parser.add_argument('-v', '--verbose', action='store_true', help="print per-class stage counters")
    args = parser.parse_args()

    seed = os.environ.get('MINISKETCH_TEST_SEED')
    seed = int(seed, 0) if seed else random.randrange(1 << 63)
    print(f"Test seed: {seed:#x} (reproduce with MINISKETCH_TEST_SEED={seed:#x})")
    rng = random.Random(seed)
    random.seed(rng.randrange(1 << 63))  # pyminisketch's own randomness

    t0 = time.time()
    conf = Conformance(Dump(args.dump), rng, args.verbose, args.big)
    try:
        print("implementations: " + ', '.join(
            f"{b}:{'/'.join(IMPL_NAMES[i] for i in conf.impls[b])}" for b in SMALL_FIELDS))
        conf.check_cantor()
        conf.small_decode(6 * args.iters)
        conf.small_roots(4 * args.iters)
        conf.large_16(args.iters)
        conf.large_other(args.iters)
        conf.large_32(args.iters)
    except Mismatch as e:
        print(f"MISMATCH (seed {seed:#x}): {e}")
        sys.exit(1)
    print(f"conformance OK in {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
