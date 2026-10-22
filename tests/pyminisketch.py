#!/usr/bin/env python3
# Copyright (c) 2020 Pieter Wuille
# Distributed under the MIT software license, see the accompanying
# file LICENSE or http://www.opensource.org/licenses/mit-license.php.

"""Native Python (slow) reimplementation of libminisketch' algorithms."""

import functools
import json
import math
import os
import random
import unittest


def setUpModule():
    """Seed the RNG reproducibly (unittest runs this before any test).

    The seed is printed and can be fixed with the MINISKETCH_TEST_SEED
    environment variable to replay a failing run.
    """
    seed = os.environ.get('MINISKETCH_TEST_SEED')
    seed = int(seed, 0) if seed else random.randrange(1 << 63)
    print(f"Test seed: {seed:#x} (reproduce with MINISKETCH_TEST_SEED={seed:#x})")
    random.seed(seed)

# Irreducible polynomials over GF(2) to use (represented as integers).
#
# Most fields can be defined by multiple such polynomials. Minisketch uses the one with the minimal
# number of nonzero coefficients, and tie-breaking by picking the lexicographically first among
# those.
#
# All polynomials for degrees 2 through 64 (inclusive) are given.
GF2_MODULI = [
    None, None,
    2**2 + 2**1 + 1,
    2**3 + 2**1 + 1,
    2**4 + 2**1 + 1,
    2**5 + 2**2 + 1,
    2**6 + 2**1 + 1,
    2**7 + 2**1 + 1,
    2**8 + 2**4 + 2**3 + 2**1 + 1,
    2**9 + 2**1 + 1,
    2**10 + 2**3 + 1,
    2**11 + 2**2 + 1,
    2**12 + 2**3 + 1,
    2**13 + 2**4 + 2**3 + 2**1 + 1,
    2**14 + 2**5 + 1,
    2**15 + 2**1 + 1,
    2**16 + 2**5 + 2**3 + 2**1 + 1,
    2**17 + 2**3 + 1,
    2**18 + 2**3 + 1,
    2**19 + 2**5 + 2**2 + 2**1 + 1,
    2**20 + 2**3 + 1,
    2**21 + 2**2 + 1,
    2**22 + 2**1 + 1,
    2**23 + 2**5 + 1,
    2**24 + 2**4 + 2**3 + 2**1 + 1,
    2**25 + 2**3 + 1,
    2**26 + 2**4 + 2**3 + 2**1 + 1,
    2**27 + 2**5 + 2**2 + 2**1 + 1,
    2**28 + 2**1 + 1,
    2**29 + 2**2 + 1,
    2**30 + 2**1 + 1,
    2**31 + 2**3 + 1,
    2**32 + 2**7 + 2**3 + 2**2 + 1,
    2**33 + 2**10 + 1,
    2**34 + 2**7 + 1,
    2**35 + 2**2 + 1,
    2**36 + 2**9 + 1,
    2**37 + 2**6 + 2**4 + 2**1 + 1,
    2**38 + 2**6 + 2**5 + 2**1 + 1,
    2**39 + 2**4 + 1,
    2**40 + 2**5 + 2**4 + 2**3 + 1,
    2**41 + 2**3 + 1,
    2**42 + 2**7 + 1,
    2**43 + 2**6 + 2**4 + 2**3 + 1,
    2**44 + 2**5 + 1,
    2**45 + 2**4 + 2**3 + 2**1 + 1,
    2**46 + 2**1 + 1,
    2**47 + 2**5 + 1,
    2**48 + 2**5 + 2**3 + 2**2 + 1,
    2**49 + 2**9 + 1,
    2**50 + 2**4 + 2**3 + 2**2 + 1,
    2**51 + 2**6 + 2**3 + 2**1 + 1,
    2**52 + 2**3 + 1,
    2**53 + 2**6 + 2**2 + 2**1 + 1,
    2**54 + 2**9 + 1,
    2**55 + 2**7 + 1,
    2**56 + 2**7 + 2**4 + 2**2 + 1,
    2**57 + 2**4 + 1,
    2**58 + 2**19 + 1,
    2**59 + 2**7 + 2**4 + 2**2 + 1,
    2**60 + 2**1 + 1,
    2**61 + 2**5 + 2**2 + 2**1 + 1,
    2**62 + 2**29 + 1,
    2**63 + 2**1 + 1,
    2**64 + 2**4 + 2**3 + 2**1 + 1
]

class GF2Ops:
    """Class to perform GF(2^field_size) operations on elements represented as integers.

    Given that elements are represented as integers, addition is simply xor, and not
    exposed here.
    """

    def __init__(self, field_size):
        """Construct a GF2Ops object for the specified field size."""
        self.field_size = field_size
        self._modulus = GF2_MODULI[field_size]
        assert self._modulus is not None

    def mul2(self, x):
        """Multiply x by 2 in GF(2^field_size)."""
        x <<= 1
        if x >> self.field_size:
            x ^= self._modulus
        return x

    def mul(self, x, y):
        """Multiply x by y in GF(2^field_size)."""
        ret = 0
        while y:
            if y & 1:
                ret ^= x
            y >>= 1
            x = self.mul2(x)
        return ret

    def sqr(self, x):
        """Square x in GF(2^field_size)."""
        return self.mul(x, x)

    def inv(self, x):
        """Compute the inverse of x in GF(2^field_size)."""
        assert x != 0
        # Use the extended polynomial Euclidean GCD algorithm on (modulus, x), over GF(2).
        # See https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor.
        t1, t2 = 0, 1
        r1, r2 = self._modulus, x
        r1l, r2l = self.field_size + 1, r2.bit_length()
        while r2:
            q = r1l - r2l
            r1 ^= r2 << q
            t1 ^= t2 << q
            r1l = r1.bit_length()
            if r1 < r2:
                t1, t2 = t2, t1
                r1, r2 = r2, r1
                r1l, r2l = r2l, r1l
        assert r1 == 1
        return t1

class TestGF2Ops(unittest.TestCase):
    """Test class for basic arithmetic properties of GF2Ops."""

    def field_size_test(self, field_size):
        """Test operations for given field_size."""

        gf = GF2Ops(field_size)
        for i in range(100):
            x = random.randrange(1 << field_size)
            y = random.randrange(1 << field_size)
            x2 = gf.mul2(x)
            xy = gf.mul(x, y)
            self.assertEqual(x2, gf.mul(x, 2)) # mul2(x) == x*2
            self.assertEqual(x2, gf.mul(2, x)) # mul2(x) == 2*x
            self.assertEqual(xy == 0, x == 0 or y == 0)
            self.assertEqual(xy == x, y == 1 or x == 0)
            self.assertEqual(xy == y, x == 1 or y == 0)
            self.assertEqual(xy, gf.mul(y, x)) # x*y == y*x
            if i < 10:
                xp = x
                for _ in range(field_size):
                    xp = gf.sqr(xp)
                self.assertEqual(xp, x) # x^(2^field_size) == x
            if y != 0:
                yi = gf.inv(y)
                self.assertEqual(y == yi, y == 1) # y==1/x iff y==1
                self.assertEqual(gf.mul(y, yi), 1) # y*(1/y) == 1
                yii = gf.inv(yi)
                self.assertEqual(y, yii) # 1/(1/y) == y
                if x != 0:
                    xi = gf.inv(x)
                    xyi = gf.inv(xy)
                    self.assertEqual(xyi, gf.mul(xi, yi)) # (1/x)*(1/y) == 1/(x*y)

    def test(self):
        """Run tests."""
        for field_size in range(2, 65):
            self.field_size_test(field_size)

# The operations below operate on polynomials over GF(2^field_size), represented as lists of
# integers:
#
#   [a, b, c, ...] = a + b*x + c*x^2 + ...
#
# As an invariant, there are never any trailing zeroes in the list representation.
#
# Examples:
# * [] = 0
# * [3] = 3
# * [0, 1] = x
# * [2, 0, 5] = 5*x^2 + 2

def poly_monic(poly, gf):
    """Return a monic version of the polynomial poly."""
    # Multiply every coefficient with the inverse of the top coefficient.
    inv = gf.inv(poly[-1])
    return [gf.mul(inv, v) for v in poly]

def poly_divmod(poly, mod, gf):
    """Return the polynomial (quotient, remainder) of poly divided by mod."""
    assert len(mod) > 0 and mod[-1] == 1 # Require monic mod.
    if len(poly) < len(mod):
        return ([], poly)
    val = list(poly)
    div = [0 for _ in range(len(val) - len(mod) + 1)]
    while len(val) >= len(mod):
        term = val[-1]
        div[len(val) - len(mod)] = term
        # If the highest coefficient in val is nonzero, subtract a multiple of mod from it.
        val.pop()
        if term != 0:
            for x in range(len(mod) - 1):
                val[1 + x - len(mod)] ^= gf.mul(term, mod[x])
    # Prune trailing zero coefficients.
    while len(val) > 0 and val[-1] == 0:
        val.pop()
    return div, val

def poly_gcd(a, b, gf):
    """Return the polynomial GCD of a and b."""
    if len(a) < len(b):
        a, b = b, a
    # Use Euclid's algorithm to find the GCD of a and b.
    # see https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Euclid's_algorithm.
    while len(b) > 0:
        b = poly_monic(b, gf)
        (_, b), a = poly_divmod(a, b, gf), b
    return a

def poly_sqr(poly, gf):
    """Return the square of polynomial poly."""
    if len(poly) == 0:
        return []
    # In characteristic-2 fields, thanks to Frobenius' endomorphism ((a + b)^2 = a^2 + b^2),
    # squaring a polynomial is easy: square all the coefficients and interleave with zeroes.
    # E.g., (3 + 5*x + 17*x^2)^2 = 3^2 + (5*x)^2 + (17*x^2)^2.
    # See https://en.wikipedia.org/wiki/Frobenius_endomorphism.
    return [0 if i & 1 else gf.sqr(poly[i // 2]) for i in range(2 * len(poly) - 1)]

def poly_tracemod(poly, param, gf):
    """Compute y + y^2 + y^4 + ... + y^(2^(field_size-1)) mod poly, where y = param*x."""
    out = [0, param]
    for _ in range(gf.field_size - 1):
        # In each loop iteration i, we start with out = y + y^2 + ... + y^(2^i). By squaring that we
        # transform it into out = y^2 + y^4 + ... + y^(2^(i+1)).
        out = poly_sqr(out, gf)
        # Thus, we just need to add y again to it to get out = y + ... + y^(2^(i+1)).
        while len(out) < 2:
            out.append(0)
        out[1] = param
        # Finally take a modulus to keep the intermediary polynomials small.
        _, out = poly_divmod(out, poly, gf)
    return out

def poly_frobeniusmod(poly, gf):
    """Compute x^(2^field_size) mod poly."""
    out = [0, 1]
    for _ in range(gf.field_size):
        _, out = poly_divmod(poly_sqr(out, gf), poly, gf)
    return out

def poly_find_roots(poly, gf):
    """Find the roots of poly if fully factorizable with unique roots, [] otherwise."""
    assert len(poly) > 0
    # If the polynomial is constant (and nonzero), it has no roots.
    if len(poly) == 1:
        return []
    # Make the polynomial monic (which doesn't change its roots).
    poly = poly_monic(poly, gf)
    # If the polynomial is of the form x+a, return a.
    if len(poly) == 2:
        return [poly[0]]
    # Otherwise, first test that poly can be completely factored into unique roots. The polynomial
    # x^(2^fieldsize)-x has every field element once as root. Thus we want to know that that is a
    # multiple of poly. Compute x^(field_size) mod poly, which needs to equal x if that is the case
    # (unless poly has degree <= 1, but that case is handled above).
    if poly_frobeniusmod(poly, gf) != [0, 1]:
        return []

    def rec_split(poly, randv):
        """Recursively split poly using the Berlekamp trace algorithm."""
        # See https://hal.archives-ouvertes.fr/hal-00626997/document.
        assert len(poly) > 1 and poly[-1] == 1 # Require a monic poly.
        # If poly is of the form x+a, its root is a.
        if len(poly) == 2:
            return [poly[0]]
        # Try consecutive randomization factors randv, until one is found that factors poly.
        while True:
            # Compute the trace of (randv*x) mod poly. This is a polynomial that maps half of the
            # domain to 0, and the other half to 1. Which half that is is controlled by randv.
            # By taking it modulo poly, we only add a multiple of poly. Thus the result has at least
            # the shared roots of the trace polynomial and poly still, but may have others.
            trace = poly_tracemod(poly, randv, gf)
            # Using the set {2^i*a for i=0..fieldsize-1} gives optimally independent randv values
            # (no more than fieldsize are ever needed).
            randv = gf.mul2(randv)
            # Now take the GCD of this trace polynomial with poly. The result is a polynomial
            # that only has the shared roots of the trace polynomial and poly as roots.
            gcd = poly_gcd(trace, poly, gf)
            # If the result has a degree higher than 1, and lower than that of poly, we found a
            # useful factorization.
            if len(gcd) != len(poly) and len(gcd) > 1:
                break
            # Otherwise, continue with another randv.
        # Find the actual factors: the monic version of the GCD above, and poly divided by it.
        factor1 = poly_monic(gcd, gf)
        factor2, _ = poly_divmod(poly, gcd, gf)
        # Recurse.
        return rec_split(factor1, randv) + rec_split(factor2, randv)

    # Invoke the recursive splitting with a random initial factor, and sort the results.
    return sorted(rec_split(poly, random.randrange(1, 1 << gf.field_size)))

class TestPolyFindRoots(unittest.TestCase):
    """Test class for poly_find_roots."""

    def field_size_test(self, field_size):
        """Run tests for given field_size."""
        gf = GF2Ops(field_size)
        for test_size in [0, 1, 2, 3, 10]:
            roots = [random.randrange(1 << field_size) for _ in range(test_size)]
            roots_set = set(roots)
            # Construct a polynomial with all elements of roots as roots (with multiplicity).
            poly = [1]
            for root in roots:
                new_poly = [0] + poly
                for n, c in enumerate(poly):
                    new_poly[n] ^= gf.mul(c, root)
                poly = new_poly
            # Invoke the root finding algorithm.
            found_roots = poly_find_roots(poly, gf)
            # The result must match the input, unless any roots were repeated.
            if len(roots) == len(roots_set):
                self.assertEqual(found_roots, sorted(roots))
            else:
                self.assertEqual(found_roots, [])

    def test(self):
        """Run tests."""
        for field_size in range(2, 65):
            self.field_size_test(field_size)

def berlekamp_massey(syndromes, gf):
    """Implement the Berlekamp-Massey algorithm.

    Takes as input a sequence of GF(2^field_size) elements, and returns the shortest LSFR
    that generates it, represented as a polynomial.
    """
    # See https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm.
    current = [1]
    prev = [1]
    b_inv = 1
    for n, discrepancy in enumerate(syndromes):
        # Compute discrepancy
        for i in range(1, len(current)):
            discrepancy ^= gf.mul(syndromes[n - i], current[i])

        # Correct if discrepancy is nonzero.
        if discrepancy:
            x = n + 1 - (len(current) - 1) - (len(prev) - 1)
            if 2 * (len(current) - 1) <= n:
                tmp = list(current)
                current.extend(0 for _ in range(len(prev) + x - len(current)))
                mul = gf.mul(discrepancy, b_inv)
                for i, v in enumerate(prev):
                    current[i + x] ^= gf.mul(mul, v)
                prev = tmp
                b_inv = gf.inv(discrepancy)
            else:
                mul = gf.mul(discrepancy, b_inv)
                for i, v in enumerate(prev):
                    current[i + x] ^= gf.mul(mul, v)
    return current

def poly_mul(a, b, gf):
    """Multiply two polynomials (schoolbook)."""
    if len(a) == 0 or len(b) == 0:
        return []
    ret = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            ret[i + j] ^= gf.mul(av, bv)
    return ret

class TestBerlekampMassey(unittest.TestCase):
    """Test class for berlekamp_massey."""

    def test_zero_syndromes(self):
        """All-zero syndromes are generated by the constant polynomial [1]."""
        gf = GF2Ops(11)
        self.assertEqual(berlekamp_massey([0] * 8, gf), [1])

    def test_known_roots(self):
        """Syndromes from a known root set yield the polynomial prod (1 + m*x)."""
        for field_size in [2, 3, 8, 11, 16, 24, 32, 48, 64]:
            gf = GF2Ops(field_size)
            for _ in range(4):
                count = random.randrange(min(8, (1 << field_size) - 1) + 1)
                roots = []
                while len(roots) < count:
                    r = random.randrange(1, 1 << field_size)
                    if r not in roots:
                        roots.append(r)
                cap = count + random.randrange(3)
                if cap == 0:
                    continue
                # syndromes[k-1] = sum of m^k over the roots, k = 1..2*cap.
                powers = {m: 1 for m in roots}
                syndromes = []
                for _ in range(2 * cap):
                    s = 0
                    for m in roots:
                        powers[m] = gf.mul(powers[m], m)
                        s ^= powers[m]
                    syndromes.append(s)
                expect = [1]
                for m in roots:
                    expect = poly_mul(expect, [1, m], gf)
                self.assertEqual(berlekamp_massey(syndromes, gf), expect)

    def test_recurrence(self):
        """On any syndrome sequence with the Frobenius structure the output
        satisfies the defining LFSR recurrence."""
        for field_size in [2, 5, 11, 32]:
            gf = GF2Ops(field_size)
            for _ in range(4):
                odd = [random.randrange(1 << field_size) for _ in range(random.randrange(1, 9))]
                syndromes = [0] * (2 * len(odd))
                for i, v in enumerate(odd):
                    syndromes[2 * i] = v
                    syndromes[2 * i + 1] = gf.sqr(syndromes[i])
                conn = berlekamp_massey(syndromes, gf)
                self.assertGreater(len(conn), 0)
                self.assertEqual(conn[0], 1)
                for n in range(len(conn) - 1, len(syndromes)):
                    acc = 0
                    for i, c in enumerate(conn):
                        acc ^= gf.mul(c, syndromes[n - i])
                    self.assertEqual(acc, 0)

class Minisketch:
    """A Minisketch sketch.

    This represents a sketch of a certain capacity, with elements of a certain bit size.
    """

    def __init__(self, field_size, capacity):
        """Initialize an empty sketch with the specified field_size size and capacity."""
        self.field_size = field_size
        self.capacity = capacity
        self.odd_syndromes = [0] * capacity
        self.gf = GF2Ops(field_size)

    def add(self, element):
        """Add an element to this sketch. 1 <= element < 2**field_size."""
        sqr = self.gf.sqr(element)
        for pos in range(self.capacity):
            self.odd_syndromes[pos] ^= element
            element = self.gf.mul(sqr, element)

    def serialized_size(self):
        """Compute how many bytes a serialization of this sketch will be in size."""
        return (self.capacity * self.field_size + 7) // 8

    def serialize(self):
        """Serialize this sketch to bytes."""
        val = 0
        for i in range(self.capacity):
            val |= self.odd_syndromes[i] << (self.field_size * i)
        return val.to_bytes(self.serialized_size(), 'little')

    def deserialize(self, byte_data):
        """Deserialize a byte array into this sketch, overwriting its contents."""
        assert len(byte_data) == self.serialized_size()
        val = int.from_bytes(byte_data, 'little')
        for i in range(self.capacity):
            self.odd_syndromes[i] = (val >> (self.field_size * i)) & ((1 << self.field_size) - 1)

    def clone(self):
        """Return a clone of this sketch."""
        ret = Minisketch(self.field_size, self.capacity)
        ret.odd_syndromes = list(self.odd_syndromes)
        ret.gf = self.gf
        return ret

    def merge(self, other):
        """Merge a sketch with another sketch. Corresponds to XOR'ing their serializations."""
        assert self.capacity == other.capacity
        assert self.field_size == other.field_size
        for i in range(self.capacity):
            self.odd_syndromes[i] ^= other.odd_syndromes[i]

    def decode(self, max_count=None):
        """Decode the contents of this sketch.

        Returns either a list of elements or None if undecodable.
        """
        # We know the odd syndromes s1=x+y+..., s3=x^3+y^3+..., s5=..., and reconstruct the even
        # syndromes from this:
        #  * s2 = x^2+y^2+.... = (x+y+...)^2 = s1^2
        #  * s4 = x^4+y^4+.... = (x^2+y^2+...)^2 = s2^2
        #  * s6 = x^6+y^6+.... = (x^3+y^3+...)^2 = s3^2
        all_syndromes = [0 for _ in range(2 * len(self.odd_syndromes))]
        for i in range(len(self.odd_syndromes)):
            all_syndromes[i * 2] = self.odd_syndromes[i]
            all_syndromes[i * 2 + 1] = self.gf.sqr(all_syndromes[i])
        # Given the syndromes, find the polynomial that generates them.
        poly = berlekamp_massey(all_syndromes, self.gf)
        # Deal with failure and trivial cases.
        if len(poly) == 0:
            return None
        if len(poly) == 1:
            return []
        if max_count is not None and len(poly) > 1 + max_count:
            return None
        # If the polynomial can be factored into (1-m1*x)*(1-m2*x)*...*(1-mn*x), then {m1,m2,...,mn}
        # is our set. As each factor (1-m*x) has 1/m as root, we're really just looking for the
        # inverses of the roots. We find these by reversing the order of the coefficients, and
        # finding the roots.
        roots = poly_find_roots(list(reversed(poly)), self.gf)
        if len(roots) == 0:
            return None
        return roots

class TestMinisketch(unittest.TestCase):
    """Test class for Minisketch."""

    @classmethod
    def construct_data(cls, field_size, num_a_only, num_b_only, num_both):
        """Construct two random lists of elements in [1..2**field_size-1].

        Each list will have unique elements that don't appear in the other (num_a_only in the first
        and num_b_only in the second), and num_both elements will appear in both."""
        sample = []
        # Simulate random.sample here (which doesn't work with ranges over 2**63).
        for _ in range(num_a_only + num_b_only + num_both):
            while True:
                r = random.randrange(1, 1 << field_size)
                if r not in sample:
                    sample.append(r)
                    break
        full_a = sample[:num_a_only + num_both]
        full_b = sample[num_a_only:]
        random.shuffle(full_a)
        random.shuffle(full_b)
        return full_a, full_b

    def field_size_capacity_test(self, field_size, capacity):
        """Test Minisketch methods for a specific field and capacity."""
        used_capacity = random.randrange(capacity + 1)
        num_a = random.randrange(used_capacity + 1)
        num_both = random.randrange(min(2 * capacity, (1 << field_size) - 1 - used_capacity) + 1)
        full_a, full_b = self.construct_data(field_size, num_a, used_capacity - num_a, num_both)
        sketch_a = Minisketch(field_size, capacity)
        sketch_b = Minisketch(field_size, capacity)
        for v in full_a:
            sketch_a.add(v)
        for v in full_b:
            sketch_b.add(v)
        sketch_combined = sketch_a.clone()
        sketch_b_ser = sketch_b.serialize()
        sketch_b_received = Minisketch(field_size, capacity)
        sketch_b_received.deserialize(sketch_b_ser)
        sketch_combined.merge(sketch_b_received)
        decode = sketch_combined.decode()
        self.assertEqual(decode, sorted(set(full_a) ^ set(full_b)))

    def test(self):
        """Run tests."""
        for field_size in range(2, 65):
            for capacity in [0, 1, 2, 5, 10, field_size]:
                self.field_size_capacity_test(field_size, min(capacity, (1 << field_size) - 1))

class TestMinisketchSerialization(unittest.TestCase):
    """Serialization, merge, and decode-limit behavior of Minisketch."""

    def test_roundtrip_and_merge(self):
        """Serialization roundtrips; merging any split equals the whole (= XOR of parts)."""
        for field_size in range(2, 65):
            capacity = random.randrange(0, 9)
            elements = [random.randrange(1, 1 << field_size)
                        for _ in range(random.randrange(0, 2 * capacity + 1))]
            sketch = Minisketch(field_size, capacity)
            for v in elements:
                sketch.add(v)
            ser = sketch.serialize()
            self.assertEqual(len(ser), sketch.serialized_size())

            sketch2 = Minisketch(field_size, capacity)
            sketch2.deserialize(ser)
            self.assertEqual(sketch2.serialize(), ser)
            self.assertEqual(sketch2.odd_syndromes, sketch.odd_syndromes)

            split = random.randrange(0, len(elements) + 1)
            part_a = Minisketch(field_size, capacity)
            part_b = Minisketch(field_size, capacity)
            for v in elements[:split]:
                part_a.add(v)
            for v in elements[split:]:
                part_b.add(v)
            ser_a, ser_b = part_a.serialize(), part_b.serialize()
            part_a.merge(part_b)
            self.assertEqual(part_a.serialize(), ser)
            self.assertEqual(bytes(x ^ y for x, y in zip(ser_a, ser_b)), ser)

    def test_decode_max_count(self):
        """Decoding with max_count below the element count fails; at the count it succeeds."""
        for field_size in [4, 8, 11, 16]:
            capacity = random.randrange(1, 9)
            count = random.randrange(1, min(capacity, (1 << field_size) - 1) + 1)
            elements = []
            while len(elements) < count:
                r = random.randrange(1, 1 << field_size)
                if r not in elements:
                    elements.append(r)
            sketch = Minisketch(field_size, capacity)
            for v in elements:
                sketch.add(v)
            self.assertEqual(sketch.decode(count), sorted(elements))
            self.assertIsNone(sketch.decode(count - 1))

################################################################################
# Optimized decode (oracle)
#
# One short function per optimization in libminisketch's decoder, each next to
# the textbook function it replaces and tested against it. The claim IDs (C13
# etc.) refer to doc/spec/claims.md, where the Lean proofs are listed. This is
# meant to be read: it keeps the structure of each technique and none of the
# C++ engineering (lazy reduction, tiers, cutoffs, scratch buffers).
################################################################################

class TableGF(GF2Ops):
    """GF2Ops with log/exp tables, for field sizes up to 16 (fast enough to
    exercise the large-degree root-finding paths in Python)."""

    def __init__(self, field_size):
        super().__init__(field_size)
        assert field_size <= 16
        order = (1 << field_size) - 1
        for gen in range(2, 1 << field_size):
            exp = [1] * (2 * order)
            seen = {1}
            for i in range(1, order):
                exp[i] = GF2Ops.mul(self, exp[i - 1], gen)
                seen.add(exp[i])
            if len(seen) == order:
                break
        for i in range(order, 2 * order):
            exp[i] = exp[i - order]
        self._exp = exp
        self._log = [0] * (1 << field_size)
        for i in range(order):
            self._log[exp[i]] = i
        self._order = order

    def mul(self, x, y):
        if x == 0 or y == 0:
            return 0
        return self._exp[self._log[x] + self._log[y]]

    def inv(self, x):
        assert x != 0
        return self._exp[(self._order - self._log[x]) % self._order]

def make_gf(field_size):
    """The fastest exact GF(2^field_size) arithmetic available here."""
    return TableGF(field_size) if field_size <= 16 else GF2Ops(field_size)

# Small polynomial helpers (same representation as above: low to high, trimmed).

def poly_trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p

def poly_add(a, b):
    if len(a) < len(b):
        a, b = b, a
    out = list(a)
    for i, v in enumerate(b):
        out[i] ^= v
    return poly_trim(out)

def poly_scale(p, c, gf):
    return poly_trim([gf.mul(c, v) for v in p])

def poly_eval(p, x, gf):
    acc = 0
    for v in reversed(p):
        acc = gf.mul(acc, x) ^ v
    return acc

def poly_deriv(p):
    """Formal derivative in characteristic 2: only odd exponents survive."""
    return poly_trim([p[i] if i & 1 else 0 for i in range(1, len(p))])

def poly_from_roots(roots, gf):
    """prod (x - r)."""
    p = [1]
    for r in roots:
        p = poly_mul(p, [r, 1], gf)
    return p

def elem_trace(a, gf):
    """Absolute trace a + a^2 + ... + a^(2^(b-1)), which is 0 or 1."""
    t, y = 0, a
    for _ in range(gf.field_size):
        t ^= y
        y = gf.sqr(y)
    return t

def gf2_rank(vectors):
    """Rank over GF(2) of a list of integers seen as bit vectors."""
    basis = []
    for v in vectors:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return len(basis)

def solve_artin_schreier(c, gf):
    """A root y of y^2 + y = c, or None (y -> y^2 + y is GF(2)-linear)."""
    b = gf.field_size
    cols = [gf.sqr(1 << i) ^ (1 << i) for i in range(b)]
    # Gaussian elimination on the b x b system sum_i y_i * cols[i] = c.
    rows = []  # (pivot value, combination of columns as a bitmask)
    for i, v in enumerate(cols):
        comb = 1 << i
        for pv, pc in rows:
            if v ^ pv < v:
                v, comb = v ^ pv, comb ^ pc
        if v:
            rows.append((v, comb))
            rows.sort(reverse=True)
    y = 0
    for pv, pc in rows:
        if c ^ pv < c:
            c, y = c ^ pv, y ^ pc
    if c:
        return None
    return y

# --- Berlekamp-Massey on the odd syndromes (C12, C13-C15, C48) ---

def verify_tail_odd(conn, odd, processed, gf):
    """C15: given conn (Lambda, conn[0] = 1, length L+1) generating
    s_1..s_processed with processed >= 2L, check that it generates all of
    s_1..s_{2c} using only the stored odd syndromes u_j = s_{2j+1} = odd[j]:
    u_j + sum_{a>=1} Lambda_a^2 u_{j-a} = 0 for every odd 2j+1 > processed."""
    sq = [gf.sqr(c) for c in conn]
    for j in range(len(odd)):
        if 2 * j + 1 <= processed:
            continue
        acc = odd[j]
        for a in range(1, len(conn)):
            acc ^= gf.mul(sq[a], odd[j - a])
        if acc:
            return False
    return True

def bm_odd(odd, max_degree, gf):
    """BerlekampMasseyOdd: textbook BM on s_1..s_{2c} that (C12) skips the even
    steps, whose discrepancy is always 0 on Frobenius-structured input, and
    (C48) stops early once a bulk check of the remaining odd syndromes passes.
    Returns the connection polynomial, or None if its degree would exceed
    max_degree (the same result as textbook BM on all 2c syndromes)."""
    n_syn = 2 * len(odd)
    syn = [0] * n_syn
    for i, v in enumerate(odd):
        syn[2 * i] = v
        syn[2 * i + 1] = gf.sqr(syn[i])
    current, prev, b_inv = [1], [1], 1
    next_check = 0
    for n in range(0, n_syn, 2):  # odd syndromes s_{n+1} only
        discrepancy = syn[n]
        for i in range(1, len(current)):
            discrepancy ^= gf.mul(syn[n - i], current[i])
        if discrepancy:
            x = n + 1 - (len(current) - 1) - (len(prev) - 1)
            mul = gf.mul(discrepancy, b_inv)
            if 2 * (len(current) - 1) <= n:
                tmp = list(current)
                current.extend(0 for _ in range(len(prev) + x - len(current)))
                for i, v in enumerate(prev):
                    current[i + x] ^= gf.mul(mul, v)
                prev, b_inv = tmp, gf.inv(discrepancy)
            else:
                for i, v in enumerate(prev):
                    current[i + x] ^= gf.mul(mul, v)
            if len(current) - 1 > max_degree:
                return None
        else:
            # Steps 1..n+1 are generated by current, and step n+2 (even) is too.
            processed = n + 2
            if processed >= 2 * (len(current) - 1) and processed >= next_check:
                if verify_tail_odd(current, odd, processed, gf):
                    return current
                next_check = processed + processed // 2 + 2
    return current

def bm_reference(odd, max_degree, gf):
    """Textbook BM on all 2c syndromes, with the max_degree failure."""
    syn = [0] * (2 * len(odd))
    for i, v in enumerate(odd):
        syn[2 * i] = v
        syn[2 * i + 1] = gf.sqr(syn[i])
    conn = berlekamp_massey(syn, gf)
    return None if len(conn) - 1 > max_degree else conn

# --- Transposed Karatsuba middle product (C17) ---

def middle_product_naive(c, w, gf):
    """MP(c, w)_i = sum_j c_j w_{i+j} for i < n = len(c), len(w) = 2n - 1."""
    n = len(c)
    return [functools.reduce(lambda acc, j: acc ^ gf.mul(c[j], w[i + j]), range(n), 0)
            for i in range(n)]

def middle_product_kara(c, w, gf):
    """Three half-size middle products instead of four (C17):
    low = alpha + MP(c0, A + B), high = alpha + MP(c1, B + D),
    alpha = MP(c0 + c1, B), with w split into overlapping windows A, B, D."""
    n = len(c)
    if n <= 2 or n & 1:
        return middle_product_naive(c, w, gf)
    h = n // 2
    c0, c1 = c[:h], c[h:]
    A, B, D = w[0:2 * h - 1], w[h:3 * h - 1], w[2 * h:4 * h - 1]
    alpha = middle_product_kara([x ^ y for x, y in zip(c0, c1)], B, gf)
    low = middle_product_kara(c0, [x ^ y for x, y in zip(A, B)], gf)
    high = middle_product_kara(c1, [x ^ y for x, y in zip(B, D)], gf)
    return [x ^ y for x, y in zip(alpha, low)] + [x ^ y for x, y in zip(alpha, high)]

# --- Shared Frobenius orbit and split parameters (C22, C26, C27) ---

def frobenius_orbit(f, gf):
    """X_i = x^(2^i) mod f for i = 0..b. One orbit at the root serves every node
    of the split tree, since (h mod f) mod g = h mod g for g | f (C26), and
    X_b = x iff f is square-free and splits (C22)."""
    orbit = [poly_divmod([0, 1], f, gf)[1]]
    for _ in range(gf.field_size):
        orbit.append(poly_divmod(poly_sqr(orbit[-1], gf), f, gf)[1])
    return orbit

def trace_from_orbit(orbit, beta, g, gf):
    """Tr(beta x) mod g = sum_i beta^(2^i) X_i, reduced mod g."""
    acc, bp = [], beta
    for i in range(gf.field_size):
        acc = poly_add(acc, poly_scale(orbit[i], bp, gf))
        bp = gf.sqr(bp)
    return poly_divmod(acc, g, gf)[1]

def subfield_generator(sub, gf):
    """An element generating GF(2^sub) inside GF(2^b) (sub a power of two
    dividing b): the relative trace of some x^k into GF(2^sub) that is not in
    the next subfield down."""
    b = gf.field_size
    xk = 2
    while True:
        w, y = 0, xk
        for _ in range(b // sub):
            w ^= y
            for _ in range(sub):
                y = gf.sqr(y)
        z = w
        for _ in range(sub // 2):
            z = gf.sqr(z)
        if z != w:
            return w
        xk = gf.mul(xk, 2)

def tower_betas(basis, gf):
    """Split parameters (C27). With 8 | b: beta_{8l+j} = basis * x^l * z^j for
    z generating GF(2^8), l < b/8, j < 8; otherwise beta_k = basis * x^k. Either
    way they are a GF(2)-basis, which is all the trace algorithm needs (C25)."""
    b = gf.field_size
    if b % 8:
        betas, v = [], basis
        for _ in range(b):
            betas.append(v)
            v = gf.mul2(v)
        return betas
    z = subfield_generator(8, gf)
    betas = []
    xl = basis
    for _ in range(b // 8):
        zj = xl
        for _ in range(8):
            betas.append(zj)
            zj = gf.mul(zj, z)
        xl = gf.mul2(xl)
    return betas

# --- Additive FFT on Cantor-basis cosets (C29, C34, C38) ---

def cantor_basis(gf):
    """beta_0 = 1, beta_{j+1}^2 + beta_{j+1} = beta_j, for b a power of two (C38:
    the equation stays solvable because Tr(beta_j) = 0 for j < b - 1)."""
    b = gf.field_size
    assert b & (b - 1) == 0
    beta = [1]
    for _ in range(b - 1):
        beta.append(solve_artin_schreier(beta[-1], gf))
        assert beta[-1] is not None
    return beta

def subspace_poly(m):
    """S^m as a polynomial with 0/1 coefficients: S^0 = x, S^{m+1} = S^m^2 + S^m,
    so S^m = sum_{j subset of m} x^(2^j) (Lucas)."""
    p = [0] * ((1 << m) + 1)
    for j in range(m + 1):
        if j & m == j:
            p[1 << j] = 1
    return p

def subspace_eval(m, x, gf):
    for _ in range(m):
        x = gf.sqr(x) ^ x
    return x

def coset_points(beta, eta, m):
    """eta + span(beta_0..beta_{m-1}), point i = eta + sum of beta_j over bits j of i."""
    pts = [eta]
    for j in range(m):
        pts += [p ^ beta[j] for p in pts]
    return pts

def to_novel(f, m):
    """Coefficients in the novel basis X_k = prod_{j in bits(k)} S^j, k < 2^m:
    f = g0 + S^{m-1} g1, recursively. S^{m-1} has 0/1 coefficients."""
    f = list(f) + [0] * ((1 << m) - len(f))
    assert len(f) == 1 << m
    if m == 0:
        return f
    s = subspace_poly(m - 1)
    exps = [e for e, v in enumerate(s) if v and e < len(s) - 1]
    top = len(s) - 1
    val, q = list(f), [0] * (len(f) - top)
    for pos in range(len(f) - 1, top - 1, -1):
        term = val[pos]
        q[pos - top] = term
        val[pos] = 0
        if term:
            for e in exps:
                val[pos - top + e] ^= term
    return to_novel(val[:top], m - 1) + to_novel(q + [0] * (top - len(q)), m - 1)

def from_novel(d, gf):
    """Inverse of to_novel."""
    m = len(d).bit_length() - 1
    if m == 0:
        return poly_trim(d)
    h = len(d) // 2
    g0, g1 = from_novel(d[:h], gf), from_novel(d[h:], gf)
    return poly_add(g0, poly_mul(subspace_poly(m - 1), g1, gf))

def novel_fft(d, beta, eta, gf):
    """Values of sum_k d_k X_k at coset_points(beta, eta, m). On the two halves
    of the coset S^{m-1} is t = S^{m-1}(eta) and t + 1, so f = g0 + S^{m-1} g1 gives
    h0 = g0 + t g1 on the first half and h0 + g1 on the second."""
    m = len(d).bit_length() - 1
    if m == 0:
        return list(d)
    h = len(d) // 2
    t = subspace_eval(m - 1, eta, gf)
    h0 = [x ^ gf.mul(t, y) for x, y in zip(d[:h], d[h:])]
    h1 = [x ^ y for x, y in zip(h0, d[h:])]
    return novel_fft(h0, beta, eta, gf) + novel_fft(h1, beta, eta ^ beta[m - 1], gf)

def novel_ifft(v, beta, eta, gf):
    """Inverse of novel_fft."""
    m = len(v).bit_length() - 1
    if m == 0:
        return list(v)
    h = len(v) // 2
    h0 = novel_ifft(v[:h], beta, eta, gf)
    h1 = novel_ifft(v[h:], beta, eta ^ beta[m - 1], gf)
    t = subspace_eval(m - 1, eta, gf)
    g1 = [x ^ y for x, y in zip(h0, h1)]
    g0 = [x ^ gf.mul(t, y) for x, y in zip(h0, g1)]
    return g0 + g1

def additive_fft(f, beta, eta, m, gf):
    """Evaluate the polynomial f on coset_points(beta, eta, m). A polynomial of
    degree >= 2^m is first reduced mod the coset's vanishing polynomial
    Z_C = S^m + S^m(eta), which is 0 on C."""
    if len(f) > 1 << m:
        f = poly_divmod(f, poly_add(subspace_poly(m), [subspace_eval(m, eta, gf)]), gf)[1]
    return novel_fft(to_novel(f, m), beta, eta, gf)

def additive_ifft(vals, beta, eta, gf):
    """The polynomial of degree < len(vals) with these values on the coset."""
    return from_novel(novel_ifft(vals, beta, eta, gf), gf)

def novel_derivative(d):
    """C34: (S^j)' = 1, so X_k' = sum_{j in bits(k)} X_{k - 2^j}."""
    out = [0] * len(d)
    for k, v in enumerate(d):
        j = 0
        while (1 << j) <= k:
            if k >> j & 1:
                out[k - (1 << j)] ^= v
            j += 1
    return out

# --- Coset Montgomery reduction (C30, C31, C33) ---

def mont_square(rv, f_vals, fd_vals, beta, eta, gf):
    """One Montgomery squaring mod f on a coset C (C30). Input: values on C of a
    polynomial r (degree < N = |C|), and of f, f' (f nonzero on C). Output: the
    values on C of r_new with r_new * Z_C = r^2 + q f, q = r^2 / f on C, i.e.
    r_new = r^2 Z_C^-1 mod f, computed as q'(w) f(w) + q(w) f'(w)."""
    q_vals = [gf.mul(gf.sqr(r), gf.inv(fv)) for r, fv in zip(rv, f_vals)]
    qd = novel_fft(novel_derivative(novel_ifft(q_vals, beta, eta, gf)), beta, eta, gf)
    return [gf.mul(a, fv) ^ gf.mul(q, fdv) for a, fv, q, fdv in zip(qd, f_vals, q_vals, fd_vals)]

def mont_carry_to_factor(v_vals, g, beta, eta, gf):
    """C31: from the values of v on C = C0 + C1 (two halves), the values on C1
    of r = (v + q g) / Z_C0 with q = v / g on C0, using Z_C0 = 1 on C1. Then
    r * Z_C0 = v (mod g): the residue moves to the child coset C1."""
    h = len(v_vals) // 2
    m = h.bit_length() - 1
    g0 = additive_fft(g, beta, eta, m, gf)
    q = additive_ifft([gf.mul(a, gf.inv(b)) for a, b in zip(v_vals[:h], g0)], beta, eta, gf)
    eta1 = eta ^ beta[m]
    g1 = additive_fft(g, beta, eta1, m, gf)
    q1 = additive_fft(q, beta, eta1, m, gf)
    return [a ^ gf.mul(b, c) for a, b, c in zip(v_vals[h:], q1, g1)]

def coset_divide(g, h, beta, eta, m, gf):
    """C33: g / h by pointwise division on a coset of 2^m >= deg g points; the
    result is accepted only if it is monic of degree deg g - deg h, which
    certifies q h = g. Returns None otherwise (or if h has a root on C)."""
    gv, hv = additive_fft(g, beta, eta, m, gf), additive_fft(h, beta, eta, m, gf)
    if 0 in hv:
        return None
    q = additive_ifft([gf.mul(a, gf.inv(b)) for a, b in zip(gv, hv)], beta, eta, gf)
    if len(q) != len(g) - len(h) + 1 or q[-1] != 1:
        return None
    return q

def montgomery_orbit(f, beta, gf):
    """The Frobenius orbit X_0..X_b of f computed with Montgomery squarings on a
    coset C of 2^m >= deg f points on which f has no root. Returns the orbit in
    normal form and the coset (eta, m) used, or None if every candidate coset
    contains a root."""
    d = len(f) - 1
    m = max(1, (d - 1).bit_length())
    if m + 2 > gf.field_size:
        return None
    fd = poly_deriv(f)
    for blk in (1, 2, 3, 0):
        eta = (beta[m] if blk & 1 else 0) ^ (beta[m + 1] if blk & 2 else 0)
        pts = coset_points(beta, eta, m)
        f_vals = [poly_eval(f, w, gf) for w in pts]
        if 0 in f_vals:
            continue
        fd_vals = [poly_eval(fd, w, gf) for w in pts]
        z = poly_add(subspace_poly(m), [subspace_eval(m, eta, gf)])
        # Montgomery form of x: x * Z_C mod f, kept as values on C.
        cur = poly_divmod(poly_mul([0, 1], z, gf), f, gf)[1]
        vals = additive_fft(cur, beta, eta, m, gf)
        z_inv = None
        orbit = []
        for i in range(gf.field_size + 1):
            # Back to normal form: multiply by Z_C^-1 mod f (textbook, for clarity).
            rep = additive_ifft(vals, beta, eta, gf)
            if z_inv is None:
                z_inv = poly_inverse_mod(z, f, gf)
            orbit.append(poly_divmod(poly_mul(rep, z_inv, gf), f, gf)[1])
            vals = mont_square(vals, f_vals, fd_vals, beta, eta, gf)
        return orbit, eta, m
    return None

def poly_inverse_mod(a, f, gf):
    """a^-1 mod f by the extended Euclidean algorithm (a coprime to f)."""
    r0, r1, t0, t1 = f, poly_divmod(a, f, gf)[1], [], [1]
    while len(r1) > 1:
        lead_inv = gf.inv(r1[-1])
        q, r = poly_divmod(r0, poly_scale(r1, lead_inv, gf), gf)
        q = poly_scale(q, lead_inv, gf)
        r0, r1 = r1, r
        t0, t1 = t1, poly_add(t0, poly_mul(q, t1, gf))
    assert len(r1) == 1
    return poly_scale(t1, gf.inv(r1[0]), gf)

# --- Power-sum splits (C35) ---

def power_sum_split(q_sums, t, n, gf):
    """C35: with t in {0,1} on the n roots of the node and Q_k their power sums
    (q_sums[k] for k = 1..), the side A = {t = 1} has power sums
    P_k = sum_j t_j Q_{j+k}. BM on P_1..P_{2M} (2M >= n) then returns A's
    locator, provided the result divides the node polynomial (checked by the
    caller). Returns (P, locator or None)."""
    two_m = n + (n & 1)
    P = [0] * (two_m + 1)
    for k in range(1, two_m + 1):
        acc = 0
        for j, tj in enumerate(t):
            acc ^= gf.mul(tj, q_sums[j + k])
        P[k] = acc
    conn = bm_reference([P[k] for k in range(1, two_m + 1, 2)], n, gf)
    return P, conn

def power_sums(roots, count, gf):
    """[0, Q_1, ..., Q_count] for a root set."""
    out = [0] * (count + 1)
    for r in roots:
        p = 1
        for k in range(1, count + 1):
            p = gf.mul(p, r)
            out[k] ^= p
    return out

# --- The whole decoder (C46) ---

def find_roots_fast(f, gf, basis=1, q_sums=None, mont=True):
    """Roots of the monic f if it splits into distinct linear factors, else None.
    The textbook Berlekamp trace algorithm (poly_find_roots) with: the full-split
    test and one Frobenius orbit at the root (C22, C26), optionally by Montgomery
    squaring on a coset (C30); tower split parameters (C27); and power-sum splits
    (C35) instead of gcds when the node's power sums are known, checked by
    CosetDivide (C33)."""
    d = len(f) - 1
    if d == 0:
        return []
    if d == 1:
        return [f[0]]
    b = gf.field_size
    orbit = None
    beta_c = None
    root_coset = None  # a coset (eta, m) with 2^m >= deg f and no root of f
    if mont and b & (b - 1) == 0 and b >= 4:
        beta_c = cantor_basis(gf)
        res = montgomery_orbit(f, beta_c, gf)
        if res is not None:
            orbit, eta, m = res
            root_coset = (eta, m)
    if orbit is None:
        orbit = frobenius_orbit(f, gf)
    if orbit[b] != poly_trim([0, 1]):
        return None
    betas = tower_betas(basis, gf)
    assert gf2_rank(betas) == b

    def split_with_sums(g, t, sums):
        """Power-sum split of node g: BM on the side t = 1, then on the other
        side (power sums Q - P), accepting a result only if CosetDivide shows
        it divides g. Returns 'trivial', None (use the gcd), or
        (side, its power sums, g / side, their power sums)."""
        n = len(g) - 1
        if sums is None or root_coset is None or len(sums) < len(t) + n + 2:
            return None
        # The node's roots are roots of f, so the first 2^m points of the root's
        # coset (a sub-coset) hold none of them; 2^m >= n points suffice (C33).
        eta, m = root_coset[0], max(1, (n - 1).bit_length())
        P, _ = power_sum_split(sums, t, n, gf)
        Pc = [x ^ y for x, y in zip(sums, P)]
        for side, other in ((P, Pc), (Pc, P)):
            conn = bm_reference([side[k] for k in range(1, len(side), 2)], n, gf)
            if conn is None:
                continue
            if len(conn) == 1:
                return 'trivial'
            loc = list(reversed(conn))
            q = coset_divide(g, loc, beta_c, eta, m, gf)
            if q is not None:
                return (loc, side, q, other)
        return None

    def rec(g, depth, sums):
        n = len(g) - 1
        if n == 1:
            return [g[0]]
        while True:
            assert depth < b  # C24/C25: independent betas separate all roots
            t = trace_from_orbit(orbit, betas[depth], g, gf)
            depth += 1
            res = split_with_sums(g, t, sums)
            if res == 'trivial':
                continue
            if res is not None:
                a, pa, other, po = res
                if len(other) == 1:
                    continue
                return rec(a, depth, pa) + rec(other, depth, po)
            a = poly_gcd(t, g, gf)
            if len(a) <= 1 or len(a) == len(g):
                continue
            a = poly_monic(a, gf)
            other = poly_divmod(g, a, gf)[0]
            return rec(a, depth, None) + rec(other, depth, None)

    return sorted(rec(f, 0, q_sums))

def decode_fast(odd, max_count, gf, basis=1, mont=True):
    """minisketch decode: bm_odd, reverse, find_roots_fast with the sketch as
    the root node's power sums. Same result as Minisketch.decode (C46)."""
    conn = bm_odd(odd, max_count, gf)
    if conn is None:
        return None
    if len(conn) == 1:
        return []
    if conn[-1] == 0:
        return None
    f = list(reversed(conn))
    n = len(f) - 1
    q_sums = None
    if 2 * len(odd) >= n + 1:
        q_sums = [0] * (2 * len(odd) + 1)
        for i, v in enumerate(odd):
            q_sums[2 * i + 1] = v
        for k in range(2, 2 * len(odd) + 1, 2):
            q_sums[k] = gf.sqr(q_sums[k // 2])
    return find_roots_fast(f, gf, basis, q_sums, mont)

class TestOracle(unittest.TestCase):
    """Each optimized-decode oracle function against its textbook counterpart."""

    def random_set(self, gf, n):
        n = min(n, (1 << gf.field_size) - 1)
        out = set()
        while len(out) < n:
            out.add(random.randrange(1, 1 << gf.field_size))
        return list(out)

    def odd_syndromes(self, gf, elements, capacity):
        sketch = Minisketch(gf.field_size, capacity)
        for e in elements:
            sketch.add(e)
        return sketch.odd_syndromes

    def test_bm_odd_random(self):
        for b in list(range(2, 17)) + [32, 64]:
            gf = make_gf(b)
            for _ in range(6):
                cap = random.randrange(1, 40 if b <= 16 else 12)
                kind = random.randrange(3)
                if kind == 0:  # a sparse set, so early termination fires
                    odd = self.odd_syndromes(gf, self.random_set(gf, random.randrange(cap // 3 + 1)), cap)
                elif kind == 1:  # a sparse set with one syndrome perturbed
                    odd = self.odd_syndromes(gf, self.random_set(gf, random.randrange(cap // 3 + 1)), cap)
                    odd[random.randrange(cap)] ^= random.randrange(1, 1 << b)
                else:  # garbage
                    odd = [random.randrange(1 << b) for _ in range(cap)]
                max_degree = random.randrange(cap + 1)
                self.assertEqual(bm_odd(odd, max_degree, gf), bm_reference(odd, max_degree, gf))

    def test_bm_odd_exhaustive(self):
        """Every odd-syndrome vector and every max_degree, GF(4)..GF(256)."""
        for b, cap in ((2, 5), (3, 4), (4, 3), (5, 2), (6, 2), (8, 2)):
            gf = make_gf(b)
            for code in range(1 << (b * cap)):
                odd = [(code >> (b * i)) & ((1 << b) - 1) for i in range(cap)]
                full = bm_reference(odd, cap, gf)
                for max_degree in range(cap + 1):
                    expect = full if full is not None and len(full) - 1 <= max_degree else None
                    self.assertEqual(bm_odd(odd, max_degree, gf), expect)

    def test_middle_product(self):
        gf = make_gf(16)
        for n in list(range(1, 20)) + [32, 48, 64]:
            c = [random.randrange(1 << 16) for _ in range(n)]
            w = [random.randrange(1 << 16) for _ in range(2 * n - 1)]
            self.assertEqual(middle_product_kara(c, w, gf), middle_product_naive(c, w, gf))

    def test_orbit_and_betas(self):
        for b in (5, 8, 12, 16, 24, 32, 64):
            gf = make_gf(b)
            betas = tower_betas(random.randrange(1, 1 << b), gf)
            self.assertEqual(gf2_rank(betas), b)
            roots = self.random_set(gf, 12)
            f = poly_from_roots(roots, gf)
            orbit = frobenius_orbit(f, gf)
            self.assertEqual(orbit[b], [0, 1])
            g = poly_from_roots(roots[:5], gf)
            for beta in betas[:6]:
                t = trace_from_orbit(orbit, beta, g, gf)
                self.assertEqual(t, poly_tracemod(g, beta, gf))
                for r in roots[:5]:
                    self.assertEqual(poly_eval(t, r, gf), elem_trace(gf.mul(beta, r), gf))
            # Not split: an irreducible quadratic factor, and a double root.
            c = next(c for c in iter(lambda: random.randrange(1, 1 << b), None) if elem_trace(c, gf) == 1)
            self.assertNotEqual(frobenius_orbit(poly_mul(f, [c, 1, 1], gf), gf)[b], [0, 1])
            self.assertNotEqual(frobenius_orbit(poly_mul(f, [roots[0], 1], gf), gf)[b], [0, 1])

    def test_additive_fft(self):
        for b in (4, 8, 16):
            gf = make_gf(b)
            beta = cantor_basis(gf)
            self.assertEqual(beta[0], 1)
            for j in range(1, b):
                self.assertEqual(gf.sqr(beta[j]) ^ beta[j], beta[j - 1])
            self.assertEqual(gf2_rank(beta), b)
            for m in range(0, min(b, 7)):
                eta = random.randrange(1 << b)
                pts = coset_points(beta, eta, m)
                self.assertEqual(len(set(pts)), 1 << m)
                # The coset's vanishing polynomial is S^m + S^m(eta), derivative 1 (C29).
                z = poly_add(subspace_poly(m), [subspace_eval(m, eta, gf)])
                self.assertEqual(z, poly_from_roots(pts, gf))
                self.assertEqual(poly_deriv(z), [1])
                f = poly_trim([random.randrange(1 << b) for _ in range(1 << m)])
                vals = additive_fft(f, beta, eta, m, gf)
                self.assertEqual(vals, [poly_eval(f, w, gf) for w in pts])
                self.assertEqual(additive_ifft(vals, beta, eta, gf), f)
                d = to_novel(f, m)
                self.assertEqual(from_novel(d, gf), f)
                self.assertEqual(from_novel(novel_derivative(d), gf), poly_deriv(f))

    def coset_without_roots(self, gf, beta, f, m):
        for blk in (1, 2, 3, 0):
            eta = (beta[m] if blk & 1 else 0) ^ (beta[m + 1] if blk & 2 else 0)
            if all(poly_eval(f, w, gf) for w in coset_points(beta, eta, m)):
                return eta
        return None

    def test_mont_square(self):
        gf = make_gf(16)
        beta = cantor_basis(gf)
        for d in (3, 8, 31, 64, 100):
            f = poly_from_roots(self.random_set(gf, d), gf)
            m = (d - 1).bit_length()
            eta = self.coset_without_roots(gf, beta, f, m)
            pts = coset_points(beta, eta, m)
            z = poly_from_roots(pts, gf)
            r = poly_trim([random.randrange(1 << 16) for _ in range(1 << m)])
            rv = additive_fft(r, beta, eta, m, gf)
            fv = [poly_eval(f, w, gf) for w in pts]
            fdv = [poly_eval(poly_deriv(f), w, gf) for w in pts]
            r_new = additive_ifft(mont_square(rv, fv, fdv, beta, eta, gf), beta, eta, gf)
            self.assertLess(len(r_new), (1 << m) + 1)
            # r_new * Z = r^2 (mod f).
            self.assertEqual(poly_divmod(poly_add(poly_mul(r_new, z, gf), poly_sqr(r, gf)), f, gf)[1], [])

    def test_mont_carry_to_factor(self):
        gf = make_gf(16)
        beta = cantor_basis(gf)
        for m in (2, 4, 6):
            n = 1 << m
            eta = random.randrange(1 << 16)
            pts = coset_points(beta, eta, m)
            while True:
                g = poly_from_roots(self.random_set(gf, n // 2), gf)
                if all(poly_eval(g, w, gf) for w in pts[:n // 2]):
                    break
            v = poly_trim([random.randrange(1 << 16) for _ in range(n)])
            vals1 = mont_carry_to_factor(additive_fft(v, beta, eta, m, gf), g, beta, eta, gf)
            r = additive_ifft(vals1, beta, eta ^ beta[m - 1], gf)
            z0 = poly_from_roots(pts[:n // 2], gf)
            self.assertEqual(poly_divmod(poly_add(poly_mul(r, z0, gf), v), g, gf)[1], [])

    def test_montgomery_orbit(self):
        """Roots planted off the first candidate coset, so the Montgomery orbit
        runs; it must equal the textbook orbit."""
        gf = make_gf(16)
        beta = cantor_basis(gf)
        for d in (3, 17, 100, 300):
            m = max(1, (d - 1).bit_length())
            avoid = set(coset_points(beta, beta[m], m))
            roots = [r for r in self.random_set(gf, d + len(avoid)) if r not in avoid][:d]
            f = poly_from_roots(roots, gf)
            orbit, eta, m2 = montgomery_orbit(f, beta, gf)
            self.assertEqual((eta, m2), (beta[m], m))
            self.assertEqual(orbit, frobenius_orbit(f, gf))

    def test_coset_divide(self):
        gf = make_gf(16)
        beta = cantor_basis(gf)
        for d in (4, 20, 64):
            roots = self.random_set(gf, d)
            g = poly_from_roots(roots, gf)
            m = (d - 1).bit_length()
            eta = self.coset_without_roots(gf, beta, g, m)
            h = poly_from_roots(roots[:d // 3], gf)
            self.assertEqual(coset_divide(g, h, beta, eta, m, gf), poly_divmod(g, h, gf)[0])
            # h not a divisor: one root outside g.
            bad = poly_from_roots(roots[:d // 3 - 1] + [next(x for x in range(1, 1 << 16) if x not in roots)], gf)
            if all(poly_eval(bad, w, gf) for w in coset_points(beta, eta, m)):
                self.assertIsNone(coset_divide(g, bad, beta, eta, m, gf))

    def test_power_sum_split(self):
        for b in (8, 16):
            gf = make_gf(b)
            for n in (2, 5, 17, 40):
                roots = self.random_set(gf, n)
                g = poly_from_roots(roots, gf)
                beta = random.randrange(1, 1 << b)
                t = poly_tracemod(g, beta, gf)
                side = [r for r in roots if elem_trace(gf.mul(beta, r), gf) == 1]
                q_sums = power_sums(roots, 3 * n + 2, gf)
                P, conn = power_sum_split(q_sums, t, n, gf)
                self.assertEqual(P[1:], power_sums(side, len(P) - 1, gf)[1:])
                # BM recovers the side when it holds at most half of the sums'
                # worth of roots; otherwise a result that divides g is still
                # the side (C35), and anything else fails CosetDivide.
                loc = list(reversed(conn)) if conn is not None and len(conn) > 1 else [1]
                if 2 * len(side) <= len(P) - 1:
                    self.assertEqual(loc, poly_from_roots(side, gf) if side else [1])
                elif conn is not None and not poly_divmod(g, loc, gf)[1]:
                    self.assertEqual(loc, poly_from_roots(side, gf))

    def test_find_roots_fast(self):
        for b, sizes in ((4, (0, 1, 2, 3, 9, 15)), (8, (2, 3, 7, 40, 200)), (16, (2, 9, 64, 300)),
                         (32, (2, 7, 20)), (64, (3, 12))):
            gf = make_gf(b)
            for n in sizes:
                roots = self.random_set(gf, n)
                f = poly_from_roots(roots, gf)
                sums = power_sums(roots, 3 * n + 2, gf)
                for mont in (True, False):
                    self.assertEqual(find_roots_fast(f, gf, random.randrange(1, 1 << b), sums, mont), sorted(roots))
                if n >= 2:
                    c = next(c for c in iter(lambda: random.randrange(1, 1 << b), None) if elem_trace(c, gf) == 1)
                    self.assertIsNone(find_roots_fast(poly_mul(f, [c, 1, 1], gf), gf))
                    self.assertIsNone(find_roots_fast(poly_mul(f, [roots[0], 1], gf), gf))

    def test_decode_fast_montgomery(self):
        """A full sketch decode down the Montgomery / power-sum / CosetDivide
        path (roots off the root coset), checking those paths actually ran."""
        global coset_divide
        gf = make_gf(16)
        beta = cantor_basis(gf)
        accepted = []
        orig = coset_divide
        def counting(*args):
            r = orig(*args)
            accepted.append(r is not None)
            return r
        coset_divide = counting
        try:
            for n in (40, 300):
                m = max(1, (n - 1).bit_length())
                avoid = set(coset_points(beta, beta[m], m))
                roots = [r for r in self.random_set(gf, n + len(avoid)) if r not in avoid][:n]
                cap = n + 8
                odd = self.odd_syndromes(gf, roots, cap)
                self.assertEqual(decode_fast(odd, cap, gf), sorted(roots))
        finally:
            coset_divide = orig
        self.assertTrue(any(accepted))

    def test_decode_fast(self):
        for b in list(range(2, 17)) + [32]:
            gf = make_gf(b)
            for _ in range(8):
                cap = random.randrange(1, 30 if b <= 16 else 8)
                n = random.randrange(2 * cap + 2)
                odd = self.odd_syndromes(gf, self.random_set(gf, n), cap)
                if random.randrange(4) == 0:
                    odd = [random.randrange(1 << b) for _ in range(cap)]
                sketch = Minisketch(b, cap)
                sketch.odd_syndromes = list(odd)
                max_count = random.randrange(cap + 1)
                self.assertEqual(decode_fast(odd, max_count, gf, random.randrange(1, 1 << b)),
                                 sketch.decode(max_count))

    def test_decode_fast_exhaustive(self):
        """Every sketch of the smallest fields and capacities (C11, C46)."""
        for b, cap in ((2, 3), (3, 3), (4, 3), (5, 2), (8, 2)):
            gf = make_gf(b)
            decodable = {}
            for code in range(1 << (b * cap)):
                odd = [(code >> (b * i)) & ((1 << b) - 1) for i in range(cap)]
                sketch = Minisketch(b, cap)
                sketch.odd_syndromes = list(odd)
                got = decode_fast(odd, cap, gf)
                self.assertEqual(got, sketch.decode(cap))
                if got is not None:
                    self.assertNotIn(0, got)
                    self.assertEqual(self.odd_syndromes(gf, got, cap), odd)  # C37: re-encodes
                    decodable[len(got)] = decodable.get(len(got), 0) + 1
            for i in range(cap + 1):
                self.assertEqual(decodable.get(i, 0), math.comb((1 << b) - 1, i))

def base_fpbits(bits, capacity):
    """floor(log2(2^(bits*capacity) / sum_{k<=capacity} C(2^bits - 1, k))), exactly.

    The number of false-positive bits a sketch with max_elements == capacity
    provides: of the 2^(bits*capacity) sketches, the denominator counts those
    that decode (one per set of at most capacity nonzero elements). The C++
    BaseFPBits approximates this, exactly while the result is below 256."""
    count = sum(math.comb((1 << bits) - 1, k) for k in range(capacity + 1))
    return bits * capacity - (count - 1).bit_length()

def compute_capacity(bits, max_elements, fpbits):
    """minisketch_compute_capacity, from the exact base_fpbits."""
    if bits == 0:
        return 0
    base = base_fpbits(bits, max_elements)
    if base >= fpbits:
        return max_elements
    return max_elements + (fpbits - base + bits - 1) // bits

def compute_max_elements(bits, capacity, fpbits):
    """minisketch_compute_max_elements: the largest max_elements whose
    compute_capacity fits in capacity, or 0 if there is none."""
    if bits == 0:
        return 0
    for max_elements in range(capacity, -1, -1):
        if compute_capacity(bits, max_elements, fpbits) <= capacity:
            return max_elements
    return 0

class TestFalsePositives(unittest.TestCase):
    """Sanity checks of the exact capacity functions."""

    def test_base_fpbits(self):
        self.assertEqual(base_fpbits(8, 0), 0)
        # One element in GF(2^b): 2^b sketches, 2^b - 1 + 1 decode.
        for bits in range(2, 65):
            self.assertEqual(base_fpbits(bits, 1), 0)
        # A full GF(4) sketch of capacity 3 decodes every 3-subset of {1,2,3}.
        self.assertEqual(base_fpbits(2, 3), 6 - 3)

    def test_inverse(self):
        for bits in (2, 5, 12, 32, 64):
            for fpbits in (0, 8, 32):
                for max_elements in range(0, 40):
                    capacity = compute_capacity(bits, max_elements, fpbits)
                    self.assertGreaterEqual(compute_max_elements(bits, capacity, fpbits), max_elements)

class TestGoldenVectors(unittest.TestCase):
    """Fixed known-good serializations, shared with the C++ test suite.

    Guards against silent serialization-format drift: the same vectors are
    checked against libminisketch by test.cpp and tests/differential.py.
    Regenerate (only if the format intentionally changes) with:
    differential.py --lib ... --generate golden_vectors.json
    """

    def test_golden_vectors(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden_vectors.json')
        with open(path, encoding='utf-8') as f:
            vectors = json.load(f)
        self.assertGreater(len(vectors), 0)
        for vec in vectors:
            sketch = Minisketch(vec['bits'], vec['capacity'])
            for v in vec['elements']:
                sketch.add(v)
            self.assertEqual(sketch.serialize().hex(), vec['serialized_hex'], vec)
            sketch2 = Minisketch(vec['bits'], vec['capacity'])
            sketch2.deserialize(bytes.fromhex(vec['serialized_hex']))
            self.assertEqual(sketch2.decode(), vec['elements'], vec)

if __name__ == '__main__':
    unittest.main()
