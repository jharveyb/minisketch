#!/usr/bin/env python3
# Copyright (c) 2026 The minisketch developers
# Distributed under the MIT software license, see the accompanying
# file LICENSE or http://www.opensource.org/licenses/mit-license.php.

"""Differential test: the C/C++ libminisketch vs. the pure-Python reimplementation.

Loads libminisketch through ctypes and cross-checks, on seeded random inputs
over small parameters, that both implementations produce byte-identical
serializations and identical decode results (in both directions), for every
implementation the library supports.

Usage: differential.py --lib path/to/libminisketch.so [--iters N] [--generate golden_vectors.json]

Reproduce a failure by setting MINISKETCH_TEST_SEED to the printed seed.
"""

import argparse
import ctypes
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pyminisketch import Minisketch as PyMinisketch  # noqa: E402


class CMinisketch:
    """Minimal ctypes wrapper around the libminisketch C API."""

    def __init__(self, lib, bits, implementation, capacity):
        self.lib = lib
        self.bits = bits
        self.capacity = capacity
        self.ptr = lib.minisketch_create(bits, implementation, capacity)
        assert self.ptr, f"minisketch_create({bits}, {implementation}, {capacity}) failed"

    def __del__(self):
        if getattr(self, 'ptr', None):
            self.lib.minisketch_destroy(self.ptr)

    def add(self, element):
        self.lib.minisketch_add_uint64(self.ptr, element)

    def serialized_size(self):
        return self.lib.minisketch_serialized_size(self.ptr)

    def serialize(self):
        buf = ctypes.create_string_buffer(self.serialized_size())
        self.lib.minisketch_serialize(self.ptr, buf)
        return buf.raw

    def deserialize(self, data):
        assert len(data) == self.serialized_size()
        self.lib.minisketch_deserialize(self.ptr, data)

    def merge(self, other):
        return self.lib.minisketch_merge(self.ptr, other.ptr)

    def decode(self, max_elements):
        out = (ctypes.c_uint64 * max(max_elements, 1))()
        ret = self.lib.minisketch_decode(self.ptr, max_elements, out)
        if ret < 0:
            return None
        return sorted(out[i] for i in range(ret))


def load_lib(path):
    lib = ctypes.CDLL(path)
    lib.minisketch_create.restype = ctypes.c_void_p
    lib.minisketch_create.argtypes = [ctypes.c_uint32, ctypes.c_uint32, ctypes.c_size_t]
    lib.minisketch_destroy.argtypes = [ctypes.c_void_p]
    lib.minisketch_bits_supported.argtypes = [ctypes.c_uint32]
    lib.minisketch_implementation_max.restype = ctypes.c_uint32
    lib.minisketch_implementation_supported.argtypes = [ctypes.c_uint32, ctypes.c_uint32]
    lib.minisketch_add_uint64.argtypes = [ctypes.c_void_p, ctypes.c_uint64]
    lib.minisketch_serialized_size.restype = ctypes.c_size_t
    lib.minisketch_serialized_size.argtypes = [ctypes.c_void_p]
    lib.minisketch_serialize.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    lib.minisketch_deserialize.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    lib.minisketch_merge.restype = ctypes.c_size_t
    lib.minisketch_merge.argtypes = [ctypes.c_void_p, ctypes.c_void_p]
    lib.minisketch_decode.restype = ctypes.c_ssize_t
    lib.minisketch_decode.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.POINTER(ctypes.c_uint64)]
    return lib


def implementations(lib, bits):
    return [impl for impl in range(lib.minisketch_implementation_max() + 1)
            if lib.minisketch_implementation_supported(bits, impl)]


def run_case(lib, bits, capacity, rng):
    """One random differential case for a given field size and capacity."""
    max_elem = (1 << bits) - 1
    count = rng.randrange(0, min(2 * capacity, max_elem) + 1)
    elements = [rng.randrange(1, max_elem + 1) for _ in range(count)]
    effective = sorted(v for v in set(elements) if elements.count(v) % 2 == 1)

    py = PyMinisketch(bits, capacity)
    for v in elements:
        py.add(v)
    py_ser = py.serialize()
    py_decode = py.decode(capacity)
    if py_decode is not None:
        py_decode = sorted(py_decode)

    for impl in implementations(lib, bits):
        c = CMinisketch(lib, bits, impl, capacity)
        for v in elements:
            c.add(v)
        # Byte-identical serialization.
        c_ser = c.serialize()
        assert c_ser == py_ser, (
            f"serialize mismatch bits={bits} impl={impl} capacity={capacity} "
            f"elements={elements}: C={c_ser.hex()} py={py_ser.hex()}")
        # Cross-deserialization: C decodes the Python serialization and vice versa.
        c2 = CMinisketch(lib, bits, impl, capacity)
        c2.deserialize(py_ser)
        assert c2.serialize() == py_ser
        c_decode = c2.decode(capacity)
        assert c_decode == py_decode, (
            f"decode mismatch bits={bits} impl={impl} capacity={capacity} "
            f"elements={elements}: C={c_decode} py={py_decode}")
        if len(effective) <= capacity:
            assert c_decode == effective
        # Python decodes the C serialization.
        py2 = PyMinisketch(bits, capacity)
        py2.deserialize(c_ser)
        py2_decode = py2.decode(capacity)
        assert (py2_decode if py2_decode is None else sorted(py2_decode)) == py_decode

        # Merge behaves identically, and equals the XOR of the parts' serializations.
        split = rng.randrange(0, len(elements) + 1)
        c_a = CMinisketch(lib, bits, impl, capacity)
        c_b = CMinisketch(lib, bits, impl, capacity)
        for v in elements[:split]:
            c_a.add(v)
        for v in elements[split:]:
            c_b.add(v)
        ser_a, ser_b = c_a.serialize(), c_b.serialize()
        c_a.merge(c_b)
        merged_ser = c_a.serialize()
        assert merged_ser == py_ser, (
            f"merge mismatch bits={bits} impl={impl} capacity={capacity} split={split}")
        assert bytes(x ^ y for x, y in zip(ser_a, ser_b)) == merged_ser


def generate_golden(lib, path):
    """Regenerate the golden vectors file from the C library (implementation 0)."""
    rng = random.Random(0x6d696e69736b6574)  # Fixed: golden vectors never change.
    vectors = []
    for bits in [2, 3, 8, 11, 16, 32, 64]:
        for capacity in [1, 2, 4, 8]:
            max_elem = (1 << bits) - 1
            count = rng.randrange(0, min(capacity, max_elem) + 1)
            elements = []
            while len(elements) < count:
                v = rng.randrange(1, max_elem + 1)
                if v not in elements:
                    elements.append(v)
            c = CMinisketch(lib, bits, 0, capacity)
            for v in elements:
                c.add(v)
            vectors.append({
                "bits": bits,
                "capacity": capacity,
                "elements": sorted(elements),
                "serialized_hex": c.serialize().hex(),
            })
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(vectors, f, indent=1)
        f.write('\n')
    print(f"Wrote {len(vectors)} golden vectors to {path}")


def check_golden(lib, path):
    """The C library must reproduce all golden vectors."""
    with open(path, encoding='utf-8') as f:
        vectors = json.load(f)
    for vec in vectors:
        c = CMinisketch(lib, vec["bits"], 0, vec["capacity"])
        for v in vec["elements"]:
            c.add(v)
        assert c.serialize().hex() == vec["serialized_hex"], f"golden vector mismatch: {vec}"
    print(f"{len(vectors)} golden vectors OK (C library)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--lib', required=True, help="path to libminisketch shared library")
    parser.add_argument('--iters', type=int, default=300, help="number of random cases")
    parser.add_argument('--generate', metavar='PATH', help="regenerate golden vectors and exit")
    args = parser.parse_args()

    lib = load_lib(args.lib)
    if args.generate:
        generate_golden(lib, args.generate)
        return

    seed = os.environ.get('MINISKETCH_TEST_SEED')
    seed = int(seed, 0) if seed else random.randrange(1 << 63)
    print(f"Test seed: {seed:#x} (reproduce with MINISKETCH_TEST_SEED={seed:#x})")
    rng = random.Random(seed)

    golden = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'golden_vectors.json')
    if os.path.exists(golden):
        check_golden(lib, golden)

    # Small parameters decode fast in pure Python; the C++ side is fully
    # random-tested elsewhere, this is about implementation agreement.
    for i in range(args.iters):
        bits = rng.randrange(2, 17)
        capacity = rng.randrange(0, 17)
        run_case(lib, bits, capacity, rng)
    print(f"{args.iters} differential cases OK")


if __name__ == '__main__':
    main()
