/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_FUZZ_FUZZ_H_
#define _MINISKETCH_FUZZ_FUZZ_H_

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstddef>

/** All fuzz targets are compiled into one binary (like Bitcoin Core's fuzz
 * harness); the target to run is selected with the FUZZ environment variable:
 *
 *   FUZZ=decode ./fuzz -max_len=512 corpus/decode
 *
 * Define a target with:
 *
 *   FUZZ_TARGET(name) { ... use `data`, `size` ... }
 */

typedef void (*FuzzTargetFn)(const uint8_t* data, size_t size);

/** Register a fuzz target under a name; returns true (for static-init use). */
bool RegisterFuzzTarget(const char* name, FuzzTargetFn fn);

#define FUZZ_TARGET(name) \
    static void name##_fuzz_target(const uint8_t* data, size_t size); \
    static const bool name##_registered [[maybe_unused]] = RegisterFuzzTarget(#name, name##_fuzz_target); \
    static void name##_fuzz_target(const uint8_t* data, size_t size)

/** Assertion for fuzz targets: aborts (which libFuzzer treats as a crash and
 *  minimizes the input for) with a message on failure. */
#define FUZZ_CHECK(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "%s:%d: fuzz check failed: %s\n", __FILE__, __LINE__, #cond); \
        abort(); \
    } \
} while(0)

#endif
