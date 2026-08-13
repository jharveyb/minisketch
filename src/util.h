/**********************************************************************
 * Copyright (c) 2018 Pieter Wuille, Greg Maxwell, Gleb Naumenko      *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_UTIL_H_
#define _MINISKETCH_UTIL_H_

#ifdef MINISKETCH_VERIFY
#include <stdio.h>
#endif

#if !defined(__GNUC_PREREQ)
# if defined(__GNUC__)&&defined(__GNUC_MINOR__)
#  define __GNUC_PREREQ(_maj,_min) \
 ((__GNUC__<<16)+__GNUC_MINOR__>=((_maj)<<16)+(_min))
# else
#  define __GNUC_PREREQ(_maj,_min) 0
# endif
#endif

#if __GNUC_PREREQ(3, 0)
#define EXPECT(x,c) __builtin_expect((x),(c))
#else
#define EXPECT(x,c) (x)
#endif

/* Decode-stage markers for the Coz causal profiler (doc/testing.md,
 * "Causal profiling"). In MINISKETCH_COZ builds these are COZ_BEGIN/COZ_END
 * latency progress points, so profiles report per-stage impact directly; in
 * normal builds they are no-ops and coz.h is not included. */
#ifdef MINISKETCH_COZ
#include <coz.h>
#define MINISKETCH_COZ_BEGIN(name) COZ_BEGIN(name)
#define MINISKETCH_COZ_END(name) COZ_END(name)
#else
#define MINISKETCH_COZ_BEGIN(name)
#define MINISKETCH_COZ_END(name)
#endif

/* Keep large cold functions out of their callers (instantiating them inline
 * inside a hot recursion measurably hurts its code locality). */
#if defined(__GNUC__)
#define NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define NOINLINE __declspec(noinline)
#else
#define NOINLINE
#endif

/* Assertion macros */

/**
 * Unconditional failure on condition failure.
 * Primarily used in testing harnesses.
 */
#define CHECK(cond) do { \
    if (EXPECT(!(cond), 0)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Check condition failed: " #cond); \
        abort(); \
    } \
} while(0)

/**
 * Check macro that does nothing in normal non-verify builds but crashes in verify builds.
 * This is used to test conditions at runtime that should always be true, but are either
 * expensive to test or in locations where returning on failure would be messy.
 */
#ifdef MINISKETCH_VERIFY
#define CHECK_SAFE(cond) CHECK(cond)
#else
#define CHECK_SAFE(cond)
#endif

/**
 * Check a condition and return on failure in non-verify builds, crash in verify builds.
 * Used for inexpensive conditions which believed to be always true in locations where
 * a graceful exit is possible.
 */
#ifdef MINISKETCH_VERIFY
#define CHECK_RETURN(cond, rvar) do { \
    if (EXPECT(!(cond), 0)) { \
        fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, "Check condition failed: " #cond); \
        abort(); \
        return rvar; /* Does nothing, but causes compile to warn on incorrect return types. */ \
    } \
} while(0)
#else
#define CHECK_RETURN(cond, rvar) do { \
    if (EXPECT(!(cond), 0)) { \
        return rvar; \
    } \
} while(0)
#endif

#endif
