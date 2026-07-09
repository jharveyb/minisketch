/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Unit tests for internal algorithm stages, instantiated with the clmul
 * (carryless multiplication) field implementations. Must be a separate TU from
 * unit_tests.cpp because the generic and clmul headers both define a `Field`
 * class in the anonymous namespace. Both are validated against the same
 * independent reference arithmetic, which also makes them transitively
 * consistent with each other. */

#ifdef HAVE_CLMUL

#include <boost/test/unit_test.hpp>

#ifdef _MSC_VER
#include <intrin.h>
#else
#include <cpuid.h>
#endif

#include "fields/clmul_1byte.cpp"
#include "fields/clmul_2bytes.cpp"
#include "fields/clmul_4bytes.cpp"
#include "fields/clmul_8bytes.cpp"

#include "test_utils.h"
#include "unit_tests_impl.h"

namespace {

/* Same runtime detection as minisketch.cpp's EnableClmul(). */
bool ClmulSupported() {
#ifdef _MSC_VER
    int regs[4];
    __cpuid(regs, 1);
    return (regs[2] & 0x2);
#else
    uint32_t eax, ebx, ecx, edx;
    return (__get_cpuid(1, &eax, &ebx, &ecx, &edx) && (ecx & 0x2));
#endif
}

} // namespace

#define MINISKETCH_CLMUL_FIELD_TEST(name, FieldType) \
    BOOST_AUTO_TEST_CASE(name) { \
        if (!ClmulSupported()) { \
            BOOST_TEST_MESSAGE("CPU lacks CLMUL support; skipping " #name); \
            return; \
        } \
        ut::RunAllFieldTests<FieldType>(__LINE__); \
    }

/* Small fields (FieldTri* are the trinomial-modulus specializations). */
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_tri_2, FieldTri2)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_tri_3, FieldTri3)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_5, Field5)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_tri_5, FieldTri5)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_8, Field8)

/* Representative larger fields, matching the generic TU's selection. */
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_11, Field11)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_tri_11, FieldTri11)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_16, Field16)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_27, Field27)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_32, Field32)
MINISKETCH_CLMUL_FIELD_TEST(clmul_field_64, Field64)

/* Full decode pipeline at capacity 1024 with 64-bit elements (the intended
 * practical configuration). */
BOOST_AUTO_TEST_CASE(clmul_field_64_capacity_1024) {
    if (!ClmulSupported()) {
        BOOST_TEST_MESSAGE("CPU lacks CLMUL support; skipping");
        return;
    }
    Field64 field;
    TestRand rng(g_test_seed + 64002);
    ut::TestDecodeStagesLarge(field, rng, 1024);
}

#endif // HAVE_CLMUL
