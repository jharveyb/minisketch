/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Unit tests for internal algorithm stages, instantiated with the generic
 * (table-based) field implementations. See unit_tests_impl.h for the actual
 * properties; unit_tests_clmul.cpp runs the same properties against the clmul
 * field implementations.
 *
 * Run with --seed=N (after --) or MINISKETCH_TEST_SEED=N to reproduce a
 * failure, e.g.: unit-tests -- --seed=0x1234
 */

#define BOOST_TEST_MODULE minisketch unit tests
#include <boost/test/unit_test.hpp>

// Bring in the generic field implementations (Field2..Field8, Field11, Field16,
// Field27, Field32, Field64, ... in this TU's anonymous namespace). This also
// pulls in sketch_impl.h. Do not link this binary against libminisketch: these
// includes already define the ConstructGeneric* symbols.
#include "fields/generic_1byte.cpp"
#include "fields/generic_2bytes.cpp"
#include "fields/generic_4bytes.cpp"
#include "fields/generic_8bytes.cpp"

#include "test_utils.h"
#include "unit_tests_impl.h"

uint64_t g_test_seed = 0;

namespace {

struct SeedFixture {
    SeedFixture() {
        auto& suite = boost::unit_test::framework::master_test_suite();
        g_test_seed = GetTestSeed(suite.argc, suite.argv);
    }
};

} // namespace

BOOST_TEST_GLOBAL_FIXTURE(SeedFixture);

/* Small fields: field ops are checked exhaustively over all element pairs. */
MINISKETCH_FIELD_TEST(generic_field_2, Field2)
MINISKETCH_FIELD_TEST(generic_field_3, Field3)
MINISKETCH_FIELD_TEST(generic_field_4, Field4)
MINISKETCH_FIELD_TEST(generic_field_5, Field5)
MINISKETCH_FIELD_TEST(generic_field_6, Field6)
MINISKETCH_FIELD_TEST(generic_field_7, Field7)
MINISKETCH_FIELD_TEST(generic_field_8, Field8)

/* Representative larger fields, including the odd-size (non-byte-aligned)
 * 11- and 27-bit fields and the widest 64-bit field. */
MINISKETCH_FIELD_TEST(generic_field_11, Field11)
MINISKETCH_FIELD_TEST(generic_field_16, Field16)
MINISKETCH_FIELD_TEST(generic_field_27, Field27)
MINISKETCH_FIELD_TEST(generic_field_32, Field32)
MINISKETCH_FIELD_TEST(generic_field_64, Field64)

/* Full decode pipeline at capacity 1024 with 64-bit elements (the intended
 * practical configuration); ~1s, by far the largest case in this suite. */
BOOST_AUTO_TEST_CASE(generic_field_64_capacity_1024) {
    Field64 field;
    TestRand rng(g_test_seed + 64001);
    ut::TestDecodeStagesLarge(field, rng, 1024);
}

/* BitWriter/BitReader roundtrip with a mix of write widths crossing byte
 * boundaries; the field Serialize/Deserialize roundtrips in the per-field
 * cases cover the remaining widths. */
BOOST_AUTO_TEST_CASE(bitrw_roundtrip) {
    TestRand rng(g_test_seed);
    for (int rep = 0; rep < 1000; ++rep) {
        uint64_t v1 = rng.RandBits(1), v2 = rng.RandBits(2), v7 = rng.RandBits(7);
        uint64_t v8 = rng.RandBits(8), v9 = rng.RandBits(9), v13 = rng.RandBits(13);
        uint64_t v31 = rng.RandBits(31), v32 = rng.RandBits(32);
        uint64_t v63 = rng.RandBits(63), v64 = rng.RandBits(64);
        std::vector<unsigned char> buf(40, 0xa5);
        BitWriter writer(buf.data());
        writer.Write<1, uint64_t>(v1);
        writer.Write<2, uint64_t>(v2);
        writer.Write<7, uint64_t>(v7);
        writer.Write<8, uint64_t>(v8);
        writer.Write<9, uint64_t>(v9);
        writer.Write<13, uint64_t>(v13);
        writer.Write<31, uint64_t>(v31);
        writer.Write<32, uint64_t>(v32);
        writer.Write<63, uint64_t>(v63);
        writer.Write<64, uint64_t>(v64);
        writer.Flush();
        BitReader reader(buf.data());
        UT_REQUIRE((reader.Read<1, uint64_t>()) == v1);
        UT_REQUIRE((reader.Read<2, uint64_t>()) == v2);
        UT_REQUIRE((reader.Read<7, uint64_t>()) == v7);
        UT_REQUIRE((reader.Read<8, uint64_t>()) == v8);
        UT_REQUIRE((reader.Read<9, uint64_t>()) == v9);
        UT_REQUIRE((reader.Read<13, uint64_t>()) == v13);
        UT_REQUIRE((reader.Read<31, uint64_t>()) == v31);
        UT_REQUIRE((reader.Read<32, uint64_t>()) == v32);
        UT_REQUIRE((reader.Read<63, uint64_t>()) == v63);
        UT_REQUIRE((reader.Read<64, uint64_t>()) == v64);
    }
}
