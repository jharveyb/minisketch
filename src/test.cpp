/**********************************************************************
 * Copyright (c) 2018,2021 Pieter Wuille, Greg Maxwell, Gleb Naumenko *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#include <algorithm>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../include/minisketch.h"
#include "test_utils.h"
#include "util.h"

namespace {

uint64_t Combination(uint64_t n, uint64_t k) {
    if (n - k < k) k = n - k;
    uint64_t ret = 1;
    for (uint64_t i = 1; i <= k; ++i) {
        ret = (ret * n) / i;
        --n;
    }
    return ret;
}

/** Create a vector with Minisketch objects, one for each implementation. */
std::vector<Minisketch> CreateSketches(uint32_t bits, size_t capacity) {
    if (!Minisketch::BitsSupported(bits)) return {};
    std::vector<Minisketch> ret;
    for (uint32_t impl = 0; impl <= Minisketch::MaxImplementation(); ++impl) {
        if (Minisketch::ImplementationSupported(bits, impl)) {
            CHECK(Minisketch::BitsSupported(bits));
            ret.push_back(Minisketch(bits, impl, capacity));
            CHECK((bool)ret.back());
        } else {
            // implementation 0 must always work unless field size is disabled
            CHECK(impl != 0 || !Minisketch::BitsSupported(bits));
        }
    }
    return ret;
}

/** Test properties by exhaustively decoding all 2**(bits*capacity) sketches
 *  with specified capacity and bits. */
void TestExhaustive(uint32_t bits, size_t capacity) {
    auto sketches = CreateSketches(bits, capacity);
    if (sketches.empty()) return;
    auto sketches_rebuild = CreateSketches(bits, capacity);

    std::vector<unsigned char> serialized;
    std::vector<unsigned char> serialized_empty;
    std::vector<uint64_t> counts; //!< counts[i] = number of results with i elements
    std::vector<uint64_t> elements_0; //!< Result vector for elements for impl=0
    std::vector<uint64_t> elements_other; //!< Result vector for elements for other impls
    std::vector<uint64_t> elements_too_small; //!< Result vector that's too small

    counts.resize(capacity + 1);
    serialized.resize(sketches[0].GetSerializedSize());
    serialized_empty.resize(sketches[0].GetSerializedSize());

    // Iterate over all (bits)-bit sketches with (capacity) syndromes.
    for (uint64_t x = 0; (x >> (bits * capacity)) == 0; ++x) {
        // Construct the serialization.
        for (size_t i = 0; i < serialized.size(); ++i) {
            serialized[i] = (x >> (i * 8)) & 0xFF;
        }

        // Compute all the solutions
        sketches[0].Deserialize(serialized);
        elements_0.resize(64);
        bool decodable_0 = sketches[0].Decode(elements_0);
        std::sort(elements_0.begin(), elements_0.end());

        // Verify that decoding with other implementations agrees.
        for (size_t impl = 1; impl < sketches.size(); ++impl) {
            sketches[impl].Deserialize(serialized);
            elements_other.resize(64);
            bool decodable_other = sketches[impl].Decode(elements_other);
            CHECK(decodable_other == decodable_0);
            std::sort(elements_other.begin(), elements_other.end());
            CHECK(elements_other == elements_0);
        }

        // If there are solutions:
        if (decodable_0) {
            if (!elements_0.empty()) {
                // Decoding with limit one less than the number of elements should fail.
                elements_too_small.resize(elements_0.size() - 1);
                for (size_t impl = 0; impl < sketches.size(); ++impl) {
                    CHECK(!sketches[impl].Decode(elements_too_small));
                }
            }

            // Reconstruct the sketch from the solutions.
            for (size_t impl = 0; impl < sketches.size(); ++impl) {
                // Clear the sketch.
                sketches_rebuild[impl].Deserialize(serialized_empty);
                // Load all decoded elements into it.
                for (uint64_t elem : elements_0) {
                    CHECK(elem != 0);
                    CHECK(elem >> bits == 0);
                    sketches_rebuild[impl].Add(elem);
                }
                // Reserialize the result
                auto serialized_rebuild = sketches_rebuild[impl].Serialize();
                // Compare
                CHECK(serialized == serialized_rebuild);
                // Count it
                if (impl == 0 && elements_0.size() <= capacity) ++counts[elements_0.size()];
            }
        }
    }

    // Verify that the number of decodable sketches with given elements is expected.
    uint64_t mask = bits == 64 ? UINT64_MAX : (uint64_t{1} << bits) - 1;
    for (uint64_t i = 0; i <= capacity && (i & mask) == i; ++i) {
        CHECK(counts[i] == Combination(mask, i));
    }
}

/** Test properties of sketches with random elements put in. */
void TestRandomized(uint32_t bits, size_t max_capacity, size_t iter, TestRand& rnd) {
    const uint64_t max_capacity_bits = std::min<uint64_t>(std::numeric_limits<uint64_t>::max() >> (64 - bits), max_capacity);
    const uint64_t max_element = std::numeric_limits<uint64_t>::max() >> (64 - bits);

    std::vector<uint64_t> decode_0;
    std::vector<uint64_t> decode_other;
    std::vector<uint64_t> decode_temp;
    std::vector<uint64_t> elements;

    for (size_t i = 0; i < iter; ++i) {
        // Determine capacity, and construct Minisketch objects for all implementations.
        uint64_t capacity = rnd.RandIncl(0, max_capacity_bits);
        auto sketches = CreateSketches(bits, capacity);
        // Sanity checks
        if (sketches.empty()) return;
        for (size_t impl = 0; impl < sketches.size(); ++impl) {
            CHECK(sketches[impl].GetBits() == bits);
            CHECK(sketches[impl].GetCapacity() == capacity);
            CHECK(sketches[impl].GetSerializedSize() == sketches[0].GetSerializedSize());
        }
        // Determine the number of elements, and create a vector to store them in.
        size_t element_count = std::max<int64_t>(0, (int64_t)capacity + (int64_t)rnd.RandIncl(0, 6) - 3);
        elements.resize(element_count);
        // Add the elements to all sketches
        for (size_t j = 0; j < element_count; ++j) {
            uint64_t elem = rnd.RandIncl(1, max_element);
            CHECK(elem != 0);
            elements[j] = elem;
            for (auto& sketch : sketches) sketch.Add(elem);
        }
        // Remove pairs of duplicates in elements, as they cancel out.
        std::sort(elements.begin(), elements.end());
        size_t real_element_count = element_count;
        for (size_t pos = 0; pos + 1 < elements.size(); ++pos) {
            if (elements[pos] == elements[pos + 1]) {
                real_element_count -= 2;
                // Set both elements to 0; afterwards we will move these to the end.
                elements[pos] = 0;
                elements[pos + 1] = 0;
                ++pos;
            }
        }
        if (real_element_count < element_count) {
            // Move all introduced zeroes (masking duplicates) to the end.
            std::sort(elements.begin(), elements.end(), [](uint64_t a, uint64_t b) { return a != b && (b == 0 || (a != 0 && a < b)); });
            CHECK(elements[real_element_count] == 0);
            elements.resize(real_element_count);
        }
        // Create and compare serializations
        auto serialized_0 = sketches[0].Serialize();
        for (size_t impl = 1; impl < sketches.size(); ++impl) {
            auto serialized_other = sketches[impl].Serialize();
            CHECK(serialized_other == serialized_0);
        }
        // Deserialize and reserialize them
        for (size_t impl = 0; impl < sketches.size(); ++impl) {
            sketches[impl].Deserialize(serialized_0);
            auto reserialized = sketches[impl].Serialize();
            CHECK(reserialized == serialized_0);
        }
        // Decode with limit set to the capacity, and compare results
        decode_0.resize(capacity);
        bool decodable_0 = sketches[0].Decode(decode_0);
        std::sort(decode_0.begin(), decode_0.end());
        for (size_t impl = 1; impl < sketches.size(); ++impl) {
            decode_other.resize(capacity);
            bool decodable_other = sketches[impl].Decode(decode_other);
            CHECK(decodable_other == decodable_0);
            std::sort(decode_other.begin(), decode_other.end());
            CHECK(decode_other == decode_0);
        }
        // If the result is decodable, it should also be decodable with limit
        // set to the actual number of elements, and not with one less.
        if (decodable_0) {
            for (auto& sketch : sketches) {
                decode_temp.resize(decode_0.size());
                bool decodable = sketch.Decode(decode_temp);
                CHECK(decodable);
                std::sort(decode_temp.begin(), decode_temp.end());
                CHECK(decode_temp == decode_0);
                if (!decode_0.empty()) {
                    decode_temp.resize(decode_0.size() - 1);
                    decodable = sketch.Decode(decode_temp);
                    CHECK(!decodable);
                }
            }
        }
        // If the actual number of elements is not higher than the capacity, the
        // result should be decodable, and the result should match what we put in.
        if (real_element_count <= capacity) {
            CHECK(decodable_0);
            CHECK(decode_0 == elements);
        }
    }
}

/** Fixed known-good serializations, shared with the Python test suite
 *  (tests/golden_vectors.json, checked by pyminisketch.py and
 *  differential.py). Guards against silent serialization-format drift.
 *  Regenerate (only if the format intentionally changes) with:
 *  tests/differential.py --lib ... --generate tests/golden_vectors.json */
void TestGoldenVectors() {
    struct Vector {
        uint32_t bits;
        size_t capacity;
        std::vector<uint64_t> elements;
        std::vector<unsigned char> serialized;
    };
    const std::vector<Vector> vectors = {
        {2, 1, {2}, {0x02}},
        {2, 2, {1, 2}, {0x03}},
        {2, 4, {}, {0x00}},
        {2, 8, {1, 3}, {0xb2, 0x2c}},
        {3, 1, {}, {0x00}},
        {3, 2, {4, 5}, {0x19}},
        {3, 4, {}, {0x00, 0x00}},
        {3, 8, {1, 2, 3, 6}, {0x4e, 0x20, 0xc6}},
        {8, 1, {}, {0x00}},
        {8, 2, {}, {0x00, 0x00}},
        {8, 4, {31, 66, 150}, {0xcb, 0x32, 0xdf, 0x3f}},
        {8, 8, {80, 145, 154, 160, 196, 211}, {0xec, 0x17, 0x0d, 0xd2, 0x6d, 0xfd, 0x05, 0xab}},
        {11, 1, {}, {0x00, 0x00}},
        {11, 2, {1316, 1575}, {0x03, 0x63, 0x1e}},
        {11, 4, {}, {0x00, 0x00, 0x00, 0x00, 0x00, 0x00}},
        {11, 8, {797}, {0x1d, 0x1b, 0x9a, 0x37, 0xec, 0x27, 0x71, 0xe7, 0x09, 0xf6, 0xa0}},
        {16, 1, {53856}, {0x60, 0xd2}},
        {16, 2, {11834, 37734}, {0x5c, 0xbd, 0x0d, 0xa2}},
        {16, 4, {15964, 19098, 24898, 46967}, {0xf3, 0xa2, 0xe0, 0x0d, 0xfe, 0x23, 0xec, 0x42}},
        {16, 8, {1257, 3279, 17862, 30004, 49457, 58430, 64787}, {0xc8, 0xe0, 0x29, 0xd5, 0x5b, 0xfc, 0x22, 0xe3, 0xf2, 0xf4, 0x76, 0xda, 0x98, 0xc6, 0x63, 0x2d}},
        {32, 1, {}, {0x00, 0x00, 0x00, 0x00}},
        {32, 2, {876209144}, {0xf8, 0xe3, 0x39, 0x34, 0x15, 0xc7, 0xf1, 0xe2}},
        {32, 4, {34494274}, {0x42, 0x57, 0x0e, 0x02, 0xcf, 0xa3, 0x41, 0x5e, 0x9f, 0x66, 0xc1, 0x4e, 0xab, 0x70, 0x35, 0x82}},
        {32, 8, {}, {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}},
        {64, 1, {}, {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}},
        {64, 2, {6232832097835246727u, 6790778162878701223u}, {0x20, 0x3a, 0xc2, 0xb2, 0x07, 0xc7, 0x42, 0x08, 0x8a, 0x1c, 0xe3, 0xbd, 0xb1, 0x8a, 0x6a, 0xa7}},
        {64, 4, {4149655823320236875u, 9344421695442453967u, 13600746181392431187u}, {0xd7, 0x3e, 0x4a, 0x33, 0x6e, 0x09, 0x87, 0x04, 0x87, 0x86, 0x2c, 0x9f, 0xb9, 0x30, 0x09, 0x72, 0x6c, 0x5e, 0x89, 0x14, 0x48, 0xdf, 0xda, 0xc6, 0xc6, 0xf5, 0x99, 0x8a, 0x11, 0x6b, 0xca, 0x52}},
        {64, 8, {1073892684874442211u, 3564563940764872595u, 7759104081498736782u, 11166360155194816610u, 12346774032542627136u, 15524082577519477263u, 15829369155827348149u, 16724752985191108517u}, {0xc3, 0x70, 0x6e, 0x1b, 0xca, 0xbc, 0x54, 0x81, 0xd6, 0xbc, 0xc2, 0x4e, 0x30, 0x78, 0x4c, 0x0e, 0x99, 0x5d, 0x0f, 0xa6, 0x82, 0x08, 0x83, 0xb3, 0x0e, 0x78, 0x3f, 0x31, 0xf2, 0x4f, 0xef, 0xa9, 0x77, 0xb4, 0x83, 0xb4, 0x91, 0x08, 0xb7, 0xf4, 0xd7, 0x68, 0x0a, 0x1e, 0x3f, 0x1f, 0x9b, 0x6d, 0x6b, 0xb0, 0x64, 0x06, 0xa3, 0xd5, 0x7e, 0x80, 0x0f, 0xb9, 0x0b, 0xa8, 0x38, 0x83, 0x40, 0x0f}},
    };
    for (const auto& vec : vectors) {
        if (!Minisketch::BitsSupported(vec.bits)) continue;
        auto sketches = CreateSketches(vec.bits, vec.capacity);
        for (auto& sketch : sketches) {
            for (uint64_t elem : vec.elements) sketch.Add(elem);
            CHECK(sketch.Serialize() == vec.serialized);
            sketch.Deserialize(vec.serialized);
            std::vector<uint64_t> decoded(vec.capacity);
            CHECK(sketch.Decode(decoded));
            std::sort(decoded.begin(), decoded.end());
            CHECK(decoded == vec.elements);
        }
    }
}

void TestComputeFunctions() {
    for (uint32_t bits = 0; bits <= 256; ++bits) {
        for (uint32_t fpbits = 0; fpbits <= 512; ++fpbits) {
            std::vector<size_t> table_max_elements(1025);
            for (size_t capacity = 0; capacity <= 1024; ++capacity) {
                table_max_elements[capacity] = minisketch_compute_max_elements(bits, capacity, fpbits);
                // Exception for bits==0
                if (bits == 0) CHECK(table_max_elements[capacity] == 0);
                // A sketch with capacity N cannot guarantee decoding more than N elements.
                CHECK(table_max_elements[capacity] <= capacity);
                // When asking for N bits of false positive protection, either no solution exists, or no more than ceil(N / bits) excess capacity should be needed.
                if (bits > 0) CHECK(table_max_elements[capacity] == 0 || capacity - table_max_elements[capacity] <= (fpbits + bits - 1) / bits);
                // Increasing capacity by one, if there is a solution, should always increment the max_elements by at least one as well.
                if (capacity > 0) CHECK(table_max_elements[capacity] == 0 || table_max_elements[capacity] > table_max_elements[capacity - 1]);
            }

            std::vector<size_t> table_capacity(513);
            for (size_t max_elements = 0; max_elements <= 512; ++max_elements) {
                table_capacity[max_elements] = minisketch_compute_capacity(bits, max_elements, fpbits);
                // Exception for bits==0
                if (bits == 0) CHECK(table_capacity[max_elements] == 0);
                // To be able to decode N elements, capacity needs to be at least N.
                if (bits > 0) CHECK(table_capacity[max_elements] >= max_elements);
                // A sketch of N bits in total cannot have more than N bits of false positive protection;
                if (bits > 0) CHECK(bits * table_capacity[max_elements] >= fpbits);
                // When asking for N bits of false positive protection, no more than ceil(N / bits) excess capacity should be needed.
                if (bits > 0) CHECK(table_capacity[max_elements] - max_elements <= (fpbits + bits - 1) / bits);
                // Increasing max_elements by one can only increment the capacity by 0 or 1.
                if (max_elements > 0 && fpbits < 256) CHECK(table_capacity[max_elements] == table_capacity[max_elements - 1] || table_capacity[max_elements] == table_capacity[max_elements - 1] + 1);
                // Check round-tripping max_elements->capacity->max_elements (only a lower bound)
                CHECK(table_capacity[max_elements] <= 1024);
                CHECK(table_max_elements[table_capacity[max_elements]] == 0 || table_max_elements[table_capacity[max_elements]] >= max_elements);
            }

            for (size_t capacity = 0; capacity <= 512; ++capacity) {
                // Check round-tripping capacity->max_elements->capacity (exact, if it exists)
                CHECK(table_max_elements[capacity] <= 512);
                CHECK(table_max_elements[capacity] == 0 || table_capacity[table_max_elements[capacity]] == capacity);
            }
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    uint64_t test_complexity = 4;
    for (int i = 1; i < argc; ++i) {
        std::string arg{argv[i]};
        if (arg.compare(0, 7, "--seed=") == 0) continue; // Handled by GetTestSeed.
        size_t len = 0;
        try {
            test_complexity = 0;
            long long complexity = std::stoll(arg, &len);
            if (complexity >= 1 && len == arg.size() && ((uint64_t)complexity <= std::numeric_limits<uint64_t>::max() >> 10)) {
                test_complexity = complexity;
            }
        } catch (const std::logic_error&) {}
        if (test_complexity == 0) {
            fprintf(stderr, "Invalid complexity specified: '%s'\n", arg.c_str());
            return 1;
        }
    }

    TestRand rng(GetTestSeed(argc, argv));

#ifdef MINISKETCH_VERIFY
    const char* mode = " in verify mode";
#else
    const char* mode = "";
#endif
    printf("Running libminisketch tests%s with complexity=%llu\n", mode, (unsigned long long)test_complexity);

    TestGoldenVectors();
    TestComputeFunctions();

    for (unsigned j = 2; j <= 64; ++j) {
        TestRandomized(j, 8, (test_complexity << 10) / j, rng);
        TestRandomized(j, 128, (test_complexity << 7) / j, rng);
        TestRandomized(j, 4096, test_complexity / j, rng);
    }

    // Test capacity==0 together with all field sizes, and then
    // all combinations of bits and capacity up to a certain bits*capacity,
    // depending on test_complexity.
    for (int weight = 0; weight <= 40; ++weight) {
        for (int bits = 2; weight == 0 ? bits <= 64 : (bits <= 32 && bits <= weight); ++bits) {
            int capacity = weight / bits;
            if (capacity * bits != weight) continue;
            TestExhaustive(bits, capacity);
        }
        if (weight >= 16 && test_complexity >> (weight - 16) == 0) break;
    }

    printf("All tests successful.\n");
    return 0;
}
