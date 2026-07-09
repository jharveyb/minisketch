/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Fuzz target: decode arbitrary serialized sketches.
 *
 * This is the coverage-guided version of TestExhaustive in test.cpp: for an
 * attacker-controlled serialization, all field implementations must agree on
 * decodability and on the decoded set, decoding with a lower limit must fail,
 * and re-encoding the decoded set must reproduce the sketch.
 */

#include <algorithm>
#include <cstdint>
#include <vector>

#include "../../include/minisketch.h"
#include "FuzzedDataProvider.h"
#include "fuzz.h"

namespace {

/** Construct sketches for all supported implementations of a field size. */
std::vector<Minisketch> CreateAll(uint32_t bits, size_t capacity) {
    std::vector<Minisketch> ret;
    for (uint32_t impl = 0; impl <= Minisketch::MaxImplementation(); ++impl) {
        if (Minisketch::ImplementationSupported(bits, impl)) {
            ret.push_back(Minisketch(bits, impl, capacity));
            FUZZ_CHECK((bool)ret.back());
        }
    }
    return ret;
}

} // namespace

FUZZ_TARGET(decode) {
    FuzzedDataProvider provider(data, size);
    uint32_t bits = provider.ConsumeIntegralInRange<uint32_t>(2, 64);
    if (!Minisketch::BitsSupported(bits)) return;
    size_t capacity = provider.ConsumeIntegralInRange<size_t>(1, 24);

    auto sketches = CreateAll(bits, capacity);
    FUZZ_CHECK(!sketches.empty()); // Implementation 0 always exists for supported field sizes.

    // Fix the root-finding basis: freshly constructed sketches seed it from
    // std::random_device, which would make runs non-reproducible.
    uint64_t seed = provider.ConsumeIntegral<uint64_t>();
    for (auto& sketch : sketches) sketch.SetSeed(seed);

    // The remaining bytes are the serialized sketch (zero-padded if short).
    std::vector<unsigned char> serialized = provider.ConsumeBytes<unsigned char>(sketches[0].GetSerializedSize());
    serialized.resize(sketches[0].GetSerializedSize(), 0);

    // Canonicalize: Deserialize ignores padding bits in the last byte, so
    // compare everything against the roundtripped serialization.
    sketches[0].Deserialize(serialized);
    auto canonical = sketches[0].Serialize();

    std::vector<uint64_t> decode0(capacity);
    bool decodable0 = sketches[0].Decode(decode0);
    std::sort(decode0.begin(), decode0.end());

    for (size_t impl = 1; impl < sketches.size(); ++impl) {
        sketches[impl].Deserialize(serialized);
        FUZZ_CHECK(sketches[impl].Serialize() == canonical);
        std::vector<uint64_t> decode_other(capacity);
        bool decodable_other = sketches[impl].Decode(decode_other);
        FUZZ_CHECK(decodable_other == decodable0);
        std::sort(decode_other.begin(), decode_other.end());
        FUZZ_CHECK(decode_other == decode0);
    }

    if (!decodable0) return;

    // Decoded elements are nonzero, in range, and distinct.
    FUZZ_CHECK(decode0.size() <= capacity);
    for (size_t i = 0; i < decode0.size(); ++i) {
        FUZZ_CHECK(decode0[i] != 0);
        FUZZ_CHECK(bits == 64 || (decode0[i] >> bits) == 0);
        if (i > 0) FUZZ_CHECK(decode0[i] != decode0[i - 1]);
    }

    for (auto& sketch : sketches) {
        // Decoding with a limit one below the element count must fail.
        if (!decode0.empty()) {
            std::vector<uint64_t> too_small(decode0.size() - 1);
            FUZZ_CHECK(!sketch.Decode(too_small));
        }
        // Re-encoding the decoded set must reproduce the sketch.
        Minisketch rebuilt(bits, 0, capacity);
        rebuilt.SetSeed(seed);
        for (uint64_t elem : decode0) rebuilt.Add(elem);
        FUZZ_CHECK(rebuilt.Serialize() == canonical);
    }
}
