/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Fuzz target: add a data-driven element multiset to sketches, then check
 * serialization roundtrips, merge additivity, and exact decode when the
 * effective set fits the capacity.
 */

#include <algorithm>
#include <cstdint>
#include <map>
#include <vector>

#include "../../include/minisketch.h"
#include "FuzzedDataProvider.h"
#include "fuzz.h"

FUZZ_TARGET(roundtrip) {
    FuzzedDataProvider provider(data, size);
    uint32_t bits = provider.ConsumeIntegralInRange<uint32_t>(2, 64);
    if (!Minisketch::BitsSupported(bits)) return;
    // Log-uniform capacity in [1, 1024] (see decode.cpp).
    size_t cap_lo = size_t{1} << provider.ConsumeIntegralInRange<int>(0, 10);
    size_t capacity = provider.ConsumeIntegralInRange<size_t>(cap_lo, std::min<size_t>(2 * cap_lo - 1, 1024));
    const uint64_t max_elem = bits == 64 ? ~uint64_t{0} : (uint64_t{1} << bits) - 1;

    std::vector<Minisketch> sketches;
    for (uint32_t impl = 0; impl <= Minisketch::MaxImplementation(); ++impl) {
        if (Minisketch::ImplementationSupported(bits, impl)) {
            sketches.push_back(Minisketch(bits, impl, capacity));
        }
    }
    FUZZ_CHECK(!sketches.empty());

    // Fix the root-finding basis: freshly constructed sketches seed it from
    // std::random_device, which would make runs non-reproducible.
    uint64_t seed = provider.ConsumeIntegral<uint64_t>();
    for (auto& sketch : sketches) sketch.SetSeed(seed);

    // Consume elements (duplicates allowed; pairs cancel out in the sketch).
    std::vector<uint64_t> elements;
    while (elements.size() < capacity + 8 && provider.remaining_bytes() > 0) {
        elements.push_back(provider.ConsumeIntegralInRange<uint64_t>(1, max_elem));
    }
    for (auto& sketch : sketches) {
        for (uint64_t elem : elements) sketch.Add(elem);
    }

    // The effective set is the elements with odd multiplicity.
    std::map<uint64_t, size_t> multiplicity;
    for (uint64_t elem : elements) ++multiplicity[elem];
    std::vector<uint64_t> effective;
    for (const auto& entry : multiplicity) {
        if (entry.second & 1) effective.push_back(entry.first);
    }

    // Serializations must agree across implementations and roundtrip.
    auto serialized = sketches[0].Serialize();
    for (auto& sketch : sketches) {
        FUZZ_CHECK(sketch.Serialize() == serialized);
        sketch.Deserialize(serialized);
        FUZZ_CHECK(sketch.Serialize() == serialized);
    }

    // Merging sketches of any split of the elements must equal the sketch of
    // the whole multiset (syndromes are additive).
    size_t split = provider.ConsumeIntegralInRange<size_t>(0, elements.size());
    Minisketch part_a(bits, 0, capacity), part_b(bits, 0, capacity);
    part_a.SetSeed(seed);
    part_b.SetSeed(seed);
    for (size_t i = 0; i < elements.size(); ++i) {
        (i < split ? part_a : part_b).Add(elements[i]);
    }
    part_a.Merge(part_b);
    FUZZ_CHECK(part_a.Serialize() == serialized);

    std::vector<uint64_t> decoded(capacity);
    bool decodable = sketches[0].Decode(decoded);
    std::sort(decoded.begin(), decoded.end());
    for (size_t impl = 1; impl < sketches.size(); ++impl) {
        std::vector<uint64_t> decoded_other(capacity);
        FUZZ_CHECK(sketches[impl].Decode(decoded_other) == decodable);
        std::sort(decoded_other.begin(), decoded_other.end());
        FUZZ_CHECK(decoded_other == decoded);
    }

    if (effective.size() <= capacity) {
        // Within capacity: decoding must succeed and recover the set exactly.
        FUZZ_CHECK(decodable);
        FUZZ_CHECK(decoded == effective);
        if (!effective.empty()) {
            std::vector<uint64_t> too_small(effective.size() - 1);
            FUZZ_CHECK(!sketches[0].Decode(too_small));
        }
    } else if (decodable) {
        // Over capacity: decoding may fail, or may return some set whose
        // sketch is identical to the input sketch (self-consistency).
        Minisketch rebuilt(bits, 0, capacity);
        rebuilt.SetSeed(seed);
        for (uint64_t elem : decoded) rebuilt.Add(elem);
        FUZZ_CHECK(rebuilt.Serialize() == serialized);
    }
}
